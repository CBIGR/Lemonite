"""
Application-level API retry + throttling — the "agent that handles rate limits".

Layered design (no LLM, fully deterministic):
  1. utils.http  -> connection-level urllib3 Retry on a shared Session.
  2. this module -> @retry_api_call decorator: per-database exponential backoff
     with jitter, explicit HTTP-429 / Retry-After handling, and detailed logging.
  3. utils.cache -> on-disk per-database checkpoints so an interrupted or
     throttled run resumes without re-querying (see process_single_database).

A function decorated with @retry_api_call returns the sentinel ``'none'`` after
exhausting retries so a single failing metabolite never aborts the whole run; the
metabolite id is logged to api_errors.log for a later --resume pass to retry.
"""

import functools
import logging
import random
import time

from requests.exceptions import Timeout, ConnectionError, RequestException, HTTPError

from config import API_RETRY_CONFIG

logger = logging.getLogger('pkn.api')

_DEFAULT = {'max_retries': 5, 'backoff_factor': 2, 'timeout': 15}


def _sleep_with_jitter(base_wait: float):
    """Sleep base_wait seconds plus up to 25% jitter (avoids thundering herd)."""
    time.sleep(base_wait + random.uniform(0, base_wait * 0.25))


def retry_api_call(db_name='default'):
    """
    Decorate an API function with per-database retry logic.

    Retries on timeouts, connection errors, HTTP 429 (respecting Retry-After) and
    5xx server errors; does not retry on other 4xx client errors. Returns ``'none'``
    when all retries are exhausted.
    """
    config = API_RETRY_CONFIG.get(db_name, _DEFAULT)
    max_retries = config['max_retries']
    backoff_factor = config['backoff_factor']

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            for attempt in range(max_retries):
                last = attempt == max_retries - 1
                try:
                    return func(*args, **kwargs)

                except (Timeout, ConnectionError) as e:
                    wait = backoff_factor ** (attempt + 2)
                    kind = 'timeout' if isinstance(e, Timeout) else 'connection error'
                    logger.warning("[%s] %s in %s (attempt %d/%d): %s. Retrying in ~%ss",
                                   db_name, kind, func.__name__, attempt + 1, max_retries,
                                   str(e)[:100], wait)
                    if not last:
                        _sleep_with_jitter(wait)

                except HTTPError as e:
                    resp = getattr(e, 'response', None)
                    status = resp.status_code if resp is not None else None
                    if status == 429:
                        wait = _retry_after(resp, backoff_factor, attempt)
                        logger.warning("[%s] rate limited (429) in %s (attempt %d/%d). Waiting ~%ss",
                                       db_name, func.__name__, attempt + 1, max_retries, wait)
                        if not last:
                            _sleep_with_jitter(wait)
                    elif status in (500, 502, 503, 504):
                        wait = backoff_factor ** (attempt + 2)
                        logger.warning("[%s] server error %s in %s (attempt %d/%d). Retrying in ~%ss",
                                       db_name, status, func.__name__, attempt + 1, max_retries, wait)
                        if not last:
                            _sleep_with_jitter(wait)
                    else:
                        logger.error("[%s] client error %s in %s: %s",
                                     db_name, status, func.__name__, str(e)[:200])
                        return 'none'

                except RequestException as e:
                    resp = getattr(e, 'response', None)
                    status = getattr(resp, 'status_code', None)
                    if status == 429:
                        wait = _retry_after(resp, backoff_factor, attempt)
                        logger.warning("[%s] rate limited (429) in %s (attempt %d/%d). Waiting ~%ss",
                                       db_name, func.__name__, attempt + 1, max_retries, wait)
                    else:
                        wait = backoff_factor ** attempt
                        logger.warning("[%s] request error in %s (attempt %d/%d): %s. Retrying in ~%ss",
                                       db_name, func.__name__, attempt + 1, max_retries, str(e)[:100], wait)
                    if not last:
                        _sleep_with_jitter(wait)

                except Exception as e:  # noqa: BLE001 - last-resort safety net
                    wait = backoff_factor ** attempt
                    logger.warning("[%s] unexpected error in %s (attempt %d/%d): %s. Retrying in ~%ss",
                                   db_name, func.__name__, attempt + 1, max_retries, str(e)[:100], wait)
                    if not last:
                        _sleep_with_jitter(wait)

            logger.error("[%s] %s failed after %d attempts", db_name, func.__name__, max_retries)
            return 'none'

        return wrapper
    return decorator


def _retry_after(response, backoff_factor, attempt):
    """Compute wait time for a 429, honouring a Retry-After header when present."""
    if response is not None:
        ra = response.headers.get('Retry-After')
        if ra:
            try:
                return int(ra) + 1
            except (ValueError, TypeError):
                pass
    return backoff_factor ** (attempt + 3)


def rate_limit_pause(call_count, db_name):
    """Periodically pause after ``pause_after`` calls to keep under rate limits."""
    config = API_RETRY_CONFIG.get(db_name, {})
    pause_after = config.get('pause_after', 100)
    pause_duration = config.get('pause_duration', 5)
    if call_count and call_count % pause_after == 0:
        logger.info("[%s] pausing %ss after %d calls", db_name, pause_duration, call_count)
        time.sleep(pause_duration)
