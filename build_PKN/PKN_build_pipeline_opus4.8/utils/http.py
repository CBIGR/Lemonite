"""
Shared HTTP session with connection-level retry.

This is the lowest layer of the rate-limit handling: a single ``requests.Session``
with a ``urllib3`` ``Retry`` mounted on the adapter. It transparently retries
transient connection resets and rate-limit / server status codes
(429, 500, 502, 503, 504), honouring the ``Retry-After`` header. The higher-level
``@retry_api_call`` decorator (api_retry.py) sits on top for per-database
backoff policy and logging; together they make the monthly API-bound run robust.
"""

import threading

import requests
from requests.adapters import HTTPAdapter

try:  # urllib3 v2 and v1 expose Retry at slightly different paths
    from urllib3.util.retry import Retry
except ImportError:  # pragma: no cover
    from requests.packages.urllib3.util.retry import Retry  # type: ignore

_RETRY_STATUS = (429, 500, 502, 503, 504)
_local = threading.local()


def _build_session(total: int = 5, backoff_factor: float = 1.0) -> requests.Session:
    session = requests.Session()
    retry = Retry(
        total=total,
        connect=total,
        read=total,
        status=total,
        backoff_factor=backoff_factor,
        status_forcelist=_RETRY_STATUS,
        allowed_methods=frozenset(['GET', 'POST']),
        respect_retry_after_header=True,
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry, pool_connections=20, pool_maxsize=20)
    session.mount('https://', adapter)
    session.mount('http://', adapter)
    session.headers.update({'User-Agent': 'LemonIte-PKN-pipeline/1.0'})
    return session


def get_session() -> requests.Session:
    """Return a thread-local session (one per worker thread)."""
    sess = getattr(_local, 'session', None)
    if sess is None:
        sess = _build_session()
        _local.session = sess
    return sess
