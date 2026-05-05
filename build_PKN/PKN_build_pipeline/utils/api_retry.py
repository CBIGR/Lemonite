"""
API retry decorator and utilities for handling flaky API calls.

Provides intelligent retry logic with exponential backoff for database APIs
that may experience timeouts, connection issues, or rate limiting.
"""

import functools
import logging
import time
from requests.exceptions import Timeout, ConnectionError, RequestException, HTTPError
from config import API_RETRY_CONFIG


def retry_api_call(db_name='default'):
    """
    Decorator that adds intelligent retry logic to any API function.
    
    Features:
    - Exponential backoff for failed requests
    - Extra delays for connection timeouts
    - HTTP 429 (Too Many Requests) handling with Retry-After header support
    - Database-specific retry configurations
    - Detailed logging of retry attempts
    
    Parameters:
    -----------
    db_name : str
        Name of database to get specific retry config
        Options: 'UniProtKB', 'IntAct', 'chEMBL', 'ChEMBL_Mapping', 'STRING', 'MyGene', 'default'
    
    Returns:
    --------
    Decorated function with retry logic
    """
    config = API_RETRY_CONFIG.get(db_name, {
        'max_retries': 5,
        'backoff_factor': 2,
        'timeout': 15
    })
    
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            max_retries = config['max_retries']
            backoff_factor = config['backoff_factor']
            
            for attempt in range(max_retries):
                try:
                    return func(*args, **kwargs)
                
                except (Timeout, ConnectionError) as e:
                    wait_time = backoff_factor ** (attempt + 2)
                    error_type = "Connection timeout" if isinstance(e, Timeout) else "Connection error"
                    logging.warning(
                        f"[{db_name}] {error_type} in {func.__name__} "
                        f"(attempt {attempt + 1}/{max_retries}): {str(e)[:100]}. "
                        f"Retrying in {wait_time}s..."
                    )
                    if attempt < max_retries - 1:
                        time.sleep(wait_time)
                
                except HTTPError as e:
                    response = getattr(e, 'response', None)
                    status_code = response.status_code if response is not None else None
                    
                    if status_code == 429:
                        # Rate limited - respect Retry-After header if present
                        retry_after = None
                        if response is not None:
                            retry_after = response.headers.get('Retry-After')
                        if retry_after:
                            try:
                                wait_time = int(retry_after) + 1
                            except ValueError:
                                wait_time = backoff_factor ** (attempt + 3)
                        else:
                            wait_time = backoff_factor ** (attempt + 3)
                        
                        logging.warning(
                            f"[{db_name}] Rate limited (HTTP 429) in {func.__name__} "
                            f"(attempt {attempt + 1}/{max_retries}). "
                            f"Waiting {wait_time}s before retry..."
                        )
                        if attempt < max_retries - 1:
                            time.sleep(wait_time)
                    elif status_code in (500, 502, 503, 504):
                        # Server errors - retry with backoff
                        wait_time = backoff_factor ** (attempt + 2)
                        logging.warning(
                            f"[{db_name}] Server error (HTTP {status_code}) in {func.__name__} "
                            f"(attempt {attempt + 1}/{max_retries}). "
                            f"Retrying in {wait_time}s..."
                        )
                        if attempt < max_retries - 1:
                            time.sleep(wait_time)
                    else:
                        # Client error (4xx other than 429) - don't retry
                        logging.error(
                            f"[{db_name}] Client error (HTTP {status_code}) in {func.__name__}: "
                            f"{str(e)[:200]}"
                        )
                        return None
                
                except RequestException as e:
                    # Check for 429 in general RequestException too
                    response = getattr(e, 'response', None)
                    status_code = getattr(response, 'status_code', None) if response is not None else None
                    
                    if status_code == 429:
                        retry_after = response.headers.get('Retry-After') if response is not None else None
                        if retry_after:
                            try:
                                wait_time = int(retry_after) + 1
                            except ValueError:
                                wait_time = backoff_factor ** (attempt + 3)
                        else:
                            wait_time = backoff_factor ** (attempt + 3)
                        logging.warning(
                            f"[{db_name}] Rate limited (HTTP 429) in {func.__name__} "
                            f"(attempt {attempt + 1}/{max_retries}). "
                            f"Waiting {wait_time}s..."
                        )
                    else:
                        wait_time = backoff_factor ** attempt
                        logging.warning(
                            f"[{db_name}] Request error in {func.__name__} "
                            f"(attempt {attempt + 1}/{max_retries}): {str(e)[:100]}. "
                            f"Retrying in {wait_time}s..."
                        )
                    if attempt < max_retries - 1:
                        time.sleep(wait_time)
                
                except Exception as e:
                    wait_time = backoff_factor ** attempt
                    logging.warning(
                        f"[{db_name}] Unexpected error in {func.__name__} "
                        f"(attempt {attempt + 1}/{max_retries}): {str(e)[:100]}. "
                        f"Retrying in {wait_time}s..."
                    )
                    if attempt < max_retries - 1:
                        time.sleep(wait_time)
            
            # All retries exhausted
            logging.error(f"[{db_name}] Failed to execute {func.__name__} after {max_retries} attempts")
            return None
        
        return wrapper
    return decorator


def rate_limit_pause(call_count, db_name):
    """
    Add periodic pauses to avoid overwhelming APIs.
    
    Parameters:
    -----------
    call_count : int
        Current number of API calls made
    db_name : str
        Database name to check pause configuration
    """
    config = API_RETRY_CONFIG.get(db_name, {})
    pause_after = config.get('pause_after', 100)
    pause_duration = config.get('pause_duration', 5)
    
    if call_count % pause_after == 0:
        logging.info(f"[{db_name}] Pausing for {pause_duration}s after {call_count} API calls...")
        time.sleep(pause_duration)
