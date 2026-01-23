"""
API retry decorator and utilities for handling flaky API calls.

Provides intelligent retry logic with exponential backoff for database APIs
that may experience timeouts or connection issues.
"""

import functools
import logging
import time
from requests.exceptions import Timeout, ConnectionError, RequestException
from config import API_RETRY_CONFIG


def retry_api_call(db_name='default'):
    """
    Decorator that adds intelligent retry logic to any API function.
    
    Features:
    - Exponential backoff for failed requests
    - Extra delays for connection timeouts
    - Database-specific retry configurations
    - Detailed logging of retry attempts
    
    Parameters:
    -----------
    db_name : str
        Name of database to get specific retry config
        Options: 'UniProtKB', 'IntAct', 'chEMBL', 'ChEMBL_Mapping', 'STRING', 'default'
    
    Returns:
    --------
    Decorated function with retry logic
    
    Example:
    --------
    >>> @retry_api_call(db_name='UniProtKB')
    ... def fetch_uniprot_data(accession):
    ...     response = requests.get(f"https://rest.uniprot.org/uniprotkb/{accession}")
    ...     return response.json()
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
                    # Connection issues - use extra exponential backoff
                    wait_time = backoff_factor ** (attempt + 2)
                    error_type = "Connection timeout" if isinstance(e, Timeout) else "Connection error"
                    logging.warning(
                        f"[{db_name}] {error_type} in {func.__name__} "
                        f"(attempt {attempt + 1}/{max_retries}): {str(e)[:100]}. "
                        f"Retrying in {wait_time}s..."
                    )
                    if attempt < max_retries - 1:
                        time.sleep(wait_time)
                
                except RequestException as e:
                    # Other request errors - standard backoff
                    wait_time = backoff_factor ** attempt
                    logging.warning(
                        f"[{db_name}] Request error in {func.__name__} "
                        f"(attempt {attempt + 1}/{max_retries}): {str(e)[:100]}. "
                        f"Retrying in {wait_time}s..."
                    )
                    if attempt < max_retries - 1:
                        time.sleep(wait_time)
                
                except Exception as e:
                    # Unexpected errors - standard backoff
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
            return 'none'  # Return 'none' for failed metabolite lookups
        
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
    
    Example:
    --------
    >>> for i, metabolite in enumerate(metabolites):
    ...     data = fetch_api_data(metabolite)
    ...     rate_limit_pause(i + 1, 'UniProtKB')
    """
    config = API_RETRY_CONFIG.get(db_name, {})
    pause_after = config.get('pause_after', 100)
    pause_duration = config.get('pause_duration', 5)
    
    if call_count % pause_after == 0:
        logging.info(f"[{db_name}] Pausing for {pause_duration}s after {call_count} API calls...")
        time.sleep(pause_duration)
