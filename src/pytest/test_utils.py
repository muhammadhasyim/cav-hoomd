"""Test utilities for CavityMD test suite."""

import pytest
import signal
from functools import wraps


def timeout(seconds=5):
    """
    Decorator to add a timeout to test functions.
    
    Parameters
    ----------
    seconds : int
        Timeout duration in seconds (default: 5)
    
    Examples
    --------
    @timeout(5)
    def test_long_simulation():
        # Test code that should complete in 5 seconds
        pass
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            def timeout_handler(signum, frame):
                raise TimeoutError(f"Test {func.__name__} exceeded {seconds} second timeout")
            
            # Set the signal alarm
            old_handler = signal.signal(signal.SIGALRM, timeout_handler)
            signal.alarm(seconds)
            
            try:
                result = func(*args, **kwargs)
            finally:
                # Reset the alarm
                signal.alarm(0)
                signal.signal(signal.SIGALRM, old_handler)
            
            return result
        
        return wrapper
    return decorator


# Pytest marker for quick tests with timeout
quick_test = pytest.mark.timeout(5)

