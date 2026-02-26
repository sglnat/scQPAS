"""
Centralized logging configuration for scQPAS.

Provides a single function to configure logging for the entire application.
Logs are written to both console and a log file by default.
"""

import logging
import sys
from pathlib import Path


def configure_logging(log_level=logging.INFO, log_file='scqpas.log'):
    """
    Configure logging for the application.
    
    Logs are written to both console (stdout) and a file by default.
    
    Parameters
    ----------
    log_level : int, optional
        Logging level (logging.DEBUG, logging.INFO, etc.). Default: logging.INFO
    log_file : str, optional
        Path to log file. Default: 'scqpas.log' in current directory.
        Set to None to disable file logging.
    
    Examples
    --------
    >>> from src.logging_config import configure_logging
    >>> configure_logging(log_level=logging.DEBUG)  # Show debug messages, write to scqpas.log
    >>> configure_logging(log_level=logging.INFO, log_file='custom_output.log')  # Custom log file
    >>> configure_logging(log_level=logging.INFO, log_file=None)  # Console only (no file)
    """
    
    # Remove any existing handlers to avoid duplicates
    root_logger = logging.getLogger()
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
    
    # Define format string
    format_string = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    formatter = logging.Formatter(format_string)
    
    # Console handler (always enabled)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(log_level)
    console_handler.setFormatter(formatter)
    
    # Create handlers list
    handlers = [console_handler]
    
    # File handler (enabled by default)
    if log_file:
        # Create parent directories if needed
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(log_level)
        file_handler.setFormatter(formatter)
        handlers.append(file_handler)
    
    # Configure root logger
    logging.basicConfig(
        level=log_level,
        handlers=handlers,
        force=True  # Force reconfiguration even if already configured
    )
