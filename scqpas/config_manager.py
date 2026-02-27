"""
Configuration manager for scQPAS.

Loads and manages configuration from YAML files with fallback to defaults.
"""

import logging
from pathlib import Path
from typing import Any, Optional

try:
    import yaml
except ImportError:
    raise ImportError(
        "PyYAML is required for configuration management. "
        "Install it with: pip install pyyaml"
    )

logger = logging.getLogger(__name__)


class ConfigManager:
    """Load and manage scQPAS configuration from YAML files."""

    def __init__(self, config_path: Optional[str] = None):
        """
        Initialize config manager.

        Parameters
        ----------
        config_path : str, optional
            Path to custom config file. If None, uses default config.
        """
        if config_path:
            self.config_path = Path(config_path)
            if not self.config_path.exists():
                raise FileNotFoundError(f"Config file not found: {config_path}")
            logger.info(f"Loaded custom config from: {config_path}")
        else:
            self.config_path = self._get_default_config_path()
            logger.debug(f"Using default config from: {self.config_path}")

        self.config = self._load_config()

    @staticmethod
    def _get_default_config_path() -> Path:
        """Get path to default config file."""
        return Path(__file__).parent / "config" / "defaults.yaml"

    def _load_config(self) -> dict:
        """
        Load YAML configuration file.

        Returns
        -------
        dict
            Configuration dictionary
        """
        try:
            with open(self.config_path) as f:
                config = yaml.safe_load(f)
                if config is None:
                    config = {}
                return config
        except yaml.YAMLError as e:
            raise ValueError(f"Invalid YAML configuration file: {e}")

    def get(self, section: str, key: str, default: Any = None) -> Any:
        """
        Get configuration value with fallback to default.

        Parameters
        ----------
        section : str
            Configuration section (e.g., 'polya_detection')
        key : str
            Configuration key (e.g., 'percentage_threshold')
        default : Any, optional
            Fallback default value if key not found

        Returns
        -------
        Any
            Configuration value or default

        Examples
        --------
        >>> config = ConfigManager()
        >>> pct_threshold = config.get('polya_detection', 'percentage_threshold')
        >>> 80
        """
        try:
            return self.config.get(section, {}).get(key, default)
        except (AttributeError, TypeError):
            return default

    def get_section(self, section: str) -> dict:
        """
        Get entire configuration section.

        Parameters
        ----------
        section : str
            Configuration section name

        Returns
        -------
        dict
            Configuration section dictionary

        Examples
        --------
        >>> config = ConfigManager()
        >>> polya_config = config.get_section('polya_detection')
        """
        return self.config.get(section, {})

    def to_dict(self) -> dict:
        """Return full configuration as dictionary."""
        return self.config.copy()

    def __repr__(self) -> str:
        return f"ConfigManager(config_path={self.config_path})"
