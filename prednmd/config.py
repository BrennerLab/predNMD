"""
Configuration management for predNMD pipeline
"""

import os
import yaml
from pathlib import Path
from typing import Optional, Dict, Any


class Config:
    """Configuration manager for predNMD pipeline"""
    
    def __init__(self, config_file: Optional[str] = None):
        """
        Initialize configuration
        
        Args:
            config_file: Path to YAML configuration file (optional)
        """
        self.config_file = config_file
        self.config = self._load_default_config()
        
        if config_file and os.path.exists(config_file):
            self._load_from_file(config_file)
    
    def _load_default_config(self) -> Dict[str, Any]:
        """Load default configuration"""
        return {
            'reference': {
                'gtf_file': None,
                'genome_fasta': None,
                'cds_fasta': None,
            },
            'annotation': {
                'gnomad_file': None,
                'phylop_bigwig': None,
                'm6a_file': None,
                'expression_file': None,
            },
            'vep': {
                'vep_path': 'vep',
                'cache_dir': None,
                'assembly': 'GRCh37',
            },
            'model': {
                'model_dir': None,
            },
            'runtime': {
                'threads': 1,
                'keep_intermediate': False,
                'canonical_only': True,
            }
        }
    
    def _load_from_file(self, config_file: str):
        """Load configuration from YAML file"""
        with open(config_file, 'r') as f:
            user_config = yaml.safe_load(f)
        
        # Merge with default config
        self._deep_update(self.config, user_config)
    
    def _deep_update(self, base_dict: dict, update_dict: dict):
        """Recursively update nested dictionaries"""
        for key, value in update_dict.items():
            if key in base_dict and isinstance(base_dict[key], dict) and isinstance(value, dict):
                self._deep_update(base_dict[key], value)
            else:
                base_dict[key] = value
    
    def get(self, *keys, default=None):
        """
        Get configuration value using dot notation
        
        Args:
            *keys: Keys to traverse (e.g., 'reference', 'gtf_file')
            default: Default value if key not found
        
        Returns:
            Configuration value or default
        """
        value = self.config
        for key in keys:
            if isinstance(value, dict) and key in value:
                value = value[key]
            else:
                return default
        return value
    
    def set(self, *keys, value):
        """
        Set configuration value using dot notation
        
        Args:
            *keys: Keys to traverse (e.g., 'reference', 'gtf_file')
            value: Value to set
        """
        config = self.config
        for key in keys[:-1]:
            if key not in config:
                config[key] = {}
            config = config[key]
        config[keys[-1]] = value
    
    def validate(self) -> bool:
        """
        Validate that required configuration is present
        
        Returns:
            True if valid, raises ValueError otherwise
        """
        required_fields = [
            ('reference', 'gtf_file'),
            ('reference', 'genome_fasta'),
            ('reference', 'cds_fasta'),
            ('annotation', 'gnomad_file'),
            ('annotation', 'phylop_bigwig'),
            ('annotation', 'm6a_file'),
            ('annotation', 'expression_file'),
            ('vep', 'cache_dir'),
            ('model', 'model_dir'),
        ]
        
        missing = []
        for keys in required_fields:
            if self.get(*keys) is None:
                missing.append('.'.join(keys))
        
        if missing:
            raise ValueError(f"Missing required configuration fields: {', '.join(missing)}")
        
        return True
    
    def save(self, output_file: str):
        """Save current configuration to YAML file"""
        with open(output_file, 'w') as f:
            yaml.dump(self.config, f, default_flow_style=False)
    
    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any]) -> 'Config':
        """Create Config from dictionary"""
        config = cls()
        config.config = config_dict
        return config
