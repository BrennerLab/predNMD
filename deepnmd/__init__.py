"""
DeepNMD: A pipeline for predicting Nonsense-Mediated Decay (NMD) in genetic variants

This package provides a comprehensive pipeline for:
- Filtering VCF files for protein-coding regions
- VEP annotation
- Feature extraction and annotation
- NMD prediction using Random Forest models
"""

__version__ = "1.0.0"
__author__ = "Your Name"
__email__ = "your.email@example.com"

from .pipeline import NMDPipeline
from .config import Config

__all__ = ['NMDPipeline', 'Config']
