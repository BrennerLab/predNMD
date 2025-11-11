"""
predNMD: A pipeline for predicting Nonsense-Mediated Decay (NMD) in genetic variants

This package provides a comprehensive pipeline for:
- Filtering VCF files for protein-coding regions
- VEP annotation
- Feature extraction and annotation
- NMD prediction using Random Forest models
"""

from .version import __version__
from .pipeline import NMDPipeline
from .config import Config

__all__ = ['NMDPipeline', 'Config', '__version__']
