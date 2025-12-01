"""
TSQVT Gravitational Wave Echoes Analysis Package

This package provides tools for searching gravitational wave echoes
predicted by Twistorial Spectral Quantum Vacuum Theory (TSQVT) in
LIGO/Virgo data.

Modules:
--------
- data_retrieval: Download and preprocess LIGO data
- echo_model: Generate TSQVT echo templates
- matched_filter: Perform matched filter searches
- statistical_tests: Compute significance and Bayes factors
- visualization: Plotting utilities

Author: Mohamed H.M. Makraini
Institution: UNED (Universidad Nacional de Educación a Distancia)
Email: mhamed34@alumno.uned.es
Date: November 2025
"""

__version__ = "0.1.0"
__author__ = "Mohamed H.M. Makraini"
__email__ = "mhamed34@alumno.uned.es"

# Import main classes for convenience
from .echo_model import TSQVTEchoModel
from .data_retrieval import LIGODataRetriever
from .matched_filter import MatchedFilterSearcher

__all__ = [
    'TSQVTEchoModel',
    'LIGODataRetriever',
    'MatchedFilterSearcher',
]
