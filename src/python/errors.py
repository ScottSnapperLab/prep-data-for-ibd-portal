#!/usr/bin/env python
"""Provide custom errors."""

# Imports


# Metadata
__author__ = "Gus Dunn"
__email__ = "w.gus.dunn@gmail.com"


class PipelineError(Exception):

    """Base class for local errors."""



class ValidationError(PipelineError):

    """Raise this when a validation/sanity check fails."""
