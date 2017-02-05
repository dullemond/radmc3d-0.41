"""
This module contains models and functions to manage the library/list of models
"""
from _libfunc import *

__all__ = ["updateModelList", "getModelNames", "getModelDesc"]

try:
    from _modellist import *
    from _modellist import _model_list

except:
    _model_list = None


if (_model_list!=None):
    for i in range(len(_model_list)):
        __all__.append(_model_list[i])

