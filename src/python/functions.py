"""Functions used in this pipeline."""

import pandas as pd
import numpy as np


def bool_to_int(x):
    x = str(x)
    valids = [
    ["Y"     ,   "N"],
    ["YES"   ,   "NO"],
    ["T"     ,   "F"],
    ["TRUE"  ,   "FALSE"],
    ['1'     ,   '0'],
    ]

    true, false = zip(*valids)

    if x.upper() in true:
        return 1
    elif x.upper() in false:
        return 0
    else:
        raise ValueError("`{x}` {x_type} not recognized as 'True' OR 'False'.".format(x_type=type(x),x=x))
