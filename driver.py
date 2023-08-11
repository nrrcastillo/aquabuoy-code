import os
from typing import List, Tuple, Dict, Any, Union
from dataclasses import dataclass, field
from collections.abc import Callable, Iterable
import timeit

import numpy as np

# TODO setting parameters for irregular waves
@dataclass
class WaveParameters:
    