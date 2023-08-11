import logging
import os
import csv
from typing import List, Dict, Tuple, Callable, TypeVar, Generic

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

# equivalent to matlab ode15s in scipy (for stiff equations)
## sp.integrate.ode(f).set_integrator('vode', method='bdf', order=15)
# equivalent to matlab ode45 in scipy (for non-stiff equations)
## sp.integrate.RK45(f, t0, y0, t_bound, max_step=0.1)

eps = TypeVar("eps", List[List[int]])


def plot_data(ys: np.array, save: bool = False) -> None:
    plt.plot(ys)
    if save:
        pass


def figs_to_eps(
    fig: plt.figimage,  # output of matplotlib plot
    path: str,
    fname: str = function.__name__,
) -> eps:
    pass


def summary_to_csv(
    df: pd.DataFrame,  # potentially change this to dict
    path: str,
    filename: str = "output",
) -> csv:
    """
    _summary_

    Args:
        df (pd.DataFrame): _description_
        path (str): _description_
        filename (str, optional): _description_. Defaults to "output".

    Returns:
        csv: on specified path, function will output dataframe as "filename.csv"
        TODO: find other way to default name

        csv output format:
        - numerical scheme used
        - duration of operation--tentatively use time.timeit()
        - accuracy (?)
    """

    pass


def plot_functions(**kwargs):
    ax, fig = plt.subplot()
    for k, v in kwargs.items():
        ax.plot(v)
    plt.savefig("output.eps", format="eps")
    plt.show()
    return fig
