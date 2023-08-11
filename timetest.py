import cProfile
import numpy as np


def ft():
    t = 3
    dt = 0.5
    test = np.array([t - n * dt for n in range(1, 1000)])
    return test


cProfile.run("ft()")
