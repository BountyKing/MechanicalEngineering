# https://stackoverflow.com/questions/52202014/how-can-i-plot-2d-fem-results-using-matplotlib
# https://scicomp.stackexchange.com/questions/31463/efficiently-plot-a-finite-element-mesh-solution-with-matplotlib



import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from mesh import make_mesh


class Study:
    def __init__(self, w, l, temperatures):
        self.W = w
        self.L = l
        self.T0 = temperatures[0]
        self.TL = temperatures[1]
        self.TT = temperatures[2]
        self.TR = temperatures[3]
        self.TB = temperatures[4]
        self.mesh = None

    def make_mesh(self, p):
        self.mesh = make_mesh(self.W, self.L, p)


if __name__ == "__main__":
    L = 10  # length of the solid (m)
    W = 10  # width of the solid (m)
    T0 = 200  # initial temperature (K)
    TL = 500  # temperature of the left side (K)
    TT, TR, TB = 300, 300, 300  # temperature of the Top, Right and Bottom side (K)

    plate = Study(L, W, )
