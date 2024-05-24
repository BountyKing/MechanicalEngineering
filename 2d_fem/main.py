# https://stackoverflow.com/questions/52202014/how-can-i-plot-2d-fem-results-using-matplotlib
# https://scicomp.stackexchange.com/questions/31463/efficiently-plot-a-finite-element-mesh-solution-with-matplotlib
# ctrl alt shift l to reformate the file

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri
from mesh import make_mesh, read_mesh, show_mesh
import utilities as util
import scipy.integrate as integrate


class Study:
    def __init__(self, _constants):
        self.k = _constants["k"]
        self.alpha = _constants["alpha"]
        self.ue = _constants["ue"]
        self.phi0 = _constants["phi0"]
        self.beta = _constants["beta"]
        self.f = _constants["f"]
        self.mesh = None
        self.nn = None
        self.ne = None
        self.nodes = None
        self.elements = None
        self.frontier = None
        self.domain = None
        self.results = None

    def make_mesh(self, _geometry, p):
        self.mesh = make_mesh(_geometry["length"], _geometry["height"], p)
        self.nn = self.mesh["nn"]
        self.ne = self.mesh["ne"]
        self.nodes = self.mesh["nodes"]
        self.elements = self.mesh["elements"]
        self.frontier = self.mesh["frontier"]
        self.domain = self.mesh["domain"]

    def read_mesh(self, filename):
        self.mesh = read_mesh(filename)
        self.nn = self.mesh["nn"]
        self.ne = self.mesh["ne"]
        self.nodes = self.mesh["nodes"]
        self.elements = self.mesh["elements"]
        self.frontier = self.mesh["frontier"]
        self.domain = self.mesh["domain"]

    def show_mesh(self):
        show_mesh(self.mesh)

    def do(self):
        A = np.zeros((self.nn, self.nn))
        B = np.zeros(self.nn)

        for elt in self.elements:
            n1, n2, n3 = self.nodes[elt[0]], self.nodes[elt[1]], self.nodes[elt[2]]
            x1, y1 = n1
            x2, y2 = n2
            x3, y3 = n3

            Area = 0.5 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))
            p = np.zeros(5, dtype=int)
            p[0:3] = elt[:]
            p[3:5] = elt[0:2]
            dN = self.nodes[np.ix_(p[1:4], [0, 1])] - self.nodes[np.ix_(p[2:5], [0, 1])]
            dN = np.array([dN[:, 1], -dN[:, 0]])
            K22 = np.dot(dN[:, 1], dN[:, 1]) / Area
            K23 = np.dot(dN[:, 1], dN[:, 2]) / Area
            K33 = np.dot(dN[:, 2], dN[:, 2]) / Area
            Ke = np.array([[K22 + K33 + 2 * K23, -K23 - K22, -K23 - K33],
                           [-K23 - K22, K22, K23],
                           [-K23 - K33, K23, K33]])

            Me = Area / 12.0 * np.array([[2, 1, 1],
                                         [1, 2, 1],
                                         [1, 1, 2]])

            Fe = util.get_Fe([n1, n2, n3], self.f)

            Be = np.dot(Me, Fe)
            A[np.ix_(elt, elt)] += self.k * Ke
            B[np.ix_(elt)] += Be

        # Boundary conditions
        edges_bc = []
        for k in range(self.ne):
            edges = [*self.elements[k, :], self.elements[k, 0]]
            for j in range(3):
                if not self.frontier[edges[j]] == 0 and not self.frontier[edges[j + 1]] == 0:
                    edges_bc.append([edges[j], edges[j + 1]])
        edges_bc_wo_double = []
        for ij in edges_bc:
            i, j = ij[0], ij[1]
            if [i, j] not in edges_bc_wo_double and [j, i] not in edges_bc_wo_double:
                edges_bc_wo_double.append([i, j])

        edges_bc = edges_bc_wo_double
        for edge in edges_bc:
            if self.frontier[edge[0]] == 4 and self.frontier[edge[1]] == 4:
                dl = np.linalg.norm(self.nodes[edge[0]] - self.nodes[edge[1]])
                A[edge[0], edge[0]] += self.beta * dl / 3
                A[edge[0], edge[1]] += self.beta * dl / 6
                A[edge[1], edge[0]] += self.beta * dl / 3
                A[edge[1], edge[1]] += self.beta * dl / 6
                B[edge[0]] -= self.phi0 * dl / 2
                B[edge[1]] -= self.phi0 * dl / 2
        for i in range(self.nn):
            if self.frontier[i] == 1:
                A[i, :] = 0
                A[:, i] = 0
                A[i, i] = 1.0
                B[i] = 0
            elif self.frontier[i] == 2:
                A[i, :] = 0
                A[:, i] = 0
                A[i, i] = 1.0
                B[i] = self.ue
        #   Dirichlet
        # A[np.where(self.frontier == 1), :] = 0
        # A[np.where(self.frontier == 1), np.where(self.frontier == 1)] = 1
        # B[np.where(self.frontier == 1)] = 0
        # A[np.where(self.frontier == 2), :] = 0
        # A[np.where(self.frontier == 2), np.where(self.frontier == 1)] = 1
        # B[np.where(self.frontier == 2)] = self.ue
        np.savetxt("A.log.txt", A, delimiter=",")
        U = np.linalg.solve(A, B)
        self.results = U

    def show_results(self):
        print(self.results)
        fig, ax = plt.subplots()
        triangles = matplotlib.tri.Triangulation(self.nodes[:, 0], self.nodes[:, 1], self.elements)
        contour = ax.tricontourf(triangles, self.results)
        ax.tricontour(triangles, self.results, colors='k')
        fig.colorbar(contour)
        plt.show()


if __name__ == "__main__":
    # c = 5  # speed of the air surrounding the solid
    # H = 10.45 - c + 10 * np.sqrt(c)
    # EPSILON = 1
    # SIGMA = 5.670374419e-8

    # geometry = {
    #     "length": 20,
    #     "height": 10,
    #     "thickness": 1
    # }

    constants = {
        "ue": 2,
        "alpha": 0.0,
        "beta": 0.0,
        "phi0": 0.0,
        "k": 1.0,
        "f": lambda x, y: 0
    }

    plate = Study(constants)
    plate.read_mesh("maillage.msh")
    plate.show_mesh()
    plate.do()
    plate.show_results()
