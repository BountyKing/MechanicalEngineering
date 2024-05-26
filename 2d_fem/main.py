# based on https://perso.univ-lyon1.fr/marc.buffat/COURS/BOOK_ELTFINIS_HTML/CoursEF/chap4.html#equation-eq4-21

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri
from mesh import make_mesh, read_mesh, show_mesh
import utilities as util
import os


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

    def make_mesh(self, corners, borders, p):
        self.mesh = make_mesh(corners, borders, p)
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
        # Constructing matrix A and B of the equation A X = B
        A = np.zeros((self.nn, self.nn))
        B = np.zeros(self.nn)

        for elt in self.elements:
            n1, n2, n3 = self.nodes[elt[0]], self.nodes[elt[1]], self.nodes[elt[2]]
            x1, y1 = n1
            x2, y2 = n2
            x3, y3 = n3

            # Calculating elementary matrix Ke
            Area = 0.5 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))
            p = np.zeros(5, dtype=int)
            p[0:3] = elt[:]
            p[3:5] = elt[0:2]
            dX = self.nodes[np.ix_(p[1:4], [0, 1])] - self.nodes[np.ix_(p[2:5], [0, 1])]

            dN = np.array([dX[:, 1], -dX[:, 0]]) / (dX[0, 0] * dX[1, 1] - dX[1, 0] * dX[0, 1])

            K22 = np.dot(dN[:, 1], dN[:, 1]) * Area
            K33 = np.dot(dN[:, 2], dN[:, 2]) * Area
            K23 = np.dot(dN[:, 1], dN[:, 2]) * Area

            Ke = np.array([[K22 + K33 + 2 * K23, -K23 - K22, -K23 - K33],
                           [-K23 - K22, K22, K23],
                           [-K23 - K33, K23, K33]])

            # Calculating elementary matrix Me
            Me = Area / 12.0 * np.array([[2, 1, 1],
                                         [1, 2, 1],
                                         [1, 1, 2]])

            # Retrieve matrix Fe, matrix given by the f function
            Fe = util.get_Fe([n1, n2, n3], self.f)

            # c.f. course cited
            Be = np.dot(Me, Fe)

            # Adding the elementary matrixes to the matrices of the main equation
            A[np.ix_(elt, elt)] += self.k * Ke + self.alpha * Me
            B[np.ix_(elt)] += Be

        # Get the edges affected by BC
        edges_bc = []
        for k in range(self.ne):
            edges = [*self.elements[k, :], self.elements[k, 0]]
            for j in range(3):
                if not self.frontier[edges[j]] == 0 and not self.frontier[edges[j + 1]] == 0:
                    edges_bc.append([edges[j], edges[j + 1]])

        # Remove the edges counted twice
        edges_bc_wo_double = []
        for ij in edges_bc:
            i, j = ij[0], ij[1]
            if [i, j] not in edges_bc_wo_double and [j, i] not in edges_bc_wo_double:
                edges_bc_wo_double.append([i, j])

        edges_bc = edges_bc_wo_double
        for edge in edges_bc:
            # Boundary driven by a Fourier type of condition
            if self.frontier[edge[0]] == 4 and self.frontier[edge[1]] == 4:
                dl = np.linalg.norm(self.nodes[edge[0]] - self.nodes[edge[1]])
                A[edge[0], edge[0]] += self.beta * dl / 3
                A[edge[0], edge[1]] += self.beta * dl / 6
                A[edge[1], edge[0]] += self.beta * dl / 3
                A[edge[1], edge[1]] += self.beta * dl / 6
                B[edge[0]] -= self.phi0 * dl / 2
                B[edge[1]] -= self.phi0 * dl / 2

        for i in range(self.nn):
            # Homogeneous Dirichlet BC
            if self.frontier[i] == 1:
                A[i, :] = 0
                A[:, i] = 0
                A[i, i] = 1.0
                B[i] = 0
            # Non-homogeneous Dirichlet BC
            elif self.frontier[i] == 2:
                A[i, :] = 0
                # A[:, i] = 0  # Error in the course cited previously
                A[i, i] = 1.0
                B[i] = self.ue

        # Solving the equation A X = B
        U = np.linalg.solve(A, B)
        self.results = U

    def show_results(self):
        # print(self.results)
        fig, ax = plt.subplots()
        fig.set_size_inches(18.5, 10.5)
        triangles = matplotlib.tri.Triangulation(self.nodes[:, 0], self.nodes[:, 1], self.elements)
        contour = ax.tricontourf(triangles, self.results, cmap="hot_r")
        ax.tricontour(triangles, self.results, colors='k', linewidths=1)
        fig.colorbar(contour)
        ax.set_aspect("equal", "box")
        ax.set_title(f"Results")
        plt.tight_layout()
        plt.show()

    def save_results(self, filename=None):
        if not filename:
            file_saving_index = str(len(os.listdir("./output/"))//2).rjust(4, "0")
            filename = f"output-{file_saving_index}"

        with open("output/" + filename + ".csv", "w") as f:
            f.write("//nodes matrix\n")
            for i in range(self.nodes.shape[0]):
                f.write(f"{self.nodes[i, 0]},{self.nodes[i, 1]}\n")
            f.write("//result matrix\n")
            for i in range(self.results.shape[0]):
                f.write(f"{self.results[i]}\n")

        fig, ax = plt.subplots()
        triangles = matplotlib.tri.Triangulation(self.nodes[:, 0], self.nodes[:, 1], self.elements)
        contour = ax.tricontourf(triangles, self.results)
        ax.tricontour(triangles, self.results, colors='k', linewidths=1)
        fig.colorbar(contour)
        plt.savefig(f"output/" + filename + ".pdf")


if __name__ == "__main__":
    constants = {
        "ue": 2.0,
        "alpha": 0.0,
        "beta": 0.0,
        "phi0": 0.0,
        "k": 1.0,
        "f": lambda x, y: 0.0
    }

    plate = Study(constants)
    plate.read_mesh("mesh/maillage.msh")
    plate.show_mesh()
    plate.do()
    plate.show_results()
