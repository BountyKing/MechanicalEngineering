import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon


# TO REWRITE !!!
def make_mesh(length, height, min_element_size):
    assert min_element_size <= min(length, height)
    numberOfNodes_l, numberOfNodes_h = int(length // min_element_size + 1), int(height // min_element_size + 1)
    step_w, step_l = length / (numberOfNodes_l - 1), height / (numberOfNodes_h - 1)

    nodes = np.empty((numberOfNodes_l * numberOfNodes_h), dtype=tuple)
    frontier = np.zeros(nodes.shape, dtype=bool)

    for i in range(numberOfNodes_l):
        for j in range(numberOfNodes_h):
            nodes[i * numberOfNodes_h + j] = (i * step_w, j * step_l)
            if (j == numberOfNodes_h or j == 0) and (i == numberOfNodes_l or i == 0):
                frontier[i, j] = True
    print(nodes.shape)
    # print(nodes)

    elements = np.empty((2 * (numberOfNodes_h - 1) * (numberOfNodes_l - 1), 3), dtype=tuple)
    for j in range(numberOfNodes_h - 1):
        for i in range(numberOfNodes_l - 1):
            elements[2 * (i + (numberOfNodes_l - 1) * j), :] = ([(i, j), (i + 1, j), (i, j + 1)])
            elements[2 * (i + (numberOfNodes_l - 1) * j) + 1, :] = ([(i + 1, j), (i + 1, j + 1), (i, j + 1)])
            # print(f"2*({i} + {numberOfNodes_l - 1} * {j}) = {2 * (i + (numberOfNodes_l - 1) * j)}")
    # print(elements)

    domain = np.ones(elements.shape[0])
    print(elements.shape)
    return {
        "nn": nodes.shape[0], "ne": elements.size,
        "nodes": nodes, "elements": elements,
        "frontier": frontier, "domain": domain
    }


def read_mesh(filename):
    with open(filename, "r") as f:
        line = f.readline().split()
        nn, ne = int(line[0]), int(line[1])
        nodes = np.empty((nn, 2), dtype=float)
        frontier = np.empty(nn, dtype=int)
        for i in range(nn):
            line = f.readline().split()
            nodes[i, 0] = float(line[0])
            nodes[i, 1] = float(line[1])
            frontier[i] = int(line[2])

        elements = np.empty((ne, 3), dtype=int)
        domain = np.empty(ne, dtype=int)
        for j in range(ne):
            line = f.readline().split()
            elements[j, :] = [int(vertex) - 1 for vertex in line[:-1]]
            domain[j] = int(line[-1])

    return {
        "nn": nn, "ne": ne,
        "nodes": nodes, "elements": elements,
        "frontier": frontier, "domain": domain
    }


def show_mesh(_mesh):
    fig, ax = plt.subplots()
    nodes = _mesh["nodes"]
    for elt in _mesh["elements"]:
        n1, n2, n3 = nodes[elt[0]], nodes[elt[1]], nodes[elt[2]]
        # print("elt:", elt)
        p = Polygon([n1, n2, n3], ec="black")
        ax.add_patch(p)

    for i, n in enumerate(nodes):
        ax.text(n[0], n[1], str(i+1))

    _w, _l = max(nodes[:, 0]), max(nodes[:, 1])
    min_w, min_l = min(nodes[:, 0]), min(nodes[:, 1])
    ax.set_xlim([min_w - 1e-2, _w + 1e-2])
    ax.set_ylim([min_l - 1e-2, _l + 1e-2])
    ax.set_aspect("equal", "box")
    plt.show()


if __name__ == "__main__":
    # m = make_mesh(10, 3, 0.5)
    m = read_mesh("maillage.msh")
    show_mesh(m)
