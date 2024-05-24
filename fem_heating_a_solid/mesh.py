import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon


def make_mesh(_width, _length, min_element_size):
    assert min_element_size <= min(_width, _length)
    numberOfNodes_w, numberOfNodes_l = int(_width // min_element_size + 1), int(_length // min_element_size + 1)
    step_w, step_l = _width / (numberOfNodes_w - 1), _length / (numberOfNodes_l - 1)

    nodes = np.empty((numberOfNodes_w, numberOfNodes_l), dtype=tuple)
    frontier = np.zeros(nodes.shape, dtype=bool)

    for i in range(numberOfNodes_w):
        for j in range(numberOfNodes_l):
            nodes[i, j] = (i * step_w, j * step_l)
            if (j == numberOfNodes_l or j == 0) and (i == numberOfNodes_w or i == 0):
                frontier[i, j] = True
    # print(nodes.shape)
    # print(nodes)

    elements = np.empty((2 * (numberOfNodes_l-1) * (numberOfNodes_w-1), 3), dtype=tuple)
    for j in range(numberOfNodes_l - 1):
        for i in range(numberOfNodes_w - 1):
            elements[2 * (i + (numberOfNodes_w - 1) * j), :] = ([(i, j), (i + 1, j), (i, j + 1)])
            elements[2 * (i + (numberOfNodes_w - 1) * j) + 1, :] = ([(i + 1, j), (i + 1, j + 1), (i, j + 1)])
            # print(f"2*({i} + {numberOfNodes_w - 1} * {j}) = {2 * (i + (numberOfNodes_w - 1) * j)}")
    # print(elements)
    domain = np.ones(elements.shape[0])
    # print(elements.shape)
    return {
            "nn": nodes.size, "ne": elements.size,
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

    _w, _l = nodes[-1, -1][0], nodes[-1, -1][1]
    ax.set_xlim([-1e-2, _w+1e-2])
    ax.set_ylim([-1e-2, _l+1e-2])
    ax.set_aspect("equal", "box")
    plt.show()


if __name__ == "__main__":
    m = make_mesh(10, 3, 0.5)
    show_mesh(m)
