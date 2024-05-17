import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon


def make_mesh(_width, _length, precision):
    assert precision <= min(_width, _length)
    numberOfNodes_w, numberOfNodes_l = _width // precision + 1, _length // precision + 1
    step_w, step_l = _width / (numberOfNodes_w - 1), _length / (numberOfNodes_l - 1)

    nodes = np.empty((numberOfNodes_w, numberOfNodes_l), dtype=tuple)

    for i in range(numberOfNodes_w):
        for j in range(numberOfNodes_l):
            nodes[i, j] = (i * step_w, j * step_l)
    print(nodes.shape)
    print(nodes)

    elements = []
    for j in range(numberOfNodes_l - 1):
        for i in range(numberOfNodes_w - 1):
            elements.append([(i, j), (i+1, j), (i, j+1)])
            elements.append([(i+1, j), (i+1, j+1), (i, j+1)])
    print(elements)
    # print(elements.shape)
    return {"nodes": nodes, "elements": elements}


def show_mesh(_mesh):
    fig, ax = plt.subplots()
    nodes = _mesh["nodes"]
    for elt in _mesh["elements"]:
        n1, n2, n3 = nodes[elt[0]], nodes[elt[1]], nodes[elt[2]]
        p = Polygon([n1, n2, n3], fc="black", ec="red")
        ax.add_patch(p)

    _w, _l = nodes[-1, -1][0], nodes[-1, -1][1]
    ax.set_xlim([-1e-2, _w+1e-2])
    ax.set_ylim([-1e-2, _l+1e-2])
    plt.show()


if __name__ == "__main__":
    m = make_mesh(3, 3, 1)
    show_mesh(m)
