import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon
from scipy.spatial import Delaunay


def make_mesh(_corners, _borders, p):
    if len(_corners) == 4:
        assert len(_borders) == 4  # assert that we have the 4 frontiers (even if one is empty)

        x0, y0 = _corners[0, 0], _corners[0, 1]
        x1, y1 = _corners[1, 0], _corners[1, 1]
        x2, y2 = _corners[2, 0], _corners[2, 1]
        x3, y3 = _corners[3, 0], _corners[3, 1]

        l1 = np.linalg.norm(_corners[0, :] - _corners[1, :])
        l2 = np.linalg.norm(_corners[1, :] - _corners[2, :])
        l3 = np.linalg.norm(_corners[2, :] - _corners[3, :])
        l4 = np.linalg.norm(_corners[3, :] - _corners[0, :])

        number_elements_L = int(max(l1, l3) // p)
        number_elements_H = int(max(l2, l4) // p)
        # print(number_elements_L, number_elements_H)

        line1 = np.empty((number_elements_L + 1, 2))
        line3 = np.empty((number_elements_L + 1, 2))

        for i in range(number_elements_L + 1):
            line1[i, :] = [x0 + i * (x1 - x0) / number_elements_L, y0 + i * (y1 - y0) / number_elements_L]
            line3[i, :] = [x3 + i * (x2 - x3) / number_elements_L, y3 + i * (y2 - y3) / number_elements_L]

        points = np.empty(((number_elements_L + 1) * (number_elements_H + 1), 2))
        frontier = np.empty(((number_elements_L + 1) * (number_elements_H + 1)), dtype=int)
        for i in range(number_elements_H + 1):
            for j in range(number_elements_L + 1):
                x0, y0 = line1[j, 0], line1[j, 1]
                x1, y1 = line3[j, 0], line3[j, 1]
                x = x0 + i * (x1 - x0) / number_elements_H
                y = y0 + i * (y1 - y0) / number_elements_H
                points[i * (number_elements_L + 1) + j, :] = [x, y]
                frontier[i * (number_elements_L + 1) + j] = is_part_of([x, y], _borders, _corners)

        tri = Delaunay(points)

        nodes = points
        nn = nodes.shape[0]
        elements = tri.simplices
        ne = len(elements)
        domain = np.ones(elements.shape[0], dtype=int)  # Only one-solid-studies

        return {
            "nn": nn, "ne": ne,
            "nodes": nodes, "elements": elements,
            "frontier": frontier, "domain": domain
        }
    elif len(_corners) == 6:
        # print(_borders, _corners)
        p0 = _corners[0]
        p1 = _corners[1]
        p2 = [_corners[2, 0], _corners[4, 1]]
        p3 = _corners[5]
        borders_sub = [_borders[i].copy() for i in range(4)]
        for i in range(4):
            j = 0
            while j < len(borders_sub[i]) // 2:
                b = borders_sub[i][2 * j:2 * j + 2]
                if (2 in b and 3 in b) or (3 in b and 4 in b):
                    borders_sub[i] = borders_sub[i][:2 * j] + borders_sub[i][2 * j + 2:]
                    j -= 1
                elif 4 in b and 5 in b:
                    borders_sub[i][2 * j:2 * j + 2] = [2, 3]
                elif 5 in b and 0 in b:
                    borders_sub[i][2 * j:2 * j + 2] = [3, 0]
                j += 1
        print("b_sub_: ", borders_sub)
        mesh_1 = make_mesh(np.array([p0, p1, p2, p3]), borders_sub, p)
        # show_mesh(mesh_1)
        p0 = _corners[2]
        p1 = _corners[3]
        p2 = _corners[4]
        p3 = [_corners[2, 0], _corners[4, 1]]
        borders_sub_2 = [_borders[i].copy() for i in range(4)]
        for i in range(4):
            j = 0
            while j < len(_borders[i]) // 2:
                b = borders_sub_2[i][2 * j:2 * j + 2]
                if (0 in b and 1 in b) or (1 in b and 2 in b) or (5 in b and 0 in b):
                    borders_sub_2[i] = borders_sub_2[i][:2 * j] + borders_sub_2[i][2 * j + 2:]
                    j -= 1
                elif 2 in b and 3 in b:
                    borders_sub_2[i][2 * j:2 * j + 2] = [0, 1]
                elif 3 in b and 4 in b:
                    borders_sub_2[i][2 * j:2 * j + 2] = [1, 2]
                elif 4 in b and 5 in b:
                    borders_sub_2[i][2 * j:2 * j + 2] = [2, 3]
                j += 1
        print(borders_sub_2)
        mesh_2 = make_mesh(np.array([p0, p1, p2, p3]), borders_sub_2, p)
        # show_mesh(mesh_2)

        mesh = mesh_1.copy()
        index_mesh_2_1 = {}

        for i, node in enumerate(mesh_2["nodes"]):
            index = index_of_this_node(node, mesh_1["nodes"])
            if not index:  # if the node isn't already in mesh_1, we append it to mesh
                mesh["nodes"] = np.append(mesh["nodes"], [node], axis=0)
                mesh["frontier"] = np.append(mesh["frontier"], [mesh_2["frontier"][i]], axis=0)
                index_mesh_2_1[i] = mesh["nodes"].shape[0] - 1
            else:  # else we update the dictionary to switch indexes in elements of mesh_2
                # print("duplicate node", node)
                index_mesh_2_1[i] = index
                mesh["frontier"][index] = is_part_of(mesh["nodes"][index], _borders, _corners)

        # print(index_mesh_2_1)
        for i in range(mesh_2["elements"].shape[0]):
            elt = mesh_2["elements"][i]
            new_elt = elt.copy()
            # print(elt, end=' ')
            for k in range(len(elt)):
                new_elt[k] = index_mesh_2_1[elt[k]]
            # print(new_elt)
            mesh["elements"] = np.append(mesh["elements"], [new_elt], axis=0)
            mesh["domain"] = np.append(mesh["domain"], [1], axis=0)
        mesh["nn"] = mesh["nodes"].shape[0]
        mesh["ne"] = mesh["elements"].shape[0]
        return mesh
    else:
        raise Exception("Not implemented yet")


def index_of_this_node(n, nodes):
    for i in range(nodes.shape[0]):
        if abs(nodes[i, 0] - n[0]) <= 1e-2 and abs(nodes[i, 1] - n[1]) <= 1e-2:
            return i


def is_part_of(_point, _frontier, _corners):
    k = 0
    x, y = _point[0], _point[1]
    print("corners: ", _corners)
    print("frontier: ", _frontier)
    while k < len(_frontier):
        for n in range(len(_frontier[k]) // 2):
            print(n, _frontier[k][2 * n:2 * n + 2])
            line = np.array([_corners[_frontier[k][2 * n]], _corners[_frontier[k][2 * n + 1]]])
            x1, y1 = line[0][0], line[0][1]
            x2, y2 = line[1][0], line[1][1]
            is_part = (abs((x2 - x1) * (y - y1) - (x - x1) * (y2 - y1)) <= 1e-6
                       and ((x1 <= x <= x2) or (x2 <= x <= x1))
                       and ((y1 <= y <= y2) or (y2 <= y <= y1)))
            # print(is_part)
            if is_part:
                return k + 1
        k += 1

    return 0


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


def save_mesh(_mesh, filename):
    if ".msh" in filename:
        filename = filename[:-4]
    if "mesh/" in filename:
        filename = filename[5:]

    with open("mesh/" + filename + ".msh", "w") as f:
        f.write(f"{_mesh['nn']} {_mesh['ne']} {filename}\n")
        for i in range(_mesh["nn"]):
            f.write(f"{_mesh['nodes'][i, 0]} {_mesh['nodes'][i, 1]} {_mesh['frontier'][i]}\n")
        for j in range(_mesh["ne"]):
            f.write(f"{_mesh['elements'][j, 0] + 1} {_mesh['elements'][j, 1] + 1} {_mesh['elements'][j, 2] + 1}"
                    f" {_mesh['domain'][j]}\n")


def show_mesh(_mesh, show_number=False):
    fig, ax = plt.subplots()
    fig.set_size_inches(18.5, 10.5)
    nodes = _mesh["nodes"]
    for elt in _mesh["elements"]:
        n1, n2, n3 = nodes[elt[0]], nodes[elt[1]], nodes[elt[2]]
        # print("elt:", elt)
        p = Polygon([n1, n2, n3], ec="black", fc="none", linewidth=0.5)
        ax.add_patch(p)

    label_showed = [0]
    for i, n in enumerate(nodes):
        if show_number:
            ax.text(n[0], n[1], str(i + 1))
        f = _mesh["frontier"][i]
        if f == 0:
            c = "black"
        elif f == 1:
            c = "red"
            label = "$\Gamma_1$"
        elif f == 2:
            c = 'green'
            label = "$\Gamma_2$"
        elif f == 3:
            c = "b"
            label = "$\Gamma_3$"
        else:
            c = "cyan"
            label = "$\Gamma_4$"

        if f not in label_showed:
            ax.plot(n[0], n[1], color=c, marker="o", label=label, markersize=5)
            label_showed.append(f)
        else:
            ax.plot(n[0], n[1], color=c, marker="o", markersize=5)

    ax.set_aspect("equal", "box")
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    ax.set_title(f"Mesh with {_mesh['nn']} nodes and {_mesh['ne']} elements")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # Build the mesh for case 3 and 4
    H = 1.
    U0 = 1 / 2
    corners = np.array([[0, 0], [4*H, 0], [4*H, H], [5*H, H], [5*H, 2*H], [0, 2*H]])
    BC = [[3, 4], [5, 0], [0, 1, 1, 2, 2, 3, 4, 5], []]

    m = make_mesh(corners, BC, 125e-3)
    # print(m["frontier"])
    save_mesh(m, "potential_flow_even_finer.msh")
    show_mesh(m)


    """# Build the mesh for case 1 and 2 (only precision changes)
    rectangle = np.array([[0, 0], [5, 0], [5, 2], [0, 2]], dtype=float)
    BC = [[], [0, 3], [], [0, 1, 1, 2, 2, 3]]
    #     G_1  G_2    G_3        G_4
    m = make_mesh(rectangle, BC, 5e-1)
    save_mesh(m, "thin_plate.msh")
    show_mesh(m)"""
