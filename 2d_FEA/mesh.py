import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon


def generate_points_in_polygon(polygon, spacing):
    n = len(polygon)
    if n != 4:
        raise NotImplementedError("Only convex quadrilaterals supported.")

    a, b, c, d = polygon
    length = int(np.linalg.norm(a - b) / spacing)
    height = int(np.linalg.norm(a - d) / spacing)

    points = []
    for i in range(height + 1):
        s = i / height
        start = (1 - s) * a + s * d
        end = (1 - s) * b + s * c
        for j in range(length + 1):
            t = j / length
            pt = (1 - t) * start + t * end
            points.append(pt)

    return np.array(points)


def is_on_edge(pt, p1, p2, tol=1e-6):
    x, y = pt
    x1, y1 = p1
    x2, y2 = p2

    cross = abs((x2 - x1) * (y - y1) - (x - x1) * (y2 - y1))
    dot = (x - x1) * (x2 - x1) + (y - y1) * (y2 - y1)
    length_sq = (x2 - x1) ** 2 + (y2 - y1) ** 2

    return cross < tol and 0 <= dot <= length_sq


def assign_groups(points, edge_groups):
    """Assign named groups to nodes based on edges."""
    groups = {}
    for name, edge_list in edge_groups.items():
        groups[name] = set()
        for p1, p2 in edge_list:
            for i, pt in enumerate(points):
                if is_on_edge(pt, p1, p2):
                    groups[name].add(i)
        groups[name] = sorted(groups[name])
    return groups


def make_mesh(polygon, group_edges, spacing):
    points = generate_points_in_polygon(polygon, spacing)
    tri = Delaunay(points)
    elements = tri.simplices

    # Build group edge dictionary
    edge_groups = {}
    for group_name, edge_indices in group_edges.items():
        edge_groups[group_name] = [(polygon[e[0]], polygon[e[1]]) for e in edge_indices]

    groups = assign_groups(points, edge_groups)

    return {
        "nodes": points,
        "elements": elements,
        "groups": groups,
        "nn": len(points),
        "ne": len(elements),
        "domain": np.ones(len(elements), dtype=int)
    }


def save_mesh(mesh, filename):
    if ".msh" in filename:
        filename = filename[:-4]
    if "mesh/" in filename:
        filename = filename[5:]

    with open("mesh/" + filename + ".msh", "w") as f:
        f.write(f"{mesh['nn']} {mesh['ne']} {filename}\n")
        for i in range(mesh["nn"]):
            x, y = mesh["nodes"][i]
            f.write(f"{x} {y}\n")
        for elt in mesh["elements"]:
            f.write(f"{elt[0] + 1} {elt[1] + 1} {elt[2] + 1} 1\n")

        # Save groups
        f.write("GROUPS\n")
        for name, indices in mesh["groups"].items():
            f.write(f"{name} {' '.join(map(str, indices))}\n")


def read_mesh(filename):
    with open(filename, "r") as f:
        header = f.readline().split()
        nn, ne = int(header[0]), int(header[1])
        nodes = np.zeros((nn, 2))
        for i in range(nn):
            nodes[i] = list(map(float, f.readline().split()))

        elements = np.zeros((ne, 3), dtype=int)
        for j in range(ne):
            line = f.readline().split()
            elements[j] = [int(x) - 1 for x in line[:3]]

        groups = {}
        if f.readline().strip() == "GROUPS":
            for line in f:
                tokens = line.strip().split()
                name = tokens[0]
                indices = list(map(int, tokens[1:]))
                groups[name] = indices

    return {
        "nodes": nodes,
        "elements": elements,
        "groups": groups,
        "nn": len(nodes),
        "ne": len(elements),
        "domain": np.ones(ne, dtype=int)
    }


def show_mesh(mesh, show_number=False):
    fig, ax = plt.subplots()
    fig.set_size_inches(18.5, 10.5)
    nodes = mesh["nodes"]
    for elt in mesh["elements"]:
        p = Polygon(nodes[elt], ec="black", fc="none", linewidth=0.5)
        ax.add_patch(p)

    colors = ['red', 'green', 'blue', 'cyan', 'magenta', 'orange']
    for i, (name, indices) in enumerate(mesh["groups"].items()):
        for idx in indices:
            x, y = nodes[idx]
            ax.plot(x, y, marker='o', color=colors[i % len(colors)], markersize=5, label=name if idx == indices[0] else "")

    if show_number:
        for i, (x, y) in enumerate(nodes):
            ax.text(x, y, str(i), fontsize=8)

    ax.set_aspect("equal", "box")
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    ax.set_title(f"Mesh with {mesh['nn']} nodes and {mesh['ne']} elements")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    rectangle = np.array([[0, 0], [5, 0], [5, 2], [0, 2]], dtype=float)

    # Define named boundary edges (by index of rectangle corners)
    group_edges = {
        "N_TXTY": [[0, 3]],
        "N_F": [[2, 1]]
    }

    mesh = make_mesh(rectangle, group_edges, spacing=0.1)
    save_mesh(mesh, "mesh/rectangle.msh")
    show_mesh(mesh, show_number=True)
