## 2D FEM project

This project is based on a course given by a professor at the University of Lyon 1, available at :
[Mr Buffat's course](https://perso.univ-lyon1.fr/marc.buffat/COURS/BOOK_ELTFINIS_HTML/CoursEF/chap4.html#conditions-aux-limites)

It is adapted to solid mechanics. The project is divided in two main components: `mesh.py` and `main.py`.

### `mesh.py`

This program build a mesh from the given inputs:

- Coordinates of a convex quadrilateral (polygon)
- Group edges : this will build a group of node to further apply BC or load.
- A mesh precision (spacing)

### `main.py`

This is the actual solver. To solve a problem, you'll need:

- Mesh : `read_mesh`
- Boundary conditions : `apply_displacement`
- Load : `apply_force`

Then you can compute the solution with `do()`. And finally plot the results with `plot_displacement`, ` plot_elelment_stress`, `plot_nodal_stress`.
