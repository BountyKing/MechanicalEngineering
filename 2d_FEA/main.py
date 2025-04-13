import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from mesh import read_mesh

class SolidMechanicsSolver:
    def __init__(self, material):
        self.E = material["E"]
        self.nu = material["nu"]
        self.f = material["f"]
        self.rho = material["rho"]
        self.mesh = None
        self.nn = None
        self.ne = None
        self.nodes = None
        self.elements = None
        self.domain = None
        self.node_groups = {}
        self.results = None
        self.strains = None
        self.stresses = None
        self.von_mises_stress = None
        self.custom_dirichlet_conditions = {}  # key: dof index, value: fixed value
        self.custom_neumann_conditions = []    # list of (group_name, fx, fy)


    def read_mesh(self, filename):
        self.mesh = read_mesh(filename)
        self.nn = self.mesh["nn"]
        self.ne = self.mesh["ne"]
        self.nodes = self.mesh["nodes"]
        self.elements = self.mesh["elements"]
        self.domain = self.mesh["domain"]
        self.node_groups = self.mesh["groups"]

    def D_matrix(self):
        E, nu = self.E, self.nu
        coeff = E / (1 - nu ** 2)
        return coeff * np.array([
            [1, nu, 0],
            [nu, 1, 0],
            [0, 0, (1 - nu) / 2]
        ])

    def B_matrix(self, n1, n2, n3):
        x1, y1 = n1
        x2, y2 = n2
        x3, y3 = n3
        A = 0.5 * ((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1))
        b = np.array([y2 - y3, y3 - y1, y1 - y2]) / (2 * A)
        c = np.array([x3 - x2, x1 - x3, x2 - x1]) / (2 * A)
        B = np.zeros((3, 6))
        for i in range(3):
            B[0, 2*i] = b[i]
            B[1, 2*i+1] = c[i]
            B[2, 2*i] = c[i]
            B[2, 2*i+1] = b[i]
        return B, A

    def apply_displacement(self, group_name, ux=0.0, uy=0.0):
        """
        Apply displacement boundary condition to a node group.
        """
        if group_name not in self.node_groups:
            raise ValueError(f"Group '{group_name}' not found in mesh groups.")

        for i in self.node_groups[group_name]:
            for d, val in zip([0, 1], [ux, uy]):
                idx = 2 * i + d
                self.custom_dirichlet_conditions[idx] = val

    def apply_force(self, group_name, fx=0.0, fy=0.0):
        if group_name not in self.node_groups:
            raise ValueError(f"Group '{group_name}' not found in mesh groups.")
        
        # Append to the list of Neumann boundary conditions
        self.custom_neumann_conditions.append((group_name, fx, fy))


    def do(self):
        dof = 2 * self.nn
        K = np.zeros((dof, dof))
        F = np.zeros(dof)
        D = self.D_matrix()

        for elt in self.elements:
            nodes_coords = [self.nodes[i] for i in elt]
            B, A = self.B_matrix(*nodes_coords)
            Ke = A * B.T @ D @ B

            fe = np.zeros(6)
            for i in range(3):
                fx, fy = self.f(*nodes_coords[i])
                fe[2*i] = fx * A / 3
                fe[2*i + 1] = fy * A / 3

            dof_idx = np.array([[2 * n, 2 * n + 1] for n in elt]).flatten()
            K[np.ix_(dof_idx, dof_idx)] += Ke
            F[dof_idx] += fe
            
        
        # --- Apply Neumann (force) boundary conditions ---
        for group_name, fx_val, fy_val in self.custom_neumann_conditions:
            group_nodes = set(self.node_groups[group_name])
            for elt in self.elements:
                for i in range(3):
                    n0, n1 = elt[i], elt[(i + 1) % 3]
                    if n0 in group_nodes and n1 in group_nodes:
                        x0, y0 = self.nodes[n0]
                        x1, y1 = self.nodes[n1]
                        length = np.linalg.norm([x1 - x0, y1 - y0])

                        fx = fx_val * length / 2
                        fy = fy_val * length / 2

                        F[2 * n0]     += fx
                        F[2 * n0 + 1] += fy
                        F[2 * n1]     += fx
                        F[2 * n1 + 1] += fy


        # --- Apply Dirichlet (displacement) boundary conditions ---
        for dof_idx, value in self.custom_dirichlet_conditions.items():
            K[dof_idx, :] = 0
            K[:, dof_idx] = 0
            K[dof_idx, dof_idx] = 1
            F[dof_idx] = value


        # --- Solve ---
        U = np.linalg.solve(K, F)
        self.results = U

        # --- Compute strains and stresses ---
        strains, stresses, von_mises_stress = [], [], []
        for elt in self.elements:
            nodes_coords = [self.nodes[i] for i in elt]
            B, _ = self.B_matrix(*nodes_coords)
            dof_idx = np.array([[2 * n, 2 * n + 1] for n in elt]).flatten()
            ue = U[dof_idx]
            strain = B @ ue
            stress = D @ strain
            strains.append(strain)
            stresses.append(stress)

            sx, sy, txy = stress
            vm = np.sqrt(0.5 * ((sx - sy)**2 + 3 * txy**2))
            von_mises_stress.append(vm)

        self.strains = np.array(strains)
        self.stresses = np.array(stresses)
        self.von_mises_stress = np.array(von_mises_stress)

    def interpolate_von_mises_to_nodes(self):
        nodal_vm = np.zeros(self.nn)
        counts = np.zeros(self.nn)

        for i, elt in enumerate(self.elements):
            vm = self.von_mises_stress[i]
            for node in elt:
                nodal_vm[node] += vm
                counts[node] += 1

        counts[counts == 0] = 1
        return nodal_vm / counts

    def plot_nodal_stress(self, stress_type="von_mises", alpha=1000):
        ux = self.results[0::2]
        uy = self.results[1::2]

        if stress_type == "von_mises":
            stress_at_nodes = self.interpolate_von_mises_to_nodes()
        else:
            raise NotImplementedError("Only von_mises stress is implemented.")

        fig, ax = plt.subplots(figsize=(10, 6))
        triang = tri.Triangulation(self.nodes[:, 0] + alpha * ux, self.nodes[:, 1] + alpha * uy, self.elements)
        tcf = ax.tricontourf(triang, stress_at_nodes, cmap='RdYlGn_r')
        ax.plot(self.nodes[:, 0], self.nodes[:, 1], 'ko', markersize=2)
        
        cbar = fig.colorbar(tcf, orientation='horizontal')
        cbar.set_label("Von Mises (Pa)", fontsize=12, labelpad=10)
        cbar.ax.xaxis.set_label_position("top")
        ax.set_aspect('equal')
        ax.set_title("Nodal Von Mises Stress")
        plt.show()

    def plot_element_stress(self, alpha=1000):
        ux = self.results[0::2]
        uy = self.results[1::2]

        fig, ax = plt.subplots(figsize=(10, 6))
        triang = tri.Triangulation(self.nodes[:, 0] + alpha * ux, self.nodes[:, 1] + alpha * uy, self.elements)
        tcf = ax.tripcolor(triang, facecolors=self.von_mises_stress, cmap='RdYlGn_r', edgecolors='k')
        cbar = fig.colorbar(tcf, orientation='horizontal')
        cbar.set_label("Von Mises (Pa)", fontsize=12, labelpad=10)
        cbar.ax.xaxis.set_label_position("top")
        ax.set_aspect('equal')
        ax.set_title("Von Mises Stress")
        plt.show()
    
    def plot_displacement(self, alpha=1000):
        ux = self.results[0::2]
        uy = self.results[1::2]
        mag = np.sqrt(ux**2 + uy**2)

        fig, ax = plt.subplots(figsize=(10, 6))
        triang = tri.Triangulation(self.nodes[:, 0] + alpha * ux, self.nodes[:, 1] + alpha * uy, self.elements)
        tcf = ax.tricontourf(triang, mag, cmap='RdYlGn_r')
        ax.plot(self.nodes[:, 0], self.nodes[:, 1], '.', markersize=2)
        cbar = fig.colorbar(tcf, orientation='horizontal')
        cbar.set_label("U sum (m)", fontsize=12, labelpad=10)
        cbar.ax.xaxis.set_label_position("top")
        ax.set_aspect('equal')
        ax.set_title("Displacement Magnitude")
        plt.show()


if __name__ == "__main__":
    material = {
        "E": 210e9,
        "nu": 0.3,
        "rho": 7850,
        "f": lambda x, y: (0, 0)
    }

    Px = 1e9  # Pa
    Py = 0  # Pa

    model = SolidMechanicsSolver(material)
    model.read_mesh("mesh/rectangle.msh")
    model.apply_displacement("N_TXTY", 0, 0)
    model.apply_force("N_F", Px, Py)
    model.do()
    model.plot_displacement(alpha=10)
    model.plot_element_stress(alpha=10)
    model.plot_nodal_stress(alpha=10)
