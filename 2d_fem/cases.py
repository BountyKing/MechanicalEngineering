from main import Study
import sys

if __name__ == "__main__":
    if len(sys.argv) == 2:
        if sys.argv[1] == '1':
            # CASE 1 : Blade heated from the left side
            constants = {
                "ue": 200.0,
                "alpha": 0.0,
                "beta": 0.0,
                "phi0": 5.0,
                "k": 1.0,
                "f": lambda x, y: 0.0
            }

            case1 = Study(constants)
            case1.read_mesh("mesh/thin_plate.msh")
            case1.show_mesh()
            case1.do()
            case1.show_results()
        elif sys.argv[1] == '2':
            # CASE 2 : Blade heated from the left side with a finer meshing
            constants = {
                "ue": 200.0,
                "alpha": 0.0,
                "beta": 0.0,
                "phi0": 5.0,
                "k": 1.0,
                "f": lambda x, y: 0.0
            }

            case2 = Study(constants)
            case2.read_mesh("mesh/thin_plate_finer.msh")
            case2.show_mesh()
            case2.do()
            case2.show_results()
        elif sys.argv[1] == '3':
            # CASE 3 : potential flow
            H = 1.
            U0 = 1 / 2

            constants = {
                "ue": 5 * U0 * H,
                "alpha": 0.0,
                "beta": 0.0,
                "phi0": 5.0,
                "k": 1.0,
                "f": lambda x, y: 0.0
            }

            case3 = Study(constants)
            case3.read_mesh("mesh/potential_flow.msh")
            case3.show_mesh()
            case3.do()
            case3.show_results()
        elif sys.argv[1] == '4':
            # CASE 4 : potential flow with a finer mesh
            H = 1.
            U0 = 1 / 2

            constants = {
                "ue": 5 * U0 * H,
                "alpha": 0.0,
                "beta": 0.0,
                "phi0": 5.0,
                "k": 1.0,
                "f": lambda x, y: 0.0
            }

            case4 = Study(constants)
            case4.read_mesh("mesh/potential_flow_even_finer.msh")
            case4.show_mesh()
            case4.do()
            case4.show_results()
        else:
            raise Exception("Not implemented yet")
