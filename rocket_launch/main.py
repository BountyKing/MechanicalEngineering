import numpy as np
import matplotlib.pyplot as plt


# Drag force according to time
def force_drag(index):
    v = np.linalg.norm(np.array([velocity[index - 1, 0], velocity[index - 1, 1]]))
    force_drag_norm = Cd * 1.225 * A * pow(v, 2) / 2
    return force_drag_norm


# Lift momentum approximation for easier calculus
def momentum_lift(index):
    angle_rocket = position[index-1, 2]
    velocity_vec = np.array([velocity[index-1, 0], velocity[index-1, 1]])
    coordinate_vec_unit = np.array([1, 0])
    if np.linalg.norm(velocity_vec) <= 1e-6:
        return 0
    else:
        angle_direction = np.arccos(np.dot(velocity_vec, coordinate_vec_unit) /
                                    (np.linalg.norm(velocity_vec) * np.linalg.norm(coordinate_vec_unit)))
        d_theta = angle_rocket - angle_direction
        return -Kl*d_theta - cl * velocity[index-1, 2]


# Force of propulsion
def force_propulsion(index):
    t = index * time_step
    if t <= 0.05 * tf:
        force_propulsion_norm = force_propulsion_threshold * t / (0.05 * tf) + force_propulsion_threshold
    elif t <= 0.1 * tf:
        force_propulsion_norm = -force_propulsion_threshold * (t - 0.05 * tf) / (
                0.05 * tf) + 2 * force_propulsion_threshold
    elif t <= 0.9 * tf:
        force_propulsion_norm = force_propulsion_threshold
    else:
        force_propulsion_norm = 0
    return force_propulsion_norm


# Rocket mass according to time
def m_rocket(index):
    t = index * time_step
    if t <= 0.9 * tf:
        return m0 - (m0 - mf) / tf * t
    else:
        return mf


def plot_graph(x, y, label):
    plt.plot(x, y, label=label)
    plt.legend()
    plt.show()


def strip_list(lists, index):
    res = lists.copy()
    for i, l in enumerate(lists):
        res[i] = l[:index]
    return res


if __name__ == "__main__":
    print("Simulation started")
    # Variable Initializing :
    m0 = 5  # initial mass (kg)
    mf = 4  # final mass (kg)
    tf = 50  # duration of the simulation (s)
    force_propulsion_threshold = 70  # Threshold of the force of propulsion (N)
    Cd = 0.75  # drag coefficient
    Kl = 4  # spring coefficient to simulate a lifting momentum (N.m.rad^-1)
    cl = 0.3  # damping coefficient

    A = np.pi * pow(10e-3, 2)  # area of reference for the drag force (m^2)
    g = 9.81  # gravity coefficient value (m.s^-2)
    v0x = 0  # initial velocity value according to x-axis (m.s^-1)
    v0y = 0  # initial velocity value according to y-axis (m.s^-1)
    v0z = 0
    x0 = 0  # initial position value according to x-axis (m)
    y0 = 0  # initial position value according to y-axis (m)
    theta0 = np.pi / 4  # initial angle value (rad)

    time_step = 1e-3
    time_stamp = np.arange(0, tf, time_step)

    force_propulsion_list = np.empty(time_stamp.shape)
    force_drag_list = np.empty(time_stamp.shape)
    gamma = np.empty((time_stamp.shape[0], 3))
    gamma[0, :] = 0, 0, 0
    velocity = np.empty((time_stamp.shape[0], 3))
    velocity[0, :] = v0x, v0y, v0z
    position = np.empty((time_stamp.shape[0], 3))
    position[0, :] = x0, y0, theta0

    vx, vy, vz = v0x, v0y, v0z
    x, y, theta = x0, y0, theta0
    i = 1
    while position[i - 1, 1] >= 0 and i < time_stamp.shape[0]:
        t = time_stamp[i]
        force_propulsion_list[i] = force_propulsion(i)
        force_drag_list[i] = force_drag(i)
        gx = ((force_propulsion(i) - force_drag(i)) * np.cos(theta)) / m_rocket(i)
        gy = -g + ((force_propulsion(i) - force_drag(i)) * np.sin(theta)) / m_rocket(i)
        gz = momentum_lift(i) / m_rocket(i)
        vx, vy, vz = vx + time_step * gx, vy + time_step * gy, vz + time_step * gz
        x, y, theta = x + time_step * vx, y + time_step * vy, theta + time_step * vz

        gamma[i, :] = [gx, gy, gz]
        velocity[i, :] = [vx, vy, vz]
        position[i, :] = [x, y, theta]

        if i <= 10:
            print("Gam:", gamma[i, :])
            print("Vel:", velocity[i, :])
            print("Pos:", position[i, :])
        i += 1
    if i >= time_stamp.shape[0]:
        print("Not enough iterations for the rocket to land")

    gamma, velocity, position = strip_list([gamma, velocity, position], i)
    time_stamp, force_drag_list, force_propulsion_list = strip_list([time_stamp, force_drag_list, force_propulsion_list], i)
    # plot_graph(time_stamp, force_propulsion_list, "Force of propulsion (N) over time (s)")
    # plot_graph(time_stamp, force_drag_list, "Force of drag (N) over time (s)")
    # plot_graph(time_stamp, force_lift_list, "Force of lift (N) over time (s)")
    # plot_graph(time_stamp, velocity_y, "Velocity along y-axis")

    plot_graph(position[:, 0], position[:, 1], "Trajectory")
    # plot_graph(time_stamp, position[:, 2]*180/np.pi, "Rotation (deg) according to time (s)")
