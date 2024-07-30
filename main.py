import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

# These masses represent the Earth-Moon system
m_1 = 5.974e24  # kg
m_2 = 7.348e22  # kg
pi_2 = m_2 / (m_1 + m_2)

def nondim_cr3bp(t, Y):
    """Solve the CR3BP in nondimensional coordinates.

    The state vector is Y, with the first three components as the
    position of $m$, and the second three components its velocity.

    The solution is parameterized on $\\pi_2$, the mass ratio.
    """
    # Get the position and velocity from the solution vector
    x, y, z = Y[:3]
    xdot, ydot, zdot = Y[3:]

    # Define the derivative vector
    Ydot = np.zeros_like(Y)
    Ydot[:3] = Y[3:]

    sigma = np.sqrt(np.sum(np.square([x + pi_2, y, z])))
    psi = np.sqrt(np.sum(np.square([x - 1 + pi_2, y, z])))
    Ydot[3] = (
        2 * ydot
        + x
        - (1 - pi_2) * (x + pi_2) / sigma**3
        - pi_2 * (x - 1 + pi_2) / psi**3
    )
    Ydot[4] = -2 * xdot + y - (1 - pi_2) * y / sigma**3 - pi_2 * y / psi**3
    Ydot[5] = -(1 - pi_2) / sigma**3 * z - pi_2 / psi**3 * z
    return Ydot

if __name__ == "__main__":
    ang = 122.7 / 180 * np.pi
    x_0 = 0.04 * np.cos(ang)
    y_0 = 0.04 * np.sin(ang)
    z_0 = 0
    speed = 8.9044
    vx_0 = np.cos(ang - np.pi/2) * speed
    vy_0 = np.sin(ang - np.pi/2) * speed
    vz_0 = 0

    # Then stack everything together into the state vector
    r_0 = np.array((x_0, y_0, z_0))
    v_0 = np.array((vx_0, vy_0, vz_0))
    Y_0 = np.hstack((r_0, v_0))

    t_0 = 0  # nondimensional time
    t_f = 1.0  # nondimensional time
    t_points = np.linspace(t_0, t_f, 10000000)
    sol = solve_ivp(nondim_cr3bp, [t_0, t_f], Y_0, t_eval=t_points, atol=1e-9, rtol=1e-6)

    Y = sol.y.T
    r = Y[:, :3]  # nondimensional distance
    v = Y[:, 3:]  # nondimensional velocity

    # Plot the speed magnitude
    speed = np.linalg.norm(v, axis=1)

    x_2 = (1 - pi_2) * np.cos(np.linspace(0, np.pi, 100))
    y_2 = (1 - pi_2) * np.sin(np.linspace(0, np.pi, 100))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5, 5), dpi=96)

    # Plot the orbits
    ax1.plot(r[:, 0], r[:, 1], "r", label="Trajectory")
    ax1.axhline(0, color="k")
    ax1.plot(np.hstack((x_2, x_2[::-1])), np.hstack((y_2, -y_2[::-1])))
    ax1.plot(-pi_2, 0, "bo", label="$m_1$")
    ax1.plot(1 - pi_2, 0, "go", label="$m_2$")
    ax1.plot(x_0, y_0, "ro")
    ax1.set_aspect("equal")
    
    ax2.plot(sol.t, speed)
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Speed")
    ax2.set_aspect("equal")
    plt.show()