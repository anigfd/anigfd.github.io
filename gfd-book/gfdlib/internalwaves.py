"""Internal gravity waves in a stratified, rotating fluid (pure NumPy).

The (x,z) vertical plane is solved with the SAME doubly-periodic spectral
machinery as the horizontal chapters (`gfdlib.spectral.Grid`), with the
grid's y-axis relabeled as the vertical z. Streamfunction convention here is
u=psi_z, w=-psi_x (2-D incompressible flow in a vertical plane) -- see
NOTATION.md. This is NOT the same sign convention as the horizontal (x,y)
streamfunction used in the barotropic/QG chapters; the two live in different
physical planes.
"""
import numpy as np


def dispersion_omega(kx, kz, N, f):
    """omega(k,m) = sqrt((N^2 k^2 + f^2 m^2)/(k^2+m^2)); bounded in [f, N]."""
    k2 = kx ** 2 + kz ** 2
    return np.sqrt((N ** 2 * kx ** 2 + f ** 2 * kz ** 2) / k2)


def beam_angle(omega, N, f):
    """Angle (rad) of the wavevector from vertical -- equivalently the energy
    beam's angle from horizontal:

        cos(theta) = sqrt((omega^2 - f^2) / (N^2 - f^2)).
    """
    cos_theta = np.sqrt(np.clip((omega ** 2 - f ** 2) / (N ** 2 - f ** 2), 0.0, 1.0))
    return np.arccos(cos_theta)


def step_leapfrog(q, q_prev, grid, N2, f, source, dt):
    """One Stormer-Verlet step of

        d^2q/dt^2 = -(N2(z)*psi_xx + f^2*psi_zz) + source,   q = nabla^2 psi

    N2 may vary with z (the grid's y-axis) -- pass a field of grid.x.shape
    (built from grid.y only, for a stratified but x-independent N^2(z)).
    Reduces to the exact per-mode SHO d^2q/dt^2 = -omega_k^2 q when N2 is
    constant (see `dispersion_omega`), since psi_xx=-kx^2 psi and
    psi_zz=-kz^2 psi then combine through the same k^2=kx^2+kz^2 that
    defines q=nabla^2 psi.
    """
    psi_hat = grid.invert_laplacian(grid.fft(q))
    psi_xx = grid.ifft(-grid.kx ** 2 * psi_hat)
    psi_zz = grid.ifft(-grid.ky ** 2 * psi_hat)
    accel = -(N2 * psi_xx + f ** 2 * psi_zz) + source
    return 2 * q - q_prev + dt ** 2 * accel
