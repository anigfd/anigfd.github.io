"""Rayleigh-Benard convection and its Lorenz-63 truncation (pure NumPy).

Boussinesq convection between free-slip, fixed-temperature plates at
z=0,1, periodic in x. Uses vorticity-streamfunction form with
(u,w)=(-psi_z,psi_x), which gives zeta = w_x - u_z = psi_xx+psi_zz =
nabla^2 psi (NOT nabla^2 psi = -zeta -- that sign is inconsistent with
this u,w convention and kills the instability entirely; see NOTATION.md).
"""
import numpy as np

RA_C = 27 * np.pi ** 4 / 4          # critical Rayleigh number, free-slip
K_C = np.pi / np.sqrt(2)            # critical horizontal wavenumber


def neutral_ra(k):
    """Ra_c(k) = (k^2+pi^2)^3 / k^2 -- marginal stability curve; minimized
    at k=K_C, where it equals RA_C."""
    return (k ** 2 + np.pi ** 2) ** 3 / k ** 2


def growth_rate(k, Ra, Pr):
    """Linear growth rate of the free-slip conduction state at wavenumber k."""
    kp2 = k ** 2 + np.pi ** 2
    return Pr / kp2 * (Ra * k ** 2 / kp2 - kp2)


class ChannelGrid:
    """Periodic in x (FFT), homogeneous Dirichlet in z (finite difference,
    interior points only) -- a rigid-lid, fixed-temperature channel.

    Poisson solves use a Thomas (tridiagonal) algorithm per horizontal
    wavenumber, vectorized across all wavenumbers at once, with the
    forward-elimination coefficients precomputed ONCE at construction
    (they depend only on kx and the grid, not on the RHS).
    """

    def __init__(self, nx, nz, Lx):
        self.nx, self.nz, self.Lx = nx, nz, Lx
        self.dx = Lx / nx
        self.dz = 1.0 / (nz + 1)
        self.x = np.arange(nx) * self.dx
        self.z = np.arange(1, nz + 1) * self.dz         # interior points only
        self.kx = 2 * np.pi * np.fft.rfftfreq(nx, d=self.dx)
        self.nkx = len(self.kx)

        # Tridiagonal (d^2/dz^2 - kx^2) with Dirichlet BCs (ghosts=0): always
        # invertible for every kx (including kx=0), no gauge fix needed.
        a = c = 1.0 / self.dz ** 2
        b = -2.0 / self.dz ** 2 - self.kx ** 2            # shape (nkx,)
        piv = np.zeros((self.nkx, nz))
        cp = np.zeros((self.nkx, nz))
        piv[:, 0] = b
        cp[:, 0] = c / piv[:, 0]
        for i in range(1, nz):
            piv[:, i] = b - a * cp[:, i - 1]
            if i < nz - 1:
                cp[:, i] = c / piv[:, i]
        self._piv, self._cp, self._a = piv, cp, a

    def poisson_solve(self, rhs):
        """Solve (d^2/dz^2 - kx^2) f_hat = rhs_hat for each kx; rhs and the
        return value are physical-space (nx,nz) fields."""
        rhs_hat = np.fft.rfft(rhs, axis=0)
        dp = np.zeros_like(rhs_hat)
        dp[:, 0] = rhs_hat[:, 0] / self._piv[:, 0]
        for i in range(1, self.nz):
            dp[:, i] = (rhs_hat[:, i] - self._a * dp[:, i - 1]) / self._piv[:, i]
        out_hat = np.zeros_like(rhs_hat)
        out_hat[:, -1] = dp[:, -1]
        for i in range(self.nz - 2, -1, -1):
            out_hat[:, i] = dp[:, i] - self._cp[:, i] * out_hat[:, i + 1]
        return np.fft.irfft(out_hat, n=self.nx, axis=0)

    def ddx(self, f):
        fh = np.fft.rfft(f, axis=0)
        return np.fft.irfft(1j * self.kx[:, None] * fh, n=self.nx, axis=0)

    def ddz(self, f):
        """Central difference, with zero ghost values at z=0,1."""
        out = np.zeros_like(f)
        out[:, 1:-1] = (f[:, 2:] - f[:, :-2]) / (2 * self.dz)
        out[:, 0] = f[:, 1] / (2 * self.dz)
        out[:, -1] = -f[:, -2] / (2 * self.dz)
        return out

    def laplacian(self, f):
        fh = np.fft.rfft(f, axis=0)
        d2x = np.fft.irfft(-self.kx[:, None] ** 2 * fh, n=self.nx, axis=0)
        d2z = np.zeros_like(f)
        d2z[:, 1:-1] = (f[:, :-2] - 2 * f[:, 1:-1] + f[:, 2:]) / self.dz ** 2
        d2z[:, 0] = (-2 * f[:, 0] + f[:, 1]) / self.dz ** 2
        d2z[:, -1] = (f[:, -2] - 2 * f[:, -1]) / self.dz ** 2
        return d2x + d2z

    def jacobian(self, psi, f):
        """J(psi,f) = psi_x f_z - psi_z f_x."""
        return self.ddx(psi) * self.ddz(f) - self.ddz(psi) * self.ddx(f)


def rhs_boussinesq(state, grid, Ra, Pr):
    """RHS of d(state)/dt for state=(zeta,theta) stacked as one (2,nx,nz)
    array:

        zeta_t  = -J(psi,zeta)  + Pr*lap(zeta)  + Pr*Ra*theta_x
        theta_t = -J(psi,theta) +    lap(theta) +       psi_x
        nabla^2 psi = zeta
    """
    zeta, theta = state
    psi = grid.poisson_solve(zeta)
    dzeta = -grid.jacobian(psi, zeta) + Pr * grid.laplacian(zeta) + Pr * Ra * grid.ddx(theta)
    dtheta = -grid.jacobian(psi, theta) + grid.laplacian(theta) + grid.ddx(psi)
    return np.stack([dzeta, dtheta])


def lorenz_rhs(state, sigma, r, b):
    """The Lorenz (1963) truncation: Xdot=sigma(Y-X), Ydot=rX-Y-XZ, Zdot=XY-bZ."""
    x, y, z = state
    return np.array([sigma * (y - x), r * x - y - x * z, x * y - b * z])
