"""Conservation / spectrum diagnostics (pure NumPy)."""
import numpy as np


def isotropic_spectrum(F, grid):
    """Shell-averaged spectrum E(k) from a spectral field F (rfft2 layout).

    Returns (k_centers, E_k) with sum(E_k) ~ total spectral variance so it can
    be overlaid against k^-5/3 / k^-3 reference slopes.
    """
    n = grid.n
    power = (np.abs(F) ** 2)
    # account for the halved real-FFT spectrum (double interior columns)
    power[:, 1:-1] *= 2
    kmag = grid.kmag.ravel()
    p = power.ravel()
    kbins = np.arange(0.5, kmag.max() + 1.0, 1.0)
    idx = np.digitize(kmag, kbins)
    E = np.zeros(len(kbins) - 1)
    for i in range(1, len(kbins)):
        E[i - 1] = p[idx == i].sum()
    kc = 0.5 * (kbins[:-1] + kbins[1:])
    return kc, E / (n ** 4)


def energy_enstrophy(psi_hat, grid):
    """From spectral streamfunction: kinetic energy and enstrophy (2D flow)."""
    u2 = grid.k2 * np.abs(psi_hat) ** 2
    z2 = grid.k2 ** 2 * np.abs(psi_hat) ** 2
    norm = grid.n ** 4
    # double interior rfft columns
    u2[:, 1:-1] *= 2; z2[:, 1:-1] *= 2
    return 0.5 * u2.sum() / norm, 0.5 * z2.sum() / norm
