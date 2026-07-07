"""Pseudo-spectral helpers on a doubly-periodic square domain (pure NumPy)."""
import numpy as np


class Grid:
    """Doubly-periodic [0,L)x[0,L) grid with cached real-FFT wavenumbers.

    Build ONCE per notebook, reuse every timestep. Uses rfft2/irfft2, so the
    spectral field has shape (n, n//2 + 1).
    """

    def __init__(self, n, L=2 * np.pi):
        self.n, self.L = n, L
        self.dx = L / n
        x = np.arange(n) * self.dx
        self.x, self.y = np.meshgrid(x, x, indexing="ij")

        k1 = 2 * np.pi * np.fft.fftfreq(n, d=self.dx)     # full wavenumbers
        kr = 2 * np.pi * np.fft.rfftfreq(n, d=self.dx)    # real-FFT (last axis)
        self.kx, self.ky = np.meshgrid(k1, kr, indexing="ij")
        self.k2 = self.kx**2 + self.ky**2
        self.k2_inv = np.where(self.k2 == 0, 0.0, 1.0 / np.where(self.k2 == 0, 1.0, self.k2))
        self.kmag = np.sqrt(self.k2)

        # 2/3-rule dealiasing mask
        kmax = (2.0 / 3.0) * np.abs(k1).max()
        self.dealias = (np.abs(self.kx) < kmax) & (np.abs(self.ky) < kmax)

    # --- transforms -----------------------------------------------------
    def fft(self, f):  return np.fft.rfft2(f)
    def ifft(self, F): return np.fft.irfft2(F, s=(self.n, self.n))

    # --- operators (spectral space) ------------------------------------
    def laplacian(self, F):        return -self.k2 * F
    def invert_laplacian(self, F): return -self.k2_inv * F          # solve ∇²ψ = F
    def ddx(self, F):              return 1j * self.kx * F
    def ddy(self, F):              return 1j * self.ky * F

    def jacobian(self, A, B):
        """J(A,B)=A_x B_y - A_y B_x, dealiased. A,B are spectral fields."""
        Ax = self.ifft(self.ddx(A)); Ay = self.ifft(self.ddy(A))
        Bx = self.ifft(self.ddx(B)); By = self.ifft(self.ddy(B))
        J = Ax * By - Ay * Bx
        return self.fft(J) * self.dealias
