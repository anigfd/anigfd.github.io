"""Wind-driven ocean circulation: the Stommel (1948) and Munk (1950) gyre
problems on a nondimensional [0,1]x[0,1] basin, with idealized single-gyre
wind forcing curl(tau)=-sin(pi*y) (pure NumPy).

Both are boundary-value problems on a BOUNDED domain -- a genuinely new
kind of numerical primitive for this book (every previous chapter used
either a doubly-periodic spectral grid or a channel with only two
boundaries). Each has an exact, separable analytic solution derived here
(not looked up) via psi=f(x)*sin(pi*y), reducing the 2D PDE to a linear
ODE for f(x) with constant coefficients, solved exactly (roots of the
characteristic polynomial via numpy.roots/quadratic formula, integration
constants via a small numpy.linalg.solve) -- and each also has a general
2D finite-difference direct solver (works for ANY forcing, not just the
separable case), cross-validated against the exact solution.

Boundary conditions: no-normal-flow (psi=0) on all four walls, always. The
Munk (biharmonic) solver additionally needs a TANGENTIAL condition: no-slip
(psi_n=0, via Thom's 1933 wall-vorticity formula, derived below) on the
EAST-WEST walls, where the western-boundary-current physics lives, and
free-slip (zeta=0, automatically satisfied by a plain Dirichlet-zero
Laplacian with no correction) on the NORTH-SOUTH walls -- the standard
textbook simplification (Pedlosky, Vallis) that keeps the problem
separable. Getting this wrong (no-slip on ALL four walls) was tried and
caught: it solves a genuinely different, non-separable problem, and
disagreed with the exact analytic solution by a stubborn ~30% at every
resolution and boundary-layer width tested, immediately ruled out as
"just needs finer resolution" by the fact that the disagreement did NOT
shrink with grid refinement (a true discretization error must vanish as
the grid is refined; a mismatched boundary condition does not).
"""
import numpy as np


def wind_stress_curl(y):
    """Idealized single-gyre wind forcing on [0,1]: curl(tau) = -sin(pi*y),
    negative everywhere -- easterlies at low y, westerlies at high y,
    driving one subtropical-gyre-like cell."""
    return -np.sin(np.pi * y)


# ---------------------------------------------------------------------
# Stommel (1948): linear bottom drag, eps*lap(psi)+psi_x=curl(tau)
# ---------------------------------------------------------------------

def stommel_f(x, eps):
    """f(x) for eps*f''+f'-eps*pi^2*f=-1, f(0)=f(1)=0 (from separating
    psi=f(x)sin(pi*y) into the Stommel PDE with the idealized forcing).
    Solved exactly: characteristic quadratic's 2 roots (quadratic formula),
    particular solution by inspection, 2 integration constants from the 2
    BCs via a 2x2 numpy.linalg.solve."""
    pi = np.pi
    a, b, c = eps, 1.0, -eps * pi ** 2
    disc = np.sqrt(b ** 2 - 4 * a * c)
    m1, m2 = (-b + disc) / (2 * a), (-b - disc) / (2 * a)
    fp = 1.0 / (eps * pi ** 2)
    M = np.array([[1.0, 1.0], [np.exp(m1), np.exp(m2)]])
    A, B = np.linalg.solve(M, np.array([-fp, -fp]))
    x = np.asarray(x)
    return A * np.exp(m1 * x) + B * np.exp(m2 * x) + fp


def stommel_analytic(x, y, eps):
    """Exact solution psi(x,y)=f(x)*sin(pi*y) of the Stommel gyre problem
    on [0,1]x[0,1] with the idealized single-gyre wind forcing."""
    return stommel_f(x, eps) * np.sin(np.pi * y)


def _lap1d_dirichlet(n, d):
    return (np.diag(-2.0 * np.ones(n)) + np.diag(np.ones(n - 1), 1)
            + np.diag(np.ones(n - 1), -1)) / d ** 2


def _d1_dirichlet(n, d):
    return (np.diag(np.ones(n - 1), 1) - np.diag(np.ones(n - 1), -1)) / (2 * d)


def solve_stommel_fd(nx, ny, eps, forcing_func=wind_stress_curl):
    """Direct 2D finite-difference solve of eps*lap(psi)+psi_x=F(x,y) on
    [0,1]x[0,1], psi=0 on all four walls: the interior Laplacian and d/dx
    operators are assembled via Kronecker sums (extending the 1D pattern
    already used in gfdlib.instability to 2D) into one dense matrix, solved
    directly with numpy.linalg.solve -- no SciPy sparse solvers. Returns
    psi on the FULL grid (including the zero boundary) and the coordinate
    arrays."""
    x_full = np.linspace(0.0, 1.0, nx)
    y_full = np.linspace(0.0, 1.0, ny)
    dx, dy = x_full[1] - x_full[0], y_full[1] - y_full[0]
    nix, niy = nx - 2, ny - 2

    Lxx, Lyy = _lap1d_dirichlet(nix, dx), _lap1d_dirichlet(niy, dy)
    Dx = _d1_dirichlet(nix, dx)
    Ix, Iy = np.eye(nix), np.eye(niy)
    LAP = np.kron(Iy, Lxx) + np.kron(Lyy, Ix)
    DX = np.kron(Iy, Dx)

    A = eps * LAP + DX
    X, Y = np.meshgrid(x_full[1:-1], y_full[1:-1], indexing="xy")
    rhs = forcing_func(Y).ravel() if forcing_func is wind_stress_curl else forcing_func(X, Y).ravel()

    psi_int = np.linalg.solve(A, rhs)
    psi_full = np.zeros((ny, nx))
    psi_full[1:-1, 1:-1] = psi_int.reshape(niy, nix)
    return psi_full, x_full, y_full


# ---------------------------------------------------------------------
# Munk (1950): lateral (biharmonic) friction,
#              -delta^3*lap(lap(psi))+psi_x=curl(tau)
# ---------------------------------------------------------------------

def _munk_coeffs(delta):
    pi = np.pi
    d3 = delta ** 3
    coeffs = [-d3, 0.0, 2 * d3 * pi ** 2, 1.0, -d3 * pi ** 4]
    roots = np.roots(coeffs)
    fp = 1.0 / (d3 * pi ** 4)
    M = np.array([
        [1.0, 1.0, 1.0, 1.0],
        roots,
        [np.exp(r) for r in roots],
        [r * np.exp(r) for r in roots],
    ])
    C = np.linalg.solve(M, np.array([-fp, 0.0, -fp, 0.0]))
    return C, roots, fp


def munk_f(x, delta):
    """f(x) for -delta^3*f''''+2*delta^3*pi^2*f''+f'-delta^3*pi^4*f=-1,
    f(0)=f'(0)=f(1)=f'(1)=0 (no-slip both ends). Solved exactly: the
    characteristic quartic's 4 roots via numpy.roots (no SciPy), 4
    integration constants from the 4 BCs via a 4x4 numpy.linalg.solve."""
    C, roots, fp = _munk_coeffs(delta)
    x = np.atleast_1d(np.asarray(x, dtype=float))
    shape, flat = x.shape, x.ravel()
    out = np.array([np.real(np.sum(C * np.exp(roots * xx)) + fp) for xx in flat])
    return out.reshape(shape)


def munk_analytic(x, y, delta):
    """Exact solution psi(x,y)=f(x)*sin(pi*y) of the Munk gyre problem on
    [0,1]x[0,1] (no-slip east/west, free-slip north/south -- see module
    docstring) with the idealized single-gyre wind forcing."""
    return munk_f(x, delta) * np.sin(np.pi * y)


def solve_munk_fd(nx, ny, delta, forcing_func=wind_stress_curl):
    """Direct 2D finite-difference solve of -delta^3*lap(lap(psi))+psi_x=
    F(x,y) on [0,1]x[0,1] via a COUPLED (psi,zeta) linear system with
    zeta=lap(psi) kept as an independent field (the same pattern already
    used in ch16/ch17/ch19), rather than eliminating it into a
    hand-derived biharmonic stencil. No-slip on east/west walls is
    supplied by Thom's (1933) wall-vorticity formula zeta_wall = 2*
    psi_adjacent/d^2 (derived in the module docstring's spirit by
    Taylor-expanding psi from a psi=psi_n=0 wall): this adds an extra
    coefficient onto the psi-block for cells adjacent to those walls,
    exactly compensating for the "missing" wall-zeta term the plain
    Dirichlet-zero Laplacian would otherwise silently drop. North/south
    walls use free-slip (zeta=0), which the plain Laplacian already gives
    with no correction needed. Solved directly with numpy.linalg.solve --
    no SciPy. Returns psi on the FULL grid and the coordinate arrays."""
    x_full = np.linspace(0.0, 1.0, nx)
    y_full = np.linspace(0.0, 1.0, ny)
    dx, dy = x_full[1] - x_full[0], y_full[1] - y_full[0]
    nix, niy = nx - 2, ny - 2
    n_int = nix * niy

    Lxx, Lyy = _lap1d_dirichlet(nix, dx), _lap1d_dirichlet(niy, dy)
    Dx = _d1_dirichlet(nix, dx)
    Ix, Iy = np.eye(nix), np.eye(niy)
    LAP = np.kron(Iy, Lxx) + np.kron(Lyy, Ix)
    DX = np.kron(Iy, Dx)

    def idx(i, j):
        return j * nix + i

    W = np.zeros((n_int, n_int))
    for j in range(niy):
        for i in range(nix):
            row = idx(i, j)
            if i == 0 or i == nix - 1:          # no-slip WEST / EAST
                W[row, idx(i, j)] += 2.0 / dx ** 4
            # free-slip NORTH/SOUTH: no correction needed (zeta=0 already)

    Iint = np.eye(n_int)
    A = np.block([
        [LAP, -Iint],
        [DX - delta ** 3 * W, -delta ** 3 * LAP],
    ])

    X, Y = np.meshgrid(x_full[1:-1], y_full[1:-1], indexing="xy")
    F = (forcing_func(Y) if forcing_func is wind_stress_curl else forcing_func(X, Y)).ravel()
    rhs = np.concatenate([np.zeros(n_int), F])

    sol = np.linalg.solve(A, rhs)
    psi_full = np.zeros((ny, nx))
    psi_full[1:-1, 1:-1] = sol[:n_int].reshape(niy, nix)
    return psi_full, x_full, y_full
