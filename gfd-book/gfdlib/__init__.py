"""
gfdlib — shared, Pyodide-safe primitives for the Interactive GFD textbook.

Pure NumPy only (no SciPy, no compiled extensions) so every notebook runs in
the browser via Pyodide/WebAssembly. Import the pieces you need:

    from gfdlib import spectral, timestep, diagnostics, plotting, shallowwater, internalwaves, convection, pv, balance, rossby, qg, instability, baroclinic, symmetric

Design rules (see CLAUDE.md):
  * one vetted implementation of each numeric primitive, reused everywhere;
  * precompute wavenumber arrays ONCE, outside the time loop;
  * sign/notation conventions fixed centrally (see NOTATION.md).
"""
from . import spectral, timestep, diagnostics, plotting, shallowwater, internalwaves, convection, pv, balance, rossby, qg, instability, baroclinic, symmetric

__all__ = ["spectral", "timestep", "diagnostics", "plotting", "shallowwater", "internalwaves", "convection", "pv", "balance", "rossby", "qg", "instability", "baroclinic", "symmetric"]
__version__ = "0.1.0"
