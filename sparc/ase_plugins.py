from __future__ import annotations

# ASE 4 experimental plugin declaration for SPARC-X-API.
# This module is loaded through the ``ase.plugins`` entry point and uses the
# ``ase._4.plugins`` plugin system introduced for ASE 4 calculator discovery.

try:
    from ase._4.plugins.calculator import CalculatorPlugin
except ImportError:
    # ASE < 4 does not provide the new plugin system.  Older ASE versions will
    # not consume the ``ase.plugins`` entry point, but keeping this module
    # importable avoids surprising failures if users/tools import it directly.
    __ase_plugins__ = set()
else:
    sparc_plugin = CalculatorPlugin(
        name="sparc",
        long_name="Python API for SPARC-X code (FileIO and SocketIO)",
        citation="SPARC-X / SPARC-X-API developers",
        implementation="sparc.calculator.SPARC",
        configurable=True,  # uses ase.config / external executable config
    )

    __ase_plugins__ = {sparc_plugin}
