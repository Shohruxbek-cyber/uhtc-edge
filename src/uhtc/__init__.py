"""uhtc: coupled thermal and oxidation response of a ZrB2-SiC stagnation region in hypersonic flight.

Modules: atmosphere (US 1976), constants, thermo (NIST-JANAF), airchem (equilibrium air),
stagnation (free stream -> stagnation state), heating (Sutton-Graves, Fay-Riddell, distributions),
materials (ZrB2-SiC and scale properties), evaporation (Hertz-Knudsen-Langmuir + boundary-layer
diffusion), oxidation (layer model), conduction1d (enthalpy solver), conduction2d (P1 finite
elements), coupled (time integration of the coupled problem), trajectory.
Every literature value carries its source in INPUTS.md and next to the code that uses it.
"""
__version__ = "0.1.0"
