# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2023 DARKCAST authors (see AUTHORS.md).
import darkcast
notes = """
"""
bibtex = """
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production("pi0_gamma")
decay = ["e_e","mu_mu"]
bounds = darkcast.Dataset("reach/SciBooNE_1.lmt")
efficiency = darkcast.Efficiency(rvals=True)