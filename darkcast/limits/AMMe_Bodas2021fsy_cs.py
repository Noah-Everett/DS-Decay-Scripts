# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2023 DARKCAST authors (see AUTHORS.md).
import darkcast
notes = """
This limit is set from the difference between the measured
electron anomalous magnetic moment and SM prediction, and was
calculated using AMM.py and the values of equation 2 from
Bodas:2021fsy. The fine-structure constant for this limit was
determined with cesium atoms and calculated at 3 sigma from the
nominal value.

Given the nature of the limit, the efficiency ratio is unity, e.g. t0
= 0 and t1 = infinity.
"""
bibtex = """
@article{Bodas:2021fsy,
 author = "Bodas, Arushi and Coy, Rupert and King, Simon J. D.",
 title = "{Solving the electron and muon $g-2$ anomalies in $Z'$ models}",
 eprint = "2102.07781",
 archivePrefix = "arXiv",
 primaryClass = "hep-ph",
 reportNumber = "ULB-TH/21-01",
 month = "2",
 year = "2021"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production("e_e")
decay = "none"
bounds = darkcast.Dataset("limits/AMMe_Bodas2021fsy_cs.lmt")
efficiency = darkcast.Efficiency(t0 = 0, t1 = float("inf"))
valid = [False, True]
