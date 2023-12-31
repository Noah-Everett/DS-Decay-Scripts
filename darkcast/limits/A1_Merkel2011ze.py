# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2023 DARKCAST authors (see AUTHORS.md).
import darkcast
notes = """
This limit was extracted from figure 5 (gray curve labeled A1) of
Blumlein:2013cua.

No detailed information is given on the prompt-like
requirements. However, since the same coupling is used in production
and decay, the efficiency ratio is assumed to be unity, e.g. t0 = 0
and t1 = infinity.
"""
bibtex = """
@article{Merkel:2014avp,
 author         = "Merkel, H. and others",
 title          = "{Search at the Mainz Microtron for Light Massive Gauge
                   Bosons Relevant for the Muon g-2 Anomaly}",
 journal        = "Phys. Rev. Lett.",
 volume         = "112",
 year           = "2014",
 number         = "22",
 pages          = "221802",
 doi            = "10.1103/PhysRevLett.112.221802",
 eprint         = "1404.5502",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ex",
 SLACcitation   = "%%CITATION = ARXIV:1404.5502;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production("e_brem")
decay = "e_e"
bounds = darkcast.Dataset("limits/A1_Merkel2011ze.lmt")
efficiency = darkcast.Efficiency(t0 = 0, t1 = float("inf"))
