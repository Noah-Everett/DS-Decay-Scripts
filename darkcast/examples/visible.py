# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2023 DARKCAST authors (see AUTHORS.md).

# This example recasts every visible limit available in
# darkcast/limits to every model available in darkcast/models. If
# matplotlib is available these limits are then plotted. The produced
# figures correspond to figures 4 through 7 of Ilten:2018crw.

# Update the system path to find the DarkCast module.
# This assumes that 'examples' is in 'darkcast/examples.'
import sys, os, inspect, itertools
sys.path.insert(1, os.path.join(os.path.dirname(os.path.realpath(
                inspect.getfile(inspect.currentframe()))), "../../"))

# Import the DarkCast module.
import darkcast

# Load all the available models in darkcast/models and DARKCAST_MODEL_PATH.
models = darkcast.Models()

# Load all the available limits in darkcast/limits and DARKCAST_LIMIT_PATH.
limits = darkcast.Limits()

# Alternatively, all the limits from a folder, '/foo/bar', could be
# loaded as:
#
# limits = darkcast.Limits("/foo/bar")
#
# Note that any limits in the directory that are not valid will not be
# loaded. A single limit with name 'foo_bar' can be loaded as:
#
# limit = darkcast.Limit("foo_bar")

# Try to load matplotlib.
try: import matplotlib.pyplot as pyplot
except: pyplot = None
colors = ["red", "green", "blue", "orange", "magenta", "cyan", "gray"]

# Create the directory for recasted limits.
if not os.path.exists("recast/limits"): os.makedirs("recast/limits")

# Loop over all the models.
for name, model in models.items():

    # Create the recasted limit directory for this model.
    print("Recasting limits to the %s model." % name)
    if not os.path.exists("recast/limits/" + name):
        os.makedirs("recast/limits/" + name)
    
    # If possible, initialize the plot.
    if pyplot:
        fig, ax = pyplot.subplots()
        icolor, lbls = itertools.cycle(colors), {}
        
    # Loop over the limits.
    for label, limit in limits.items():
        if (limit.model.width("invisible", 1) != 0 and
            not limit.production.name.endswith("_scat")): continue
        else: print(label)
        
        # Recast the limit, this returns an object of type 'Datasets' if valid.
        recast = limit.recast(model)
        if recast == None: continue

        # Save the limit to a text file. This is done with the
        # 'Datasets.write' method.
        recast.write("recast/limits/%s/%s.lmt" % (name, label))
            
        # Plot. The 'Datasets.plots' method returns formatted
        # lists of x and y points which can be easily passed to a
        # plotting package.
        if pyplot:
            for x, y in recast.plots():
                lbl = darkcast.utils.latex(limit)
                if not lbl in lbls: c = next(icolor); lbls[lbl] = c
                else: c = lbls[lbl]; lbl = None
                f = not (label.endswith("_2l") or label.endswith("_2u"))
                ax.fill(x, y, label = lbl, alpha = 0.3, fill = f, color = c)

    # Save the plot.
    if pyplot:
        legend = ax.legend(loc = "best", ncol = 2, fontsize = 10)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim([2e-3, 1e2])
        ax.set_ylim([1e-8, 1e1])
        ax.set_xlabel("mass [GeV]")
        ax.set_ylabel("g")
        ax.set_title(darkcast.utils.latex(model.name))
        darkcast.utils.logo()
        fig.savefig("visible_%s.pdf" % name)
