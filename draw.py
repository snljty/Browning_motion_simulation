#! /usr/bin/env python3
# -*- Coding: UTF-8 -*-

r"""
draw MSD(dt) with dt graph of Browning motion simulation.
the data is already prepared by another program.
"""

import numpy as np
import matplotlib.pyplot as plt

with open("Browning_motion_result.txt") as f:
    x, y = np.loadtxt(f, comments = "#", unpack = True)
    
fig, ax = plt.subplots(figsize = (9.6, 4.8))
ax.plot(x, y, "k-")
ax.set_xlabel("dt")
ax.set_ylabel("MSD")
ax.set_title("MSD(dt) with dt")
fig.savefig("MSD_with_dt.png")
plt.show()

