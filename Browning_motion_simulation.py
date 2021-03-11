#! python3
# -*- Coding: UTF-8 -*-

r"""
Simulation of Browning motion.

The MSD(t) is defined as:
MSD(t) = <|r(t) - r(0)|^2> .

There should be a relationship like:
MSD(t) / t = const .

"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import tqdm

Get_SD = lambda arr, pos, interval: np.sum((arr[pos, :] - \
                                            arr[pos + interval, :]) ** 2)

Get_MSD = lambda arr, interval: np.sum([Get_SD(arr, pos, interval) \
                                     for pos in np.arange(nsize - interval)]) \
                                / (nsize - interval)

def Get_all_MSD(arr):
    msd = np.empty(arr.shape[0])
    msd[0] = 0.0
    for dt in np.arange(1, arr.shape[0]):
        msd[dt] = Get_MSD(arr, dt)
    return msd

def Auto_Corr_FFT(val):
    N = val.shape[0]
    F = np.fft.fft(val, n = 2 * N)
    PSD = F * F.conjugate()
    res = np.fft.ifft(PSD)
    res = (res[:N]).real
    n = N * np.ones(N) - np.arange(0, N)
    return res / n

def Get_all_MSD_FFT(arr):
    N = arr.shape[0]
    D = np.square(arr).sum(axis = 1) 
    D = np.append(D, 0)
    S2 = sum([Auto_Corr_FFT(arr[:, i]) for i in np.arange(arr.shape[1])])
    Q = 2 * D.sum()
    S1 = np.zeros(N)
    for m in np.arange(N):
        Q = Q - D[m - 1] - D[N - m]
        S1[m] = Q / (N - m)
    return S1 - 2 * S2

parser = argparse.ArgumentParser(prog = "Browning_motion_simulation", 
    description = "simulates Browning motion and calculates msd", epilog = 
"""
This program uses fft to accelerate the MSD calculation.

There should be a proportional relationship between MSD(dt) and dt.

The default
""")

nsize    =  1000
ndim     =     2
stepsize =    10
nsimu    = 10000

parser.add_argument("--nsize", type = int, default = nsize, 
    help = "steps recorded by the simulation, default to {nsize:d}".format(nsize = nsize))
parser.add_argument("--ndim", type = int, default = ndim, 
    help = "amount of dimensions,             default to {ndim:d}".format(ndim = ndim))
parser.add_argument("--stepsize", type = int, default = stepsize, 
    help = "steps between each recored steps, default to {stepsize:d}".format(stepsize = stepsize))
parser.add_argument("--nsimu", type = int, default = nsimu, 
    help = "amount of simulations,            default to {nsimu:d}".format(nsimu = nsimu))

args = parser.parse_args()
nsize = args.nsize
ndim = args.ndim
stepsize = args.stepsize
nsimu = args.nsimu

print("nsize = {nsize:d}, ndim = {ndim:d}, stepsize = {stepsize:d}, nsimu = {nsimu:d}".format \
    (nsize = nsize, ndim = ndim, stepsize = stepsize, nsimu = nsimu))
x = np.zeros((nsize, ndim))
MSD = np.zeros(nsize, dtype = np.float64)

for isimu in tqdm.trange(1, 1 + nsimu):
    x[:, :] = 0.0
    for istep in np.arange(stepsize):
        x[1:, :] += np.cumsum(np.random.normal(size = (nsize - 1, ndim)), axis = 0)
    # MSD += Get_all_MSD(x)
    MSD += Get_all_MSD_FFT(x)

MSD /= nsimu

def Linear_through_origin_fit(arr: np.ndarray, start_pos = None, end_pos = None) -> (np.float64, np.float64):
    r"""
fit a line through origin
"""
    if (arr.ndim != 1): raise ValueError("Dimension of array to be fitted is no 1.")
    if not start_pos: start_pos = round(arr.size * 0.1)
    if not   end_pos:   end_pos = round(arr.size * 0.9)

    var = np.arange(start_pos, end_pos)
    slope = sum(arr[start_pos:end_pos] * var) / sum(var ** 2)
    arr_average = np.average(arr[start_pos:end_pos])
    cod = 1.0 - sum((arr[start_pos:end_pos] - var * slope) ** 2) / \
                sum((arr[start_pos:end_pos] - arr_average) ** 2)
    return slope, cod

slope, cod = Linear_through_origin_fit(MSD)
co_diffu = slope / (2 * ndim)
print("Statistic is from 10% to 90%")
print("R**2 =  %7.4lf" % cod)
print("Diffusivity = %11.4le" % co_diffu)

fig, ax = plt.subplots(figsize = (9.6, 4.8))

ax.plot(np.arange(nsize), MSD, "k-")
ax.set_title("MSD with \u0394t, ndim = %d, stepsize = %d, Diffusivity = %.4lf" % (ndim, stepsize, co_diffu))
ax.set_xlabel(u"\u0394t")
ax.set_ylabel("MSD(\u0394t)")
fig.savefig(u"MSD_with_dt.png")
plt.show()

