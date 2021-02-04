#! /usr/bin/env python
import argparse
from scipy.stats import chi2
import matplotlib.pyplot as plt
from scipy.stats import norm
import os
import math
import numpy as np
import random

def sampleFair():
    return random.randrange(0,6,1)

def sampleCrucket():
    return abs(random.randrange(-3,6,1))

def accumulate(n, f):
    r = [0] * 6
    for i in range(n):
        r[f()] += 1
    return r

def getChiSquared(n, f):
    r = accumulate(n, f)
    v = 0
    for x in r:
        c = (x - n/6)
        v += c * c / (n/6)
    return v

def getQuantile(alpha, df=5):
    return chi2.ppf(1-alpha, df)

def plotSampleCluster(m, n, f, outdir, alpha, desc, df=5):
    r = [getChiSquared(n, f) for x in range(m)]
    r.sort()
    fig, ax = plt.subplots(1, 1)
    x = np.linspace(0, r[-1], m)
    ax.plot(x, chi2.cdf(x, df), 'r', label='chi2 cdf')
    ax.plot(r, np.linspace(0,1,m), 'b', label='{} m={} n={} cdf'.format(desc, m, n))
    ax.plot([0,r[-1]], [1-alpha,1-alpha], 'y')
    ax.plot([getQuantile(alpha),getQuantile(alpha)], [0,1], 'b')
    ax.legend(loc='best', frameon=False)
    plt.savefig(outdir + "/chi2_{}_{}_{}.png".format(desc, m, n))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Empirical look into Chi Squared distributions')
    parser.add_argument('m', type=int, help='Number of samples to perform')
    parser.add_argument('alpha', nargs='?', type=float, help='significance value to mark in the distribution', default=0.05)
    parser.add_argument('outdir', nargs='?', type=str, help='output directory', default='out')
    args = parser.parse_args()
    N = [6, 30, 150]
    for n in N:
        plotSampleCluster(args.m, n, sampleFair, args.outdir, args.alpha, 'fair')
        plotSampleCluster(args.m, n, sampleCrucket, args.outdir, args.alpha, 'broken')

