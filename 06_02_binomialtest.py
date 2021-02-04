#! /usr/bin/env python
import argparse
from scipy.stats import binom
import matplotlib.pyplot as plt
from scipy.stats import norm
import os
import math
import numpy as np

def testScipy(N, p, outdir):
    fig, ax = plt.subplots(1, 1)
    x = list(range(0,N+1))
    ax.plot(x, binom.pmf(x, N, p), 'bo', ms=8, label='binom pmf')
    ax.plot(x, binom.cdf(x, N, p), 'ro', ms=8, label='binom cdf')
    ax.plot(x, binom.sf(x, N, p), 'go', ms=8, label='binom sf')
    ax.legend(loc='best', frameon=False)
    plt.savefig(outdir + "/distributions.png")

def main(n, p_0, alpha=1/10, outdir='out'):
    testScipy(n, p_0, outdir)
    k_cutoff = binom.ppf(alpha, n, p_0) - 1
    print('k*: {}'.format(k_cutoff))
    print('P(binom < k*): {}'.format(binom.cdf(k_cutoff, n, p_0)))

    fig, ax = plt.subplots(1, 1)
    x = list(range(0,n+1))
    ax.plot(x, binom.cdf(x, n, p_0), 'ro', ms=5)
    ax.plot([0,n], [alpha,alpha], 'g')
    ax.plot([k_cutoff,k_cutoff], [0,1], 'b')
    # TODO alpha and k*
    plt.savefig(outdir + "/{}_{}.png".format(n, p_0))

if __name__ == "__main__":
    main(20, 1/2)
    main(100, 1/7)

