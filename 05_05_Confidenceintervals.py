#! /usr/bin/env python
import argparse
from scipy.stats import binom
import matplotlib.pyplot as plt
from scipy.stats import norm
import os
import math
import numpy as np

def f_range(start, stop, step):
    while start <= stop:
        yield start
        start += step

def q_u_norm(n, k, beta):
    return k/n - getZfromBeta(beta) * math.sqrt(1/n * k/n * (1 - k/n))

def q_o_norm(n, k, beta):
    return k/n + getZfromBeta(beta) * math.sqrt(1/n * k/n * (1 - k/n))

def get_epsilon(n,beta):
    return 1/(2 * math.sqrt(2*n*beta))

def q_o_tscheb(n,k,beta):
    return k/n + get_epsilon(n,beta)

def q_u_tscheb(n,k,beta):
    return k/n - get_epsilon(n,beta)

def get_P_p_gt(n, p, q_o):
    vals = []
    for k in range(len(q_o)):
        if p > q_o[k]:
            vals.append(binom.pmf(k,n,p))
        else:
            vals.append(0)
    return np.sum(vals)

def get_P_p_lt(n, p, q_u):
    vals = []
    for k in range(len(q_u)):
        if p < q_u[k]:
            vals.append(binom.pmf(k,n,p))
        else:
            vals.append(0)
    return np.sum(vals)

def get_P_p_outs(n, p, q_u, q_o):
    vals = []
    for k in range(len(q_o)):
        if p < q_u[k] or p > q_o[k]:
            vals.append(binom.pmf(k,n,p))
        else:
            vals.append(0)
    return np.sum(vals)

def testScipy(N, p, outdir):
    fig, ax = plt.subplots(1, 1)
    x = list(range(0,N+1))
    ax.plot(x, binom.pmf(x, N, p), 'bo', ms=8, label='binom pmf')
    ax.plot(x, binom.cdf(x, N, p), 'ro', ms=8, label='binom cdf')
    ax.plot(x, binom.sf(x, N, p), 'go', ms=8, label='binom sf')
    ax.legend(loc='best', frameon=False)
    plt.savefig(outdir + "/distributions.png")

def getZfromBeta(beta):
    return norm.ppf(1-beta)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot several confidence intervals')
    parser.add_argument('N', type=int, help='Number of element for binomial distribution')
    parser.add_argument('beta', type=float, help='Expected 1 - confidence')
    parser.add_argument('precision', nargs='?', type=float, help='Step size to use for p', default=0.01)
    parser.add_argument('outdir', nargs='?', type=str, help='output directory', default='out')
    args = parser.parse_args()

    testScipy(10,0.2,args.outdir)
    print("z: {}".format(getZfromBeta(args.beta)))

    x = list(range(0, args.N+1))
    q_u_norm_v = [q_u_norm(args.N, i, args.beta) for i in x]
    q_o_norm_v = [q_o_norm(args.N, i, args.beta) for i in x]
    p_range = list(f_range(0,1,args.precision))
    P_p_gt_norm = [get_P_p_gt(args.N, p, q_o_norm_v) for p in p_range]
    P_p_lt_norm = [get_P_p_lt(args.N, p, q_u_norm_v) for p in p_range]
    P_p_outs_norm = [get_P_p_outs(args.N, p, q_u_norm_v, q_o_norm_v) for p in p_range]
    fig, ax = plt.subplots(1, 1)
    ax.plot(p_range, [args.beta]*len(p_range), 'k', ms=1, label='beta')
    ax.plot(p_range, [2 * args.beta]*len(p_range), 'k', ms=1, label='2 beta')
    ax.plot(p_range, P_p_gt_norm, 'b', label='P_p_gt_norm')
    ax.plot(p_range, P_p_lt_norm, 'r', label='P_p_lt_norm')
    ax.plot(p_range, P_p_outs_norm, 'g', label='P_p_outs_norm')
    ax.legend(loc='best', frameon=False)
    plt.savefig(args.outdir + "/norm_{}_{}.png".format(args.N, args.beta))

    q_u_tscheb_v = [q_u_tscheb(args.N, i, args.beta) for i in x]
    q_o_tscheb_v = [q_o_tscheb(args.N, i, args.beta) for i in x]
    P_p_gt_tscheb = [get_P_p_gt(args.N, p, q_o_tscheb_v) for p in p_range]
    P_p_lt_tscheb = [get_P_p_lt(args.N, p, q_u_tscheb_v) for p in p_range]
    P_p_outs_tscheb = [get_P_p_outs(args.N, p, q_u_tscheb_v, q_o_tscheb_v) for p in p_range]
    fig, ax = plt.subplots(1, 1)
    ax.plot(p_range, [args.beta]*len(p_range), 'k', ms=1, label='beta')
    ax.plot(p_range, [2 * args.beta]*len(p_range), 'k', ms=1, label='2 beta')
    ax.plot(p_range, P_p_gt_tscheb, 'b', label='P_p_gt_tscheb')
    ax.plot(p_range, P_p_lt_tscheb, 'r', label='P_p_lt_tscheb')
    ax.plot(p_range, P_p_outs_tscheb, 'g', label='P_p_outs_tscheb')
    ax.legend(loc='best', frameon=False)
    plt.savefig(args.outdir + "/tscheb_{}_{}.png".format(args.N, args.beta))


    try:
        os.makedirs(args.outdir)
    except FileExistsError:
        pass


