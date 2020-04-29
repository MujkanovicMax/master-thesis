#!/usr/bin/env python

import argparse
import numpy as np

def gen_phis(Nphi):
    phis = np.linspace(0,360,Nphi+1)[:-1]
    phis = np.round(phis, decimals=2)
    wgts = 2*np.pi*np.ones_like(phis) / Nphi
    return phis, wgts

def gen_mus(Nmus,a,b):
    mus,wgts = np.polynomial.legendre.leggauss(int(Nmus))
    mus = (b - a)/2. * mus + (b + a)/2.
    wgts = wgts * (b - a)/2.
    return mus,wgts


def _main():
    parser = argparse.ArgumentParser(description="generate sampling points for radiance computations")
    parser.add_argument('-Nphi', type=int)
    parser.add_argument('-Nmu', nargs="+", type=float)
    

    args = parser.parse_args()
    
    angles = None
    wgts = None

    if args.Nphi is not None:
        angles, wgts = gen_phis(args.Nphi)

    if args.Nmu is not None:
        angles, wgts = gen_mus(args.Nmu[0], args.Nmu[1], args.Nmu[2])

    print(' '.join(map(str, angles)))
    print(' '.join(map(str, wgts)))

if __name__ == "__main__":
    _main()
