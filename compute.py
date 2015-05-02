from __future__ import print_function

from itertools import product

import shell
import obara_saika as os

za = 1.1
zb = 1.2
zc = 1.3
zd = 1.4

ra = [1.0, 0.0, 1.0]
rb = [0.0, 1.0, 2.0]
rc = [0.0, 0.0, 3.0]
rd = [0.0, 0.0, 4.0]

q_base  = [0, 0, 0, 0, 0, 0]
q_small = [0, 1, 0, 0, 0, 1]
q_large = [1, 1, 1, 2, 2, 2]


def compute_coulomb():
    # get all integrals up to pppp
    for p in product('01', repeat=4):
        for c in shell.get_shell4(int(p[0]), int(p[1]), int(p[2]), int(p[3])):
            if sum(c) > 0:
                print(c, os.get_coulomb(za, zb, zc, zd, ra, rb, rc, rd, c))


def compute_overlap():
    print('# compute_overlap_base')
    S = os.get_overlap(za, zb, ra, rb, q_base)
    print(S)
    print('# compute_overlap_small')
    S = os.get_overlap(za, zb, ra, rb, q_small)
    print(S)
    print('# compute_overlap_large')
    S = os.get_overlap(za, zb, ra, rb, q_large)
    print(S)


def compute_kinetic():
    print('# compute_kinetic_base')
    T = os.get_kinetic(za, zb, ra, rb, q_base)
    print(T)
    print('# compute_kinetic_small')
    T = os.get_kinetic(za, zb, ra, rb, q_small)
    print(T)
    print('# compute_kinetic_large')
    T = os.get_kinetic(za, zb, ra, rb, q_large)
    print(T)


def compute_nuclear():
    print('# compute_nuclear_base')
    V = os.get_nuclear(za, zb, ra, rb, rc, q_base)
    print(V)
    print('# compute_nuclear_small')
    V = os.get_nuclear(za, zb, ra, rb, rc, q_small)
    print(V)
    print('# compute_nuclear_large')
    V = os.get_nuclear(za, zb, ra, rb, rc, q_large)
    print(V)

def compute_moment():
    print('# compute_moment_base')
    M = os.get_moment(za, zb, ra, rb, rc, q_base, [1, 1, 1])
    print(M)
    print('# compute_moment_small')
    M = os.get_moment(za, zb, ra, rb, rc, q_small, [1, 1, 1])
    print(M)
    print('# compute_moment_large')
    M = os.get_moment(za, zb, ra, rb, rc, q_large, [1, 1, 1])
    print(M)


if __name__ == '__main__':
    compute_overlap()
    compute_kinetic()
    compute_nuclear()
    compute_moment()
