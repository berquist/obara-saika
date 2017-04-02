from __future__ import print_function

import timeit

import obarasaika.obarasaika.obara_saika as os


za = 1.1
zb = 1.2

ra = [1.0, 0.0, 1.0]
rb = [0.0, 1.0, 2.0]
rc = [0.0, 0.0, 3.0]

q_base  = [0, 0, 0, 0, 0, 0]
q_small = [0, 1, 0, 0, 1, 0]
q_large = [4, 4, 4, 4, 4, 4]


def compute_overlap_base():
    os.get_overlap(za, zb, ra, rb, q_base)


def compute_kinetic_base():
    os.get_overlap(za, zb, ra, rb, q_base)


def compute_nuclear_base():
    os.get_nuclear(za, zb, ra, rb, rc, q_base)


if __name__ == '__main__':

    timer_overlap = timeit.Timer('compute_overlap_base()')
    timer_kinetic = timeit.Timer('compute_kinetic_base()')
    timer_nuclear = timeit.Timer('compute_nuclear_base()')

    # results_overlap = timer_overlap.repeat()
    # results_kinetic = timer_kinetic.repeat()
    # results_nuclear = timer_nuclear.repeat()
