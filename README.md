[![Build Status](https://travis-ci.org/berquist/obarasaika.svg?branch=master)](https://travis-ci.org/berquist/obarasaika/builds)

# Obara-Saika

Obara-Saika molecular integral code sandbox.

This is a terribly slow but in principle open-ended recursive implementation of the [Obara-Saika algorithm](http://dx.doi.org/10.1063/1.450106) using Python.

It may be useful for educational purposes or to generate and/or test a more efficient code.

The computed integrals are not normalized.

Licensed under [BSD-3](../master/LICENSE).

Originally forked from [here](https://github.com/bast/obara-saika) with many thanks for providing the building blocks.

## Dependencies

* Python (tested 2.6, 2.7, 3.2, 3.3, 3.4, 3.5, 3.6)
* `mpmath`
* For testing, either nose or pytest should work; check the Makefile.
