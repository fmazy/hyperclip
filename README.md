# hyperclip

This Python 3.5+ package implements volume computation of hypercubes clipped by hyperplanes.
All methods implemented here have been proposed by Yunhi Cho and Seonhwa Kim (2020) in the article [Volume of Hypercubes Clipped by Hyperplanes and Combinatorial Identities](https://doi.org/10.13001/ela.2020.5085). An arxiv paper is available [here](https://arxiv.org/pdf/1512.07768.pdf).
The documentation is available on [Read the Doc](https://hyperclip.readthedocs.io/en/latest/).

## Installation

Hyperclip is available through [PyPI](https://pypi.org/project/hyperclip/), and may be installed using `pip` :
   
    pip install hyperclip

## Introduction

The package is essentially composed by two classes : `hyperclip.Hyperplane` and `hyperclip.Hyperclip`.

* `hyperclip.Hyperplane` allows users to create a n-dimensional hyperplane defined as `a.x + r \geq 0`. It is possible to directly set `a` and `r` or to provide `n` distinct points which belongs to the hyperplane, i.e `a.x + r = 0`.
* `hyperclip.Hyperclip` allows users to create an hyperclip object. It aims to compute the volume of `A.X+R \leq 0` for `X` inside the uniform hypercube `[0,1]^n`. It is possible to directly set `A` and `R` or to set a list of `hyperclip.Hyperplane` objects.

## Example code

Here's an example showing the usage of `hyperclip.Hyperclip` for a 2-dimensional case.
The result provided by Hyperclip is compared to a MonteCarlo volume estimation.
The source code is available on [the documentation main page](https://hyperclip.readthedocs.io/en/latest/).

![example_figure](docs/source/figures/example_2d.png)

## Contact

Please, send me an email at [francois-remi.mazy@inria.fr](mailto:francois-remi.mazy@inria.fr).