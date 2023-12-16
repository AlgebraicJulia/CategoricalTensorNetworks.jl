# CategoricalTensorNetworks.jl

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://AlgebraicJulia.github.io/CategoricalTensorNetworks.jl/stable)
[![Development Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://AlgebraicJulia.github.io/CategoricalTensorNetworks.jl/dev)
[![Code Coverage](https://codecov.io/gh/AlgebraicJulia/CategoricalTensorNetworks.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/AlgebraicJulia/CategoricalTensorNetworks.jl)
[![CI/CD](https://github.com/AlgebraicJulia/CategoricalTensorNetworks.jl/actions/workflows/julia_ci.yml/badge.svg)](https://github.com/AlgebraicJulia/CategoricalTensorNetworks.jl/actions/workflows/julia_ci.yml)

A Julia package for tensor networks from a categorical point of view. Tensor
networks are viewed as an algebra of the operad of undirected wiring diagrams
(UWDs).

**Warning**: This repository consists of experimental code plus a few
illustrative examples. It is not ready for serious use but, with some additional
work, could become the basis of a useful package.

## Relevant papers

- Spivak et al, 2016: Pixel Arrays: A fast and elementary method for solving
  nonlinear systems ([arXiv:1609.00061](https://arxiv.org/abs/1609.00061))
  
  Introduces unconventional numerical methods based on the idea that relations
  ("pixel arrays") form an algebra of the operad of UWDs. In this package,
  pixel arrays are tensors valued in the boolean rig.

- Biamonte et al, 2011: Categorical Tensor Network States
  ([arXiv:1012.0531](https://arxiv.org/abs/1012.0531))
  
  Explains tensor networks in terms of string diagrams and other categorical
  concepts in a way that should be intelligble to physicists
