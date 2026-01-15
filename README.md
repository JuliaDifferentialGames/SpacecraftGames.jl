# SpacecraftGameBenchmarks.jl

[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)
[![Build Status](https://github.com/BennetOutland/SpacecraftGameBenchmarks.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/BennetOutland/SpacecraftGameBenchmarks.jl/actions/workflows/CI.yml?query=branch%3Amain)


A Julia package for benchmarking differential game solvers applied to space robotics built on [DifferentialGames.jl](https://github.com/JuliaDifferentialGames/DifferentialGames.jl).

## Purpose

This package puts forth a set of benchmark, or gym, environments for evaluating differential game solvers applied to problems in space robotics. 

## Installation

```julia
using Pkg
Pkg.add("https://github.com/JuliaDifferentialGames/SpacecraftGameBenchmarks.jl.git")
```

## Development Status

**Planned:**
- 📋 Hill-Clohessy-Wiltshire Pursuit-Evasion: N vs. M Player
- 📋 Pursuit-Evasion: N vs. M Player
- 📋 Attacker-Defender Race-to-Target with Keep-out Zones: 1 vs. 1 Player
- 📋 Sun Blocking: 1 vs. 1 Player
- 📋 Lady-Bandit-Guard Problem: 1 vs. 1 Player


## License

MIT License - see LICENSE file for details.


## Disclosure of Generative AI Usage

Generative AI, Claude, was used in the creation of this library as a programming aid including guided code generation, assistance with performance optimization, and for assistance in writing documentation. All code and documentation included in this repository, whether written by the author(s) or generative AI, has been reviewed by the author(s) for accuracy.
