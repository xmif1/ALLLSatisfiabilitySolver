# Algorithmic Lovász Local Lemma (ALLL) Satisfiability Solver

---

Xandru Mifsud (2021)

---

## Requirements

This program is intended to run on Linux platforms, in particular on Debian-based systems such as Ubuntu and
Lubuntu, on which we have tested the implementation with relative success.

For installation, ```cmake``` version 3.17+ is required.

## Installation Instructions

Clone the repository, and ```cd``` into the project directory. Then run:

1. ```cmake .```
2. ```make```

## Execution Instructions

Simply ```cd``` into the directory containing the compiled executable, and run ```./ALLLSatisfiabilitySolver [-p] <path/to/sat_instance.cnf>```,
where the required filepath ```path/to/sat_instance.cnf``` is the path to the CNF SAT instance to be loaded, represented
in the DIMACS format, and ```-p``` is an optional positional argument specifying that parallelisation should be
used wherever possible.