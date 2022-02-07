# Algorithmic Lovász Local Lemma (ALLL) Satisfiability Solver

---

Xandru Mifsud (2021)

---

## Requirements

This program is intended to run on 64-bit Linux and macOS (Intel and Apple Silicon) based systems, on which it has been tested thoroughly.

For installation, ```cmake``` version 3.17+ is required, as well as ```libomp```.

## Installation Instructions

Clone the repository, and ```cd``` into the project directory. Then run:

1. ```mkdir ./build```
1. ```cmake ..```
2. ```make```

## Execution Instructions

Simply ```cd``` into the directory containing the compiled executable, and run ```./ALLLSatisfiabilitySolver [-p [n_threads]] <path/to/sat_instance.cnf>```, where the required filepath ```path/to/sat_instance.cnf``` is the path to the CNF SAT instance to be loaded, represented in the DIMACS format, and ```-p``` is an optional positional argument specifying that parallelisation should be used wherever possible. By default, if ```n_threads``` is not specified, this will use all available (physical) threads on the system. Otherwise, ```n_threads``` are used for parallelisation.
