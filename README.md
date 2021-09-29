# Algorithmic Lovász Local Lemma (ALLL) Satisfiability Solver

---

Xandru Mifsud (2021)

---

## Requirements

This program is intended to run on Linux platforms, in particular on Debian-based systems such as Ubuntu and
Lubuntu, on which we have tested the implementation with relative success.

For installation, ```cmake``` version 3.17+ is required, as well as ```Eigen``` version 3.3 and above, and the ```CNF_IO```
package available at: https://people.sc.fsu.edu/~jburkardt/cpp_src/cnf_io/cnf_io.html

## Installation Instructions

Clone the repository, and ```cd``` into the project directory. Then run:

1. ```cmake .```
2. ```make```

## Execution Instructions

Simply ```cd``` into the directory containing the compiled executable, and run ```./ALLLSatisfiabilitySolver <path/to/sat_instance.cnf>```,
where the required filepath ```path/to/sat_instance.cnf``` is the path to the CNF SAT instance to be loaded, represented
in the DIMACS format.