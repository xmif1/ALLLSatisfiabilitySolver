# Algorithmic Lovász Local Lemma (ALLL) Satisfiability Solver

Xandru Mifsud (2022)

Undergraduate CS APT 'Random Boolean Satisfiability' (University of Malta)

Supervisor: Dr. Sandro Spina

## Description

The following repository consists of serial and parallel implementations of the SAT solver described by [1], implemented
with consultation of [2] and [3].

__The code here forms the main code artefact forming part of the CS APT, as required for the awarding of the B.Sc.
(Hons) in Computer Science and Mathematics.__

## Requirements

This program is intended to run on 64-bit Linux and macOS (Intel and Apple Silicon) based systems, on which it has been tested thoroughly.

For installation, ```cmake``` version 3.17+ is required, as well as ```libomp```.

## Library Installation Instructions

Clone the repository, and ```cd``` into the 'library' directory at the project root. Then run:

1. ```mkdir ./build```
2. ```cd ./build```
3. ```cmake ..```
4. ```sudo make install```

## Minimal Working Example Installation Instructions

After installing the library, ```cd``` into the 'example' directory at the project root. Then run:

1. ```mkdir ./build```
2. ```cd ./build```
3. ```cmake ..```
4. ```make```

## Minimal Working Example Execution Instructions

Simply ```cd``` into the directory containing the compiled executable, and run 
```./ALLLSatisfiabilitySolver [-o] [-p [n_threads]] <path/to/sat_instance.cnf>```. 

The required filepath ```path/to/sat_instance.cnf``` is the path to the CNF SAT instance to be loaded, represented in the DIMACS format.

The optional positional argument```-o``` specifies that command line output should be persisted to disk as a ```.out```
text file, along with any statistics gathered during the solve saved as a ```.csv``` file.

The optional positional argument```-p``` specifies that parallelisation should be used wherever possible.

By default, if ```n_threads``` is not specified, this will use all available (physical) threads on the system. Otherwise, ```n_threads``` are used for parallelisation.

## References

[1] R. A. Moser and G. Tardos, 'A constructive proof of the general Lovász Local Lemma' (2009), availble at: https://arxiv.org/pdf/0903.0544.pdf

[2] D. E. Knuth, 'The Art of Computer Programming, Vol. 4A' (2011), Addision-Wesley Professional.

[3] D. E. Knuth, 'The Art of Computer Programming, Vol. 4B, Fascicle 6' (2015), Addision-Wesley Professional.
