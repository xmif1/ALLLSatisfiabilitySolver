//
// Created by Xandru Mifsud on 27/09/2021.
//

#ifndef ALLLSATISFIABILITYSOLVER_VARIABLESARRAY_H
#define ALLLSATISFIABILITYSOLVER_VARIABLESARRAY_H

#include <iostream>
#include <random>

#include "RandomBoolGenerator.h"

using namespace std;

// Convenience class to represent the variables of a SAT instance, along with any related meta data. May be expanded
// with additional functionality later on, such as partitioning or thread-safety related mechanism.
template <typename tV>
class VariablesArray{
    public:
        tV n_vars;
        bool* vars;

        explicit VariablesArray (tV n_vars) {
            default_random_engine engine(std::random_device{}()); // Get the system default random generator with a random seed
            RBG<default_random_engine> rbg(engine); // Initialise an instance of a random boolean generator...

            this->n_vars = n_vars;
            vars = new bool[n_vars];

            // Populate each entry in the vars array with a randomly generated boolean value...
            for (tV i = 0; i < n_vars; i++) {
                vars[i] = rbg.sample();
            }
        }
};

#endif //ALLLSATISFIABILITYSOLVER_VARIABLESARRAY_H
