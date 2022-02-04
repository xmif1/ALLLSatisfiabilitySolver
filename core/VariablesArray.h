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
class VariablesArray{
    public:
        uint32_t n_vars;
        bool* vars;

        explicit VariablesArray(uint32_t n_vars);
};

#endif //ALLLSATISFIABILITYSOLVER_VARIABLESARRAY_H
