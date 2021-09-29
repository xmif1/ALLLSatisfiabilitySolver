//
// Created by Xandru Mifsud on 27/09/2021.
//

#ifndef ALLLSATISFIABILITYSOLVER_VARIABLESARRAY_H
#define ALLLSATISFIABILITYSOLVER_VARIABLESARRAY_H

#include <iostream>
#include <random>

#include "RandomBoolGenerator.h"

using namespace std;

typedef unsigned long long ull;

// Convenience class to represent the variables of a SAT instance, along with any related meta data. May be expanded
// with additional functionality later on, such as partitioning or thread-safety related mechanism.
class VariablesArray{
    public:
        ull n_vars;
        bool* vars;

        explicit VariablesArray(ull n_vars);
};

#endif //ALLLSATISFIABILITYSOLVER_VARIABLESARRAY_H
