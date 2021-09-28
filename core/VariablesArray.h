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

class VariablesArray{
    public:
        ull n_vars;
        bool* vars;

        explicit VariablesArray(ull n_vars);
};

#endif //ALLLSATISFIABILITYSOLVER_VARIABLESARRAY_H
