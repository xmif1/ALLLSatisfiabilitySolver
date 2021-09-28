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
typedef long long int lli;

class VariablesArray{
    public:
        ull n_vars;
        bool* vars;

        explicit VariablesArray(ull n_vars);

        bool& operator[](lli index) const;

    //private:

};

#endif //ALLLSATISFIABILITYSOLVER_VARIABLESARRAY_H
