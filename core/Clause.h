//
// Created by Xandru Mifsud on 27/09/2021.
//

#ifndef ALLLSATISFIABILITYSOLVER_CLAUSE_H
#define ALLLSATISFIABILITYSOLVER_CLAUSE_H

#include <vector>

#include "VariablesArray.h"

class Clause{
    public:
        lli* literals;
        ull n_literals;

        Clause(lli* literals, ull n_literals);

        bool is_not_satisfied(VariablesArray* var_arr) const;
};


#endif //ALLLSATISFIABILITYSOLVER_CLAUSE_H
