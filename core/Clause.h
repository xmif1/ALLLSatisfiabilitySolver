//
// Created by Xandru Mifsud on 27/09/2021.
//

#ifndef ALLLSATISFIABILITYSOLVER_CLAUSE_H
#define ALLLSATISFIABILITYSOLVER_CLAUSE_H

#include <vector>

#include "VariablesArray.h"

/* Convenience class to represent the clauses of a SAT instance, along with any related meta data. Provides a means of
 * checking whether a clause is satisfied or not for a given variable assignment. May be expanded with additional
 * functionality later on, such as partitioning or thread-safety related mechanism.
 */
class Clause{
    public:
        ull* literals;
        ull n_literals;

        Clause(ull* literals, ull n_literals);

        bool is_not_satisfied(VariablesArray* var_arr) const;
};


#endif //ALLLSATISFIABILITYSOLVER_CLAUSE_H
