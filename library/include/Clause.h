//
// Created by Xandru Mifsud on 27/09/2021.
//

#ifndef ALLLSATISFIABILITYSOLVER_CLAUSE_H
#define ALLLSATISFIABILITYSOLVER_CLAUSE_H

#include <vector>

#include "VariablesArray.h"

/* Convenience class to represent the clauses of a SAT instance, along with any related meta-data. Provides a means of
 * checking whether a clause is satisfied or not for a given variable assignment. May be expanded with additional
 * functionality later on, such as partitioning or thread-safety related mechanism.
 */

template <typename tV>
class Clause{
    public:
        typedef vector<Clause<tV>*> ClauseArray;

        vector<tV>* literals;
        unsigned short int t_id{};

        explicit Clause (vector<tV>* literals, unsigned short int t_id) {
            this->literals = literals;
            this->t_id = t_id;
        }

        /* Utility function for checking whether a variable assignment of a SAT instance (represented as a VariablesArray
         * instance) satisfies or not the SAT instance. Returns true if the clause is NOT satisfied, and returns false
         * otherwise.
         */
        bool is_not_satisfied (const bool* var_arr) const {
            for(auto& l : *literals){ // For every literal in the clause...
                /* Let v = literals[i] >> 1 (i.e. v is the variable associated with the literal).
                 * If the i^th literal is encoded to an odd number (literals[i] & 1) then check if !(var_arr->vars)[v] is true,
                 * else check if (var_arr->vars)[v] is true.
                 */
                if ((l & 1) ? !(var_arr[l >> 1]) : var_arr[l >> 1]) {
                    return false; // Clause is satisfied, hence we answer 'false' to the question 'Is the clause NOT satisfied?'
                }
            }

            return true; // Clause is not satisfied, hence we answer 'true' to the question 'Is the clause NOT satisfied?'
        }
};


#endif //ALLLSATISFIABILITYSOLVER_CLAUSE_H
