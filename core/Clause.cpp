//
// Created by Xandru Mifsud on 27/09/2021.
//

#include "Clause.h"

/* Constructor for a clause instance, setting the internal state using the passed parameters:
 *  i. ull n_literals is the number of literals in the clause
 * ii. ull* literals is an array of size n_literals where each element corresponds to a literal in the clause
 */
Clause::Clause(ull* literals, ull n_literals){
    this->literals = literals;
    this->n_literals = n_literals;
}

/* Utility function for checking whether a variable assignment of a SAT instance (represented as a VariablesArray
 * instance) satisfies or not the SAT instance. Returns true if the clause is NOT satisfied, and returns false otherwise.
 */
bool Clause::is_not_satisfied(VariablesArray* var_arr) const{
    for(ull i = 0; i < n_literals; i++){ // For every literal in the clause...
        /* Let v = literals[i] >> 1 (i.e. v is the variable associated with the literal).
         * If the i^th literal is encoded to an odd number (literals[i] & 1) then check if !(var_arr->vars)[v] is true,
         * else check if (var_arr->vars)[v] is true.
         */
        if((literals[i] & 1) ? !((var_arr->vars)[literals[i] >> 1]) : (var_arr->vars)[literals[i] >> 1]){
            return false; // Clause is satisfied, hence we answer 'false' to the question 'Is the clause NOT satisfied?'
        }
    }

    return true; // Clause is not satisfied, hence we answer 'true' to the question 'Is the clause NOT satisfied?'
}