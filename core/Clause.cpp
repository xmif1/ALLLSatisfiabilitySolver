//
// Created by Xandru Mifsud on 27/09/2021.
//

#include "Clause.h"

/* Constructor for a clause instance, setting the internal state using the passed parameters:
 *  i. ull n_literals is the number of literals in the clause
 * ii. ull* literals is an array of size n_literals where each element corresponds to a literal in the clause
 */
Clause::Clause(uint32_t* literals, int n_literals){
    this->literals = literals;
    this->n_literals = n_literals;
}

/* Utility function for checking whether a variable assignment of a SAT instance (represented as a VariablesArray
 * instance) satisfies or not the SAT instance. Returns true if the clause is NOT satisfied, and returns false otherwise.
 */
bool Clause::is_not_satisfied(const bool* var_arr) const{
    for(int i = 0; i < n_literals; i++){ // For every literal in the clause...
        /* Let v = literals[i] >> 1 (i.e. v is the variable associated with the literal).
         * If the i^th literal is encoded to an odd number (literals[i] & 1) then check if !(var_arr->vars)[v] is true,
         * else check if (var_arr->vars)[v] is true.
         */
        if((literals[i] & 1) ? !(var_arr[literals[i] >> 1]) : var_arr[literals[i] >> 1]){
            return false; // Clause is satisfied, hence we answer 'false' to the question 'Is the clause NOT satisfied?'
        }
    }

    return true; // Clause is not satisfied, hence we answer 'true' to the question 'Is the clause NOT satisfied?'
}

// Given another Clause instance, establishes whether the two are dependent or not. The criteria of dependency is having
// one or more variables in common between the literals of the two clauses.
bool Clause::dependent_clauses(Clause* c) const{
    // For every pair (x, y) of literals between the two clauses...
    for(int i = 0; i < n_literals; i++){
        for(int j = 0; j < c->n_literals; j++){
            // If the variable associated with x is equal to the variable associated with y, then the clauses are
            // dependent and hence we return true (no need to check further - we require AT LEAST one).
            if((literals[i] >> 1) == ((c->literals)[j] >> 1)){
                return true;
            }
        }
    }

    // Otherwise, we return false (independent if FOR ALL literal pairs, the underlying variables are distinct)
    return false;
}