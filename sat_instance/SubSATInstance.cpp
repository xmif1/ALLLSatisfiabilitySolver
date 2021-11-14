//
// Created by Xandru Mifsud on 04/10/2021.
//

#include "SubSATInstance.h"

SubSATInstance::SubSATInstance(VariablesArray* variables, vector<Clause*>* clauses, ull parallel_resample){
    this->var_arr = variables;
    this->clauses = clauses;
    this->parallel_resample = parallel_resample;

    ull max_clause_literals = 0;
    for(auto c: *clauses){
        if(c->n_literals > max_clause_literals){
            max_clause_literals = c->n_literals;
        }
    }

    for(ull i = 0; i < max_clause_literals; i++){
        auto engine = new default_random_engine(std::random_device{}()); // Get the system default random generator with a random seed
        rbg_ensemble.push_back(new RBG<default_random_engine>(*engine)); // Initialise an instance of a random boolean generator...
    }

    n_clauses = (this->clauses)->size(); // Number of clauses in SAT instance
    C = P_2e64_m59 % n_clauses; // Starting seed for an LCG based index to the this.clauses vector, which pseudo-randomly
                                // visits all the clauses exactly once with a period of n_clauses, allowing for efficient
                                // 'shuffling' of the clauses vector.
}

/* Checks whether the SAT instance is satisfied by a given variable assignment (specified as a VariablesArray instance),
 * by checking if each Clause instance in the clauses vector is satisfied. We iterate over the clauses vector in a pseudo
 * random manner, by means of an LCG. Note that we maintain the state of the LCG as a class variable.
 *
 * In the case that a clause is not satisfied, we return a pointer to the Clause instance. Otherwise, if every clause is
 * satisfied i.e. the SAT instance is satisfied, we return a nullptr.
 */
Clause* SubSATInstance::is_satisfied(){
    for(ull i = 0; i < n_clauses; i++){ // For each clause...
        // Fetch the clause specified at the (class variable) index C, and check if it is satisfied...
        if((clauses->at(C))->is_not_satisfied(var_arr)){ // If it is not satisfied, return a ptr to the Clause instance
            return clauses->at(C);
        }
        else{ /* Otherwise, calculate the next index C by means of the LCG; note that the LCG has a period of n_clauses
               * and hence we iterate over each clause exactly once, in a non-linear fashion; this is equivalent to a
               * pseudo-random shuffling of the clauses vector, and is required for optimal convergence of the Algorithmic
               * Lovasz Local Lemma method.
               */
            C = (C + P_2e64_m59) % n_clauses;
        }
    }

    return nullptr;
}

void SubSATInstance::solve(){
    Clause* c = is_satisfied();
    if(parallel_resample){
        while(c){ // While there exists a clause c which is not satisfied...
            // For every variable in the clause (obtained by left shifting by 1 the literal encoding), randomly re-sample
            #pragma omp parallel for default(none) if (parallel_resample <= c->n_literals) shared(c)
            for(ull i = 0; i < c->n_literals; i++){
                (var_arr->vars)[(c->literals)[i] >> 1] = (rbg_ensemble.at(i))->sample();
            }

            c = is_satisfied();
        }
    }
    else{
        while(c){ // While there exists a clause c which is not satisfied...
            // For every variable in the clause (obtained by left shifting by 1 the literal encoding), randomly re-sample
            for(ull i = 0; i < c->n_literals; i++){
                (var_arr->vars)[(c->literals)[i] >> 1] = (rbg_ensemble.at(i))->sample();
            }

            c = is_satisfied();
        }
    }
}

bool SubSATInstance::is_ALLL_compatible() const{
    for(auto c: *clauses){
        if(floor(pow(2, c->n_literals) - 1) < c->degree){
            return false;
        }
    }

    return true;
}