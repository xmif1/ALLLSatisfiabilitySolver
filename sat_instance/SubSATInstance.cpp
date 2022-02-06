//
// Created by Xandru Mifsud on 04/10/2021.
//

#include "SubSATInstance.h"

SubSATInstance::SubSATInstance(VariablesArray* variables, vector<Clause*>* clauses, int n_threads){
    this->var_arr = variables;
    this->clauses = clauses;
    this->n_threads = n_threads;

    int rbg_ensemble_size = 1;
    if(n_threads){
        rbg_ensemble_size = n_threads;
    }

    for(int i = 0; i < rbg_ensemble_size; i++){
        auto engine = new default_random_engine(std::random_device{}()); // Get the system default random generator with a random seed
        rbg_ensemble.push_back(new RBG<default_random_engine>(*engine)); // Initialise an instance of a random boolean generator...
    }

    n_clauses = (this->clauses)->size(); // Number of clauses in SAT instance

    // Starting seed for an LCG based index to the this.clauses vector, which pseudo-randomly visits all the clauses
    // exactly once with a period of n_clauses, allowing for efficient 'shuffling' of the clauses vector; in case of
    // multi--threading, we use one for each chunk.
    if(n_threads < 2){
        C.push_back(P_9223372036854775783 % n_clauses);
    }
    else{
        for(int t = 0; t < n_threads; t++){
            C.push_back(P_9223372036854775783 % n_clauses);
        }
    }
}

/* Checks whether the SAT instance is satisfied by a given variable assignment (specified as a VariablesArray instance),
 * by checking if each Clause instance in the clauses vector is satisfied. We iterate over the clauses vector in a pseudo
 * random manner, by means of an LCG. Note that we maintain the state of the LCG as a class variable.
 *
 * In the case that a clause is not satisfied, we return a pointer to the Clause instance. Otherwise, if every clause is
 * satisfied i.e. the SAT instance is satisfied, we return a nullptr.
 */
Clause* SubSATInstance::sequential_is_satisfied(){
    for(uint32_t i = 0; i < n_clauses; i++){ // For each clause...
        // Fetch the clause specified at the (class variable) index C, and check if it is satisfied...
        if((clauses->at(C[0]))->is_not_satisfied(var_arr->vars)){ // If it is not satisfied, return a ptr to the Clause instance
            return clauses->at(C[0]);
        }
        else{ /* Otherwise, calculate the next index C by means of the LCG; note that the LCG has a period of n_clauses
               * and hence we iterate over each clause exactly once, in a non-linear fashion; this is equivalent to a
               * pseudo-random shuffling of the clauses vector, and is required for optimal convergence of the Algorithmic
               * Lovasz Local Lemma method.
               */
            C[0] = (C[0] + P_9223372036854775783) % n_clauses;
        }
    }

    return nullptr;
}

void SubSATInstance::solve(){
    if(n_threads){
        omp_set_num_threads(n_threads);
        while(true){
            auto unsat_clauses = new vector<vector<Clause*>*>;
            for(int t = 0; t < n_threads; t++){
                unsat_clauses->push_back(new vector<Clause*>);
            }

            #pragma omp parallel for schedule(dynamic) default(none) shared(unsat_clauses)
            for(uint32_t c = 0; c < n_clauses; c++){
                if((clauses->at(c))->is_not_satisfied(var_arr->vars)){
                    (unsat_clauses->at(omp_get_thread_num()))->push_back(clauses->at(c));
                }
            }

            auto indep_unsat_clauses = new vector<vector<Clause*>*>;
            for(int t = 0; t < n_threads; t++){
                indep_unsat_clauses->push_back(new vector<Clause*>);
            }

            #pragma omp parallel for schedule(static, 1) default(none) shared(unsat_clauses, indep_unsat_clauses)
            for(int t = 0; t < n_threads; t++){
                uint32_t chunk_size = unsat_clauses->at(t)->size();
                C[t] = C[t] % chunk_size;

                for(uint32_t i = 0; i < chunk_size; i++){
                    bool independent = true;
                    for(auto c: *(indep_unsat_clauses->at(t))){
                        if(((unsat_clauses->at(t))->at(C[t]))->dependent_clauses(c)){
                            independent = false;
                            break;
                        }
                    }

                    if(independent){
                        (indep_unsat_clauses->at(t))->push_back((unsat_clauses->at(t))->at(C[t]));
                    }

                    C[t] = (C[t] + P_9223372036854775783) % chunk_size;
                }

                delete unsat_clauses->at(t);
            }

            delete unsat_clauses;

            auto max_indep_unsat_clauses = parallel_k_partite_mis(indep_unsat_clauses);
            if(max_indep_unsat_clauses->empty()){
                break;
            }

            #pragma omp parallel for schedule(dynamic) default(none) shared(max_indep_unsat_clauses)
            for(uint32_t c = 0; c < max_indep_unsat_clauses->size(); c++){
                for(int i = 0; i < (max_indep_unsat_clauses->at(c))->n_literals; i++){
                    (var_arr->vars)[((max_indep_unsat_clauses->at(c))->literals)[i] >> 1] = (rbg_ensemble[omp_get_thread_num()])->sample();
                }
            }
        }
    }
    else{
        Clause* c = sequential_is_satisfied();

        while(c){ // While there exists a clause c which is not satisfied...
            // For every variable in the clause (obtained by left shifting by 1 the literal encoding), randomly re-sample
            for(int i = 0; i < c->n_literals; i++){
                (var_arr->vars)[(c->literals)[i] >> 1] = (rbg_ensemble[0])->sample();
            }

            c = sequential_is_satisfied();
        }
    }
}

bool SubSATInstance::is_ALLL_compatible() const{
    for(auto c: *clauses){
        if(floor((pow(2, c->n_literals) / exp(1.0)) - 1.0) < c->degree){
            return false;
        }
    }

    return true;
}

vector<Clause*>* SubSATInstance::bipartite_mis(vector<Clause*>* set1, vector<Clause*>* set2){
    auto result = new vector<Clause*>;

    auto idx1 = set1->begin();
    while(idx1 != set1->end()){
        auto idx2 = set2->begin();
        while(idx2 != set2->end()){
            if((*idx1)->dependent_clauses(*idx2)){
                idx2 = set2->erase(idx2);
            }
            else{
                idx2++;
            }
        }

        result->push_back(*idx1);
        idx1 = set1->erase(idx1);
    }

    for(auto c: *set2){
        result->push_back(c);
    }

    return result;
}

vector<Clause*>* SubSATInstance::parallel_k_partite_mis(vector<vector<Clause*>*>* sets){
    if(sets->size() == 1){
        return sets->at(0);
    }
    else{
        auto joined_sets = new vector<vector<Clause*>*>;
        int offset = 0;

        if(sets->size() % 2){
            offset = 1;

            auto max_idx = sets->begin();
            for(auto idx = sets->begin(); idx != sets->end(); idx++){
                if((*max_idx)->size() < (*idx)->size()){
                    max_idx = idx;
                }
            }

            joined_sets->push_back(*max_idx);
            sets->erase(max_idx);
        }

        auto n_pairs = (int) (sets->size() / 2);
        for(int t = 0; t < n_pairs; t++){
            joined_sets->push_back(nullptr);
        }

        #pragma omp parallel for schedule(static, 1) default(none) shared(n_pairs, offset, sets, joined_sets)
        for(int t = 0; t < n_pairs; t++){
            joined_sets->at(offset + t) = bipartite_mis(sets->at(2 * t), sets->at((2 * t) + 1));

            delete sets->at(2*t);
            delete sets->at((2*t) + 1);
        }

        return parallel_k_partite_mis(joined_sets);
    }
}