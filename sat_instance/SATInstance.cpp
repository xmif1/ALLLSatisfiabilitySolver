//
// Created by Xandru Mifsud on 28/09/2021.
//

#include "SATInstance.h"

// =====================================================================================================================
// --------------------------------------------------- INITIALISATION --------------------------------------------------
// =====================================================================================================================

/* Constructor for a SATInstance object, with the file path to a DIMACS formatted .cnf file accepted as input.
 * This constructor is responsible for:
 *  i. Loading meta-data (number of variables, clauses etc) as well as reading the clauses from the specified .cnf file,
 *     using the CNF_IO package (available at https://people.sc.fsu.edu/~jburkardt/cpp_src/cnf_io/cnf_io.html)
 * ii. Encoding literals such that a variable x, where x is a non-negative integer, is mapped to 2x while its negation
 *     is mapped to 2x + 1. For this encoding, we can obtain the variable associated with a literal by a single left
 *     shift i.e. the variable v associated with a literal l is: v = l >> 1.
 */
SATInstance::SATInstance(const string& cnf_file_name, int n_threads){
    int c_num, v_num, l_num;
    int* l_c_num;
    int* l_val;

    // Read the meta-data from the .cnf file
    if(cnf_header_read(cnf_file_name, &v_num, &c_num, &l_num)){
        cout << "The header information could not be read. Exiting..." << endl;
        exit(1);
    }

    l_c_num = new int[c_num];
    l_val = new int[l_num];

    // Read the clause data from the .cnf file
    cnf_data_read(cnf_file_name, v_num, c_num, l_num, l_c_num, l_val);

    int c, l, l_c;

    l = 0;
    for(c = 0; c < c_num; c++){ // For each read clauses, create a new Clause instance with its literals encoded in the
                                // aforementioned manner
        // l_c_num[c] = # of signed literals in clause c
        uint32_t* literals; literals = new uint32_t[l_c_num[c]];

        for(l_c = 0; l_c < l_c_num[c]; l_c++){ // For every literal in the clause...
            /* Note that the DIMACS format uses 0 as a special character, hence the variables are labelled as positive
             * (and hence non-zero) integers; Since we wish to use the variable x as an index to an array (where array
             * index starts from zero), we wish to shift the DIMACS variable label x to x - 1. Hence the DIMACS variable
             * x is encoded as 2x - 2. The negation of a variable x in DIMACS is the negative integer -x, hence applying
             * the shift and encoding, it is mapped to -2x - 1.
             */
            literals[l_c] = (0 < l_val[l]) ? (2 * l_val[l]) - 2 : ((-2) * l_val[l]) - 1;

            l += 1;
        }

        clauses->push_back(new Clause(literals, l_c_num[c])); // Populate clauses array
    }

    // Initialise internal state variables...
    n_vars = v_num;                                         // Number of variables in SAT instance
    n_clauses = c_num;                                      // Number of clauses in SAT instance
    n_literals = l_num;                                     // Number of literals in SAT instance
    this->n_threads = n_threads;                            // Number of threads for solving (if using parallel solver)
    var_arr = new VariablesArray(n_vars);                   // Encoded variables array
}

// =====================================================================================================================
// --------------------------------------------------- ALLL SOLVERS ----------------------------------------------------
// =====================================================================================================================

// SAT solver based on the Algorithmic Lovasz Local Lemma of Moser and Tardos (2010); wrapper to either sequential or
// parallel solver depending on whether 1 or more threads specified during instantiation.
Statistics* SATInstance::solve(){
    if(n_threads){
        return parallel_solve();
    }
    else{
        return sequential_solve();
    }
}

Statistics* SATInstance::sequential_solve() const{
    auto statistics = new Statistics;

    Clause* unsat_clause;

    auto engine = new default_random_engine(std::random_device{}()); // Get the system default random generator with a random seed
    auto rbg = new RBG<default_random_engine>(*engine); // Initialise an instance of a random boolean generator...
    ull iterator = P_9223372036854775783 % n_clauses; // Starting seed for an LCG based index to the this.clauses vector,
                                                      // which pseudo-randomly visits all the clauses exactly once with
                                                      // a period of n_clauses, allowing for efficient 'shuffling' of the
                                                      // clauses vector;
    bool solved = false;
    while(!solved){ // While there exists a clause c which is not satisfied...
        bool unsat_exists = false;
        statistics->n_iterations += 1;

        for(uint32_t i = 0; i < n_clauses; i++){ // For each clause...
            // Fetch the clause specified at the (class variable) index C, and check if it is satisfied...
            if((clauses->at(iterator))->is_not_satisfied(var_arr->vars)){ // If it is not satisfied, return a ptr to the Clause instance
                unsat_clause = clauses->at(iterator);
                unsat_exists = true;
                break;
            }
            else{ /* Otherwise, calculate the next index C by means of the LCG; note that the LCG has a period of n_clauses
               * and hence we iterate over each clause exactly once, in a non-linear fashion; this is equivalent to a
               * pseudo-random shuffling of the clauses vector, and is required for optimal convergence of the Algorithmic
               * Lovasz Local Lemma method.
               */
                iterator = (iterator + P_9223372036854775783) % n_clauses;
            }
        }

        if(unsat_exists){
            // For every variable in the clause (obtained by left shifting by 1 the literal encoding), randomly re-sample
            for(int i = 0; i < unsat_clause->n_literals; i++){
                (var_arr->vars)[(unsat_clause->literals)[i] >> 1] = rbg->sample();
            }

            statistics->n_resamples += unsat_clause->n_literals;
        }
        else{
            solved = true;
        }
    }

    return statistics;
}

Statistics* SATInstance::parallel_solve(){
    auto statistics = new Statistics;

    vector<RBG<default_random_engine>*> rbg_ensemble;
    vector<ull> iterators;

    for(int i = 0; i < n_threads; i++){
        auto engine = new default_random_engine(std::random_device{}()); // Get the system default random generator with a random seed
        rbg_ensemble.push_back(new RBG<default_random_engine>(*engine)); // Initialise an instance of a random boolean generator...
        iterators.push_back(P_9223372036854775783 % n_clauses);
        statistics->n_thread_resamples.push_back(0);
    }

    omp_set_num_threads(n_threads);
    while(true){
        auto unsat_clauses = new vector<ClausesArray*>;
        for(int t = 0; t < n_threads; t++){
            unsat_clauses->push_back(new ClausesArray);
        }

        statistics->n_iterations += 1;

        #pragma omp parallel for schedule(dynamic) default(none) shared(unsat_clauses)
        for(uint32_t c = 0; c < n_clauses; c++){
            if((clauses->at(c))->is_not_satisfied(var_arr->vars)){
                (unsat_clauses->at(omp_get_thread_num()))->push_back(clauses->at(c));
            }
        }

        auto indep_unsat_clauses = new vector<ClausesArray*>;
        for(int t = 0; t < n_threads; t++){
            indep_unsat_clauses->push_back(new ClausesArray);
        }

        #pragma omp parallel for schedule(static, 1) default(none) shared(unsat_clauses, indep_unsat_clauses, iterators)
        for(int t = 0; t < n_threads; t++){
            uint32_t chunk_size = unsat_clauses->at(t)->size();
            iterators[t] = iterators[t] % chunk_size;

            for(uint32_t i = 0; i < chunk_size; i++){
                bool independent = true;
                for(auto c: *(indep_unsat_clauses->at(t))){
                    if(((unsat_clauses->at(t))->at(iterators[t]))->dependent_clauses(c)){
                        independent = false;
                        break;
                    }
                }

                if(independent){
                    (indep_unsat_clauses->at(t))->push_back((unsat_clauses->at(t))->at(iterators[t]));
                }

                iterators[t] = (iterators[t] + P_9223372036854775783) % chunk_size;
            }

            delete unsat_clauses->at(t);
        }

        delete unsat_clauses;

        auto max_indep_unsat_clauses = parallel_k_partite_mis(indep_unsat_clauses);
        if(max_indep_unsat_clauses->empty()){
            delete max_indep_unsat_clauses;
            break;
        }

        statistics->avg_mis_size += max_indep_unsat_clauses->size();

        #pragma omp parallel for schedule(dynamic) default(none) shared(max_indep_unsat_clauses, rbg_ensemble, statistics)
        for(uint32_t c = 0; c < max_indep_unsat_clauses->size(); c++){
            int tid = omp_get_thread_num();

            for(int i = 0; i < (max_indep_unsat_clauses->at(c))->n_literals; i++){
                (var_arr->vars)[((max_indep_unsat_clauses->at(c))->literals)[i] >> 1] = (rbg_ensemble[tid])->sample();
            }

            statistics->n_thread_resamples.at(tid) += (max_indep_unsat_clauses->at(c))->n_literals;
        }

        delete max_indep_unsat_clauses;
    }

    for(int t = 0; t < n_threads; t++){
        statistics->n_resamples += (statistics->n_thread_resamples).at(t);
    }
    statistics->avg_mis_size = (ull) (statistics->avg_mis_size / statistics->n_iterations);

    return statistics;
}

// =====================================================================================================================
// ----------------------------------------------------- UTILITIES -----------------------------------------------------
// =====================================================================================================================

// Convenience function for checking whether the assignments in var_arr represent a valid solution or not.
bool SATInstance::verify_validity() const{
    for(auto c: *clauses){
        if(c->is_not_satisfied(var_arr->vars)){
            return false;
        }
    }

    return true;
}

// Utility function for constructing a maximally independent set (MIS) from two other such sets
ClausesArray* SATInstance::bipartite_mis(ClausesArray* set1, ClausesArray* set2){
    auto result = new ClausesArray; // The resulting MIS from joining set1 and set2

    /* We begin by iterating across all the elements in set1 and checking it against every element in set2. If for some
     * element x in set1 there is an element y in set2 such that x and y are independent, then we delete y from set2 (on
     * the fly). Hence, the resulting MIS will always contain set1 as a subset (where set1 is assumed to be a MIS).
     */
    auto idx1 = set1->begin();
    while(idx1 != set1->end()){ // Consider an element x in set1...
        auto idx2 = set2->begin();
        while(idx2 != set2->end()){ // for every element y in set2
            if((*idx1)->dependent_clauses(*idx2)){ // if y and x are dependent, remove y from set2
                idx2 = set2->erase(idx2);
            }
            else{ // else maintain y in set2 and check the next element y' (if any) in set2
                idx2++;
            }
        }

        result->push_back(*idx1); // In any case, add x from set1 to the resulting MIS
        idx1 = set1->erase(idx1); // Erase x from set1 (simply for memory management)
    }

    /* As a result, all the remaining elements in set2 are independent from all the elements in set1, i.e. from all the
     * elements in the resulting MIS thus far. Hence we simply take the union of the two sets (the current resulting set
     * and the augmented set2).
     */
    for(auto c: *set2){ // Added every element in the augmented set2 to the resulting MIS
        result->push_back(c);
    }

    return result;
}

ClausesArray* SATInstance::parallel_k_partite_mis(vector<ClausesArray*>* sets){
    if(sets->size() == 1){
        auto ret_set = sets->at(0);
        delete sets;

        return ret_set;
    }
    else{
        auto joined_sets = new vector<ClausesArray*>;
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

        delete sets;

        return parallel_k_partite_mis(joined_sets);
    }
}