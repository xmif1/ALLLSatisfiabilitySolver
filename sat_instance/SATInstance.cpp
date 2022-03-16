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
             * the shift and encoding, it is mapped to 2x - 1.
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

/* Serial Algorithmic Lovasz Local Lemma of Moser and Tardos (2010), with 'shuffling' of the clauses vector to satisfy
 * the 'there exists an unsatisfied clause' condition stochasticallly.
 */
Statistics* SATInstance::sequential_solve() const{
    auto statistics = new Statistics;

    Clause* unsat_clause;

    // Get the system default random generator with a random seed
    auto engine = new default_random_engine(std::random_device{}());

    // Initialise an instance of a random boolean generator...
    auto rbg = new RBG<default_random_engine>(*engine);

    ull iterator = P_9223372036854775783 % n_clauses; // Starting seed for an LCG based index to the this.clauses vector,
                                                      // which pseudo-randomly visits all the clauses exactly once with
                                                      // a period of n_clauses, allowing for efficient 'shuffling' of the
                                                      // clauses vector;

    /* Note that by shuffling through the clauses array, we are effectively breaking any 'local minima' in the search
     * space; in between iterations, we do not get stuck successively in the same clause until it becomes satisfied, but
     * rather we visit other unsatisfied clauses and satisfy those (which in turn may influence out original unsatisfied
     * clause).
     */

    bool solved = false;
    while(!solved){ // While there exists a clause c which is not satisfied...
        bool unsat_exists = false;
        statistics->n_iterations += 1; // Update statistics

        for(uint32_t i = 0; i < n_clauses; i++){ // For each clause...
            // Fetch the clause specified at the (class variable) index C, and check if it is satisfied...

            // If it is not satisfied, return a ptr to the Clause instance
            if((clauses->at(iterator))->is_not_satisfied(var_arr->vars)){
                unsat_clause = clauses->at(iterator);
                unsat_exists = true;
                break;
            }
            else{ /* Otherwise, calculate the next index C by means of the LCG; note that the LCG has a period of
                   * n_clauses and hence we iterate over each clause exactly once, in a non-linear fashion; this is
                   * equivalent to a pseudo-random shuffling of the clauses vector, and is required for optimal
                   * convergence of the Algorithmic Lovasz Local Lemma method.
                   */
                iterator = (iterator + P_9223372036854775783) % n_clauses;
            }
        }

        if(unsat_exists){
            // For every variable in the clause (obtained by left shifting by 1 the literal encoding), randomly re-sample
            for(int i = 0; i < unsat_clause->n_literals; i++){
                (var_arr->vars)[(unsat_clause->literals)[i] >> 1] = rbg->sample();
            }

            statistics->n_resamples += unsat_clause->n_literals; // Update statistics
        }
        else{
            solved = true;
        }
    }

    return statistics;
}

/* Parallel Algorithmic Lovasz Local Lemma of Moser and Tardos (2010), with parallel unsatisfied clause checking and
 * parallel divide-and-conquer generation of non-trivial independent sets of unsatisfied clauses; the sets of independent
 * unsatisfied clauses chosen by each thread is done in a pseudo-random manner.
 */
Statistics* SATInstance::parallel_solve(){
    auto statistics = new Statistics;

    /* Each thread will be allocated a batch of clauses, divided dynamically amongst the threads. For each of these
     * batches, we check which clauses are unsatisfied and only maintain those. Hence each thread t will have a set U_t
     * of unsatisfied clauses associated with it, where for i != j then U_i and U_j are disjoint.
     *
     * For each set U_t, we must greedily choose a subset which is maximally independent I_t (ideally of largest size -
     * however this is not guaranteed and is indeed another problem in NP). To do so, we 'shuffle' U_t such that, if in
     * between iterations U_t is generated twice, the resulting independent set I_t is (probably) not. Note that by
     * shuffling through the clauses array, we are effectively breaking any 'local minima' in the search space.
     *
     * Now, to shuffle the unsatisfied clauses in U_t, we require an LCG. In particular, since our LCG depends on the
     * size of U_t and the current clause processed by a thread, we require an LCG for each thread ie. we require an
     * ensemble of LCG--based iterators, one for each thread.
     *
     * Lastly, all the maximally independent sets (MIS) I_t will be joined (through a parallel and divide-and-conquer
     * based join algorithm) into a single MIS S. The clauses in S are all independent and hence can have all their
     * variables re-sampled. For maximum utilisation, variable resampling is carried out in parallel as well, by dynamic
     * allocation of clauses in S to all the available threads.
     *
     * Note that the conditions on S we consider are more relaxed than those required by Moser and Tardos (2010); namely
     * S is an independent set in the set of unsatisfied clauses but not a MIS, however S must be a MIS in the union of
     * all the I_ts.
     *
     * Since re-sampling is stochastically carried out using our random boolean generator (RBG) implementation, and an
     * RBG instance maintains state, we require an RBG for each thread - hence we must also initialise one for each.
     *
     * SUMMARY: 1. Split clauses in n_threads batches {C_j}.
     *          2. For each C_j, in parallel find the subset U_j of C_j where U_j is all the unsatisfied clauses.
     *          3. For each U_j, in parallel greedily find a MIS I_j in U_j, using shuffling to break any local minima.
     *          4. Join the family of MISs {I_j} into a single MIS S (in parallel using divide-and-conquer).
     *          5. For each clause C in S, in parallel resample the variables in C using an RBG.
     */

    vector<RBG<default_random_engine>*> rbg_ensemble;
    vector<ull> iterators;

    // Initialise ensembles of RBGs and iterators...
    for(int i = 0; i < n_threads; i++){
        // Get the system default random generator with a random seed
        auto engine = new default_random_engine(std::random_device{}());

        // Initialise an instance of a random boolean generator
        rbg_ensemble.push_back(new RBG<default_random_engine>(*engine));

        // Initialise LCG-based iterators
        iterators.push_back(P_9223372036854775783 % n_clauses);

        // Initialise statistics
        statistics->n_thread_resamples.push_back(0);
    }

    omp_set_num_threads(n_threads);
    while(true){ // While there exists a clause c which is not satisfied...

        // Initialise an empty vector U_t for the unsatisfied clauses found by each thread
        auto unsat_clauses = new vector<ClausesArray*>;
        for(int t = 0; t < n_threads; t++){
            unsat_clauses->push_back(new ClausesArray);
        }

        statistics->n_iterations += 1; // Update statistics

        // Split clauses in n_threads batches {C_j}, in parallel through dynamic allocation for maximum utilisation
        #pragma omp parallel for schedule(dynamic) default(none) shared(unsat_clauses)
        for(uint32_t c = 0; c < n_clauses; c++){
            if((clauses->at(c))->is_not_satisfied(var_arr->vars)){
                (unsat_clauses->at(omp_get_thread_num()))->push_back(clauses->at(c));
            }
        }

        // Initialise an empty vector I_t for the MIS of unsatisfied clauses associated with U_t of each thread
        auto indep_unsat_clauses = new vector<ClausesArray*>;
        for(int t = 0; t < n_threads; t++){
            indep_unsat_clauses->push_back(new ClausesArray);
        }

        // Find a subset MIS I_t in each set of unsatisfied clauses U_t, in parallel (statically, one for each thread);
        // As described earlier, we visit the clauses in each U_t in a pseudo-random order (effectively shuffling)
        #pragma omp parallel for schedule(static, 1) default(none) shared(unsat_clauses, indep_unsat_clauses, iterators)
        for(int t = 0; t < n_threads; t++){
            // Initialise iterator LCG such that it does not exceed chunk_size (ie number of clauses in U_t)
            uint32_t chunk_size = unsat_clauses->at(t)->size();
            iterators[t] = iterators[t] % chunk_size;

            for(uint32_t i = 0; i < chunk_size; i++){ // Visit each clause in U_t; note that the iterator i is not used
                // as an index; rather the LCG corresponding to the thread is. Since the LCG has period chunk_size, each
                // time we update it with every iteration i, we ensure that no clause is visited twice. Hence since i
                // has a range from 0 to chunk_size - 1, we visit all the clauses in U_t in a non-linear fashion.

                bool independent = true; // If current clause in U_t is dependent on some clause in the current state of
                                         // the independent set I_t, then change flag to false.
                for(auto c: *(indep_unsat_clauses->at(t))){ // Check if dependent on any clause in I_t...
                    if(((unsat_clauses->at(t))->at(iterators[t]))->dependent_clauses(c)){
                        // If dependent, update independent flag to false and break (no need to continue checking -
                        // the current clause in U_t will not be added to the independent set)
                        independent = false;
                        break;
                    }
                }

                if(independent){ // If independent, then add to independent set I_t
                    (indep_unsat_clauses->at(t))->push_back((unsat_clauses->at(t))->at(iterators[t]));
                }

                iterators[t] = (iterators[t] + P_9223372036854775783) % chunk_size; // Update LCG
            } // At the end, I_t will be a MIS since all clauses in U_t have been exhausted

            delete unsat_clauses->at(t); // Memory management; we no longer require U_t, only I_t
        }

        delete unsat_clauses; // Memory management

        // Join the n_thread MISs {I_t} into a single MIS (through the greedy_parallel_mis_join() function)
        auto max_indep_unsat_clauses = greedy_parallel_mis_join(indep_unsat_clauses);

        // If max_indep_unsat_clauses is empty, then no unsat clauses have been found across all threads - SATISFIED
        if(max_indep_unsat_clauses->empty()){
            delete max_indep_unsat_clauses; // Memory management
            break;
        }

        statistics->avg_mis_size += max_indep_unsat_clauses->size(); // Update statistics

        // In parallel and dynamically, re-sample the variables appearing in the clauses of the MIS constructed
        #pragma omp parallel for schedule(dynamic) default(none) shared(max_indep_unsat_clauses, rbg_ensemble, statistics)
        for(uint32_t c = 0; c < max_indep_unsat_clauses->size(); c++){
            int tid = omp_get_thread_num(); // Get thread identifier

            // Re--sample every variable in the clause using the RBG associated with the thread tid
            for(int i = 0; i < (max_indep_unsat_clauses->at(c))->n_literals; i++){
                (var_arr->vars)[((max_indep_unsat_clauses->at(c))->literals)[i] >> 1] = (rbg_ensemble[tid])->sample();
            }

            // Update statistics
            statistics->n_thread_resamples.at(tid) += (max_indep_unsat_clauses->at(c))->n_literals;
        }

        delete max_indep_unsat_clauses; // Memory management
    }

    for(int t = 0; t < n_threads; t++){
        statistics->n_resamples += (statistics->n_thread_resamples).at(t); // Update statistics
    }
    statistics->avg_mis_size = (ull) (statistics->avg_mis_size / statistics->n_iterations); // Update statistics

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
ClausesArray* SATInstance::greedy_mis_join(ClausesArray* set1, ClausesArray* set2){
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

// Utility function for recursively and in parallel joining k disjoint maximally independent sets into a single one
ClausesArray* SATInstance::greedy_parallel_mis_join(vector<ClausesArray*>* sets){
    // BASE CASE
    if(sets->size() == 1){ // If only a single MIS is passed as input, then simply return it
        auto ret_set = sets->at(0);
        delete sets;

        return ret_set;
    }
    else{ // RECURSIVE CASE
        /* The join works by pairing the independent sets passed as input, and (in parallel) each pair is joined into a
         * single MIS using the greedy_mis_join() algorithm. Note that in the case of an odd number of sets, one is not
         * paired.
         *
         * The greedy_parallel_mis_join is then called once again on the resulting sets, until the base case is satisfied.
         */

        auto joined_sets = new vector<ClausesArray*>; // Resulting MIS from pair--wise joins
        int offset = 0; // Index offset from which to start populating joined_sets (if odd then index 0 is populated
        // with the largest input MIS, and the remaining are paired off and joined_sets is populated with the joined
        // sets from index 1 onwards ie. offset is set to 1.

        if(sets->size() % 2){ // If odd number of initial sets...
            offset = 1; // Set offset to 1

            // Find the input MIS with the largest size
            auto max_idx = sets->begin();
            for(auto idx = sets->begin(); idx != sets->end(); idx++){
                if((*max_idx)->size() < (*idx)->size()){
                    max_idx = idx;
                }
            }

            joined_sets->push_back(*max_idx); // Append largest MIS to joined_sets
            sets->erase(max_idx); // Memory management
        } // The remaining input MISs are paired off and joined into a single MIS

        auto n_pairs = (int) (sets->size() / 2); // Number of sets resulting from pairing
        for(int t = 0; t < n_pairs; t++){ // Initialisation
            joined_sets->push_back(nullptr);
        }

        // In parallel join each pair of input MISs using the greedy_mis_join() algorithm (one for each thread; number of
        // pairs does not exceed n_threads)
        #pragma omp parallel for schedule(static, 1) default(none) shared(n_pairs, offset, sets, joined_sets)
        for(int t = 0; t < n_pairs; t++){
            // Join pairs greedily into a single independent set...
            if(sets->at((2 * t) + 1)->size() <= sets->at(2 * t)->size()) {
                joined_sets->at(offset + t) = greedy_mis_join(sets->at(2 * t), sets->at((2 * t) + 1));
            }
            else{
                joined_sets->at(offset + t) = greedy_mis_join(sets->at((2 * t) + 1), sets->at(2 * t));
            }

            // Memory management (we only require the resulting joined set from two initial sets)
            delete sets->at(2*t);
            delete sets->at((2*t) + 1);
        }

        delete sets; // Memory management

        return greedy_parallel_mis_join(joined_sets); // Recursive join of the resulting MISs
    }
}