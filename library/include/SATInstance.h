//
// Created by Xandru Mifsud on 28/09/2021.
//

#ifndef ALLLSATISFIABILITYSOLVER_SATINSTANCE_H
#define ALLLSATISFIABILITYSOLVER_SATINSTANCE_H

#include <unordered_map>
#include <iostream>
#include <utility>
#include <vector>
#include <random>
#include <set>

#include <omp.h>
#include <google/dense_hash_set>
#include <boost/functional/hash.hpp>

#include "RandomBoolGenerator.h"
#include "VariablesArray.h"
#include "Clause.h"

using namespace std;

typedef struct Statistics{ // Structure for maintaining some simple statistics on the instance solve, such as the number
                           // of iterations, variable resamples (on a per-thread basis) and the average maximal indep.
                           // set size (per-thread, in the parallel ALLL algorithm case)
    ull n_iterations = 0;
    ull n_resamples = 0;
    ull avg_mis_size = 0;
    vector<ull> n_thread_resamples;
} Statistics;

/* Represents a CNF SAT instance loaded from a DIMACS-formatted file, providing a number of functions to:
 *  i. Solve, by means of the Algorithmic Lovasz Local Lemma (Moser & Tardos, 2010), the SAT instance, in either a
 *     sequential or parallel manner.
 * ii. Check whether a variable assignment satisfies the SAT instance.
 */
template <typename tV>
class SATInstance{
    typedef vector<Clause<uint32_t>*> ClausesArray;
    typedef pair<Clause<uint32_t>*, Clause<uint32_t>*> ClausePair;
    typedef google::dense_hash_set<Clause<uint32_t>*, boost::hash<Clause<uint32_t>*>> ClauseHashSet;
    typedef google::dense_hash_set<ClausePair, boost::hash<ClausePair>> ClauseCache;

    public:
        tV n_vars;
        ull n_clauses = 0;

        VariablesArray<tV>* var_arr;
        vector<ClauseHashSet*>* clauses{};

        // Constructor for a SATInstance object
        SATInstance(VariablesArray<tV>* var_arr, vector<ClauseHashSet*>* clauses, int n_threads, int cache_depth){
            this->var_arr = var_arr;                        // Encoded variables array
            this->clauses = clauses;                        // Clauses composing the SAT instance
            this->n_threads = n_threads;                    // Number of threads for solving (if using parallel solver)
            this->cache_depth = cache_depth;                // Size of clause dependency hash map

            n_vars = var_arr->n_vars;                       // Number of variables in SAT instance

            for(auto c : *clauses) {    // Number of clauses in SAT instance
                n_clauses += c->size();
            }
        }

        // SAT solver based on the Algorithmic Lovasz Local Lemma of Moser and Tardos (2010); wrapper to either sequential
        // or parallel solver depending on whether 1 or more threads specified during instantiation.
        Statistics* solve(){
            return parallel_solve();
//            if(n_threads){
//                return parallel_solve();
//            }
//            else{
//                return sequential_solve();
//            }
        }

        // Convenience function for checking whether the assignments in var_arr represent a valid solution or not.
        bool verify_validity() const{
            for(auto clauseArr : *clauses){
                for(auto clause = clauseArr->begin(); clause != clauseArr->end(); clause++){
                    if((*clause)->is_not_satisfied(var_arr->vars)) {
                        return false;
                    }
                }
            }

            return true;
        }

    private:
        // Prime number for LCG over the clauses array (which should be reasonably large enough...we hope...)
        int n_threads{};
        int cache_depth{};

        // =============================================================================================================
        // ----------------------------------------------- ALLL SOLVERS ------------------------------------------------
        // =============================================================================================================

        /* Parallel Algorithmic Lovasz Local Lemma of Moser and Tardos (2010), with parallel unsatisfied clause checking
         * and parallel divide-and-conquer generation of non-trivial independent sets of unsatisfied clauses; the sets
         * of independent unsatisfied clauses chosen by each thread is done in a pseudo-random manner.
         */
        Statistics* parallel_solve(){
            auto statistics = new Statistics;

            /* Each thread will be allocated a batch of clauses, divided dynamically amongst the threads. For each of
             * these batches, we check which clauses are unsatisfied and only maintain those. Hence, each thread t will
             * have a set U_t of unsatisfied clauses associated with it, where for i != j then U_i and U_j are disjoint.
             *
             * For each set U_t, we must greedily choose a subset which is maximally independent I_t (ideally of the
             * largest size - however this is not guaranteed and is indeed another problem in NP). To do so, we 'shuffle'
             * U_t such that, if in between iterations U_t is generated twice, the resulting independent set I_t is
             * (probably) not. Note that by shuffling through the clauses array, we are effectively breaking any 'local
             * minima' in the search space.
             *
             * Now, to shuffle the unsatisfied clauses in U_t, we require an LCG. In particular, since our LCG depends
             * on the size of U_t and the current clause processed by a thread, we require an LCG for each thread i.e.
             * we require an ensemble of LCG--based clause_iterators, one for each thread.
             *
             * Lastly, all the maximally independent sets (MIS) I_t will be joined (through a parallel and divide-and-
             * conquer based join algorithm) into a single MIS S. The clauses in S are all independent and hence can
             * have all their variables re-sampled. For maximum utilisation, variable resampling is carried out in
             * parallel as well, by dynamic allocation of clauses in S to all the available threads.
             *
             * Note that the conditions on S we consider are more relaxed than those required by Moser and Tardos (2010);
             * namely S is an independent set in the set of unsatisfied clauses but not a MIS, however S must be a MIS
             * in the union of all the I_ts.
             *
             * Since re-sampling is stochastically carried out using our random boolean generator (RBG) implementation,
             * and an RBG instance maintains state, we require an RBG for each thread - hence we must also initialise
             * one for each.
             *
             * SUMMARY: 1. Split clauses in n_threads batches {C_j}.
             *          2. For each C_j, in parallel find the subset U_j of C_j where U_j is all the unsatisfied clauses.
             *          3. For each U_j, in parallel greedily find a MIS I_j in U_j, using shuffling to break any local minima.
             *          4. Join the family of MISs {I_j} into a single MIS S (in parallel using divide-and-conquer).
             *          5. For each clause C in S, in parallel resample the variables in C using an RBG.
             */

            vector<ClauseCache*> caches;
            vector<RBG<default_random_engine>*> rbg_ensemble;

            // Initialise ensembles of RBGs and clause_iterators...
            for(int i = 0; i < n_threads; i++){
                // Get the system default random generator with a random seed
                auto engine = new default_random_engine(std::random_device{}());

                // Initialise an instance of a random boolean generator
                rbg_ensemble.push_back(new RBG<default_random_engine>(*engine));

                // Initialise statistics
                statistics->n_thread_resamples.push_back(0);

                // Initialise caching structures
                caches.push_back(new ClauseCache());
                caches.at(i)->set_empty_key({new Clause<uint32_t>(nullptr, i), new Clause<uint32_t>(nullptr, i)});
                caches.at(i)->set_deleted_key({new Clause<uint32_t>(nullptr, i), new Clause<uint32_t>(nullptr, i)});
            }

            auto unsat_clauses = new vector<ClauseHashSet*>;

            // Initialise an empty vector U_t for the unsatisfied clauses found by each thread
            for(int t = 0; t < n_threads; t++){
                unsat_clauses->push_back(new ClauseHashSet);
                unsat_clauses->at(t)->set_empty_key(new Clause<uint32_t>(nullptr, t));
                unsat_clauses->at(t)->set_deleted_key(new Clause<uint32_t>(nullptr, t));
            }

            omp_set_num_threads(n_threads);
            while(true){ // While there exists a clause c which is not satisfied...
                statistics->n_iterations += 1; // Update statistics

                bool allEmpty = true;
                for(auto unsatClauses : *unsat_clauses){
                    if(!unsatClauses->empty()){
                        allEmpty = false;
                        break;
                    }
                }

                if(allEmpty){
                    // Split clauses in n_threads batches {C_j}, in parallel through dynamic allocation
                    #pragma omp parallel for schedule(static, 1) default(none) shared(unsat_clauses)
                    for(int t = 0; t < n_threads; t++){
                        for(ClauseHashSet::iterator clause = clauses->at(t)->begin(); clause != clauses->at(t)->end(); ++clause){
                            if((*clause)->is_not_satisfied(var_arr->vars)){
                                (*unsat_clauses->at(t)).insert(*clause);
                            }
                        }
                    }

                    bool solved = true;
                    for(auto unsatClauses : *unsat_clauses){
                        if(!unsatClauses->empty()){
                            solved = false;
                            break;
                        }
                    }

                    if(solved){
                        // Memory management
                        unsat_clauses->clear();
                        delete unsat_clauses;

                        break;
                    }
                }

                // Initialise an empty vector I_t for the MIS of unsatisfied clauses associated with U_t of each thread
                auto indep_unsat_clauses = new vector<ClausesArray*>;
                for(int t = 0; t < n_threads; t++){
                    indep_unsat_clauses->push_back(new ClausesArray);
                }

                // Find a subset MIS I_t in each set of unsatisfied clauses U_t, in parallel (statically, one for each thread);
                // As described earlier, we visit the clauses in each U_t in a pseudo-random order (effectively shuffling)
                #pragma omp parallel for schedule(static, 1) default(none) shared(unsat_clauses, indep_unsat_clauses, caches)
                for(int t = 0; t < n_threads; t++){

                    for(auto clause = unsat_clauses->at(t)->begin();
                        clause != unsat_clauses->at(t)->end(); ++clause){ // Visit each clause in U_t; note that the iterator i is not used
                        // as an index; rather the LCG corresponding to the thread is. Since the LCG has period chunk_size,
                        // each time we update it with every iteration i, we ensure that no clause is visited twice. Hence,
                        // since i has a range from 0 to chunk_size - 1, we visit all the clauses in U_t in a non-linear
                        // fashion.

                        bool independent = true; // If current clause in U_t is dependent on some clause in the current
                        // state of the independent set I_t, then change flag to false.
                        for(auto& c: *(indep_unsat_clauses->at(t))){ // Check if dependent on any clause in I_t...
                            if(dependent_clauses((*clause), c, cache_depth,caches.at(t))){

                                // If dependent, update independent flag to false and break (no need to continue checking -
                                // the current clause in U_t will not be added to the independent set)
                                independent = false;

                                break;
                            }
                        }

                        if(independent){ // If independent, then add to independent set I_t
                            (*indep_unsat_clauses->at(t)).push_back(*clause);
                        }
                    } // At the end, I_t will be a MIS since all clauses in U_t have been exhausted
                }

                // Join the n_thread MISs {I_t} into a single MIS (through the greedy_parallel_mis_join() function)
                auto max_indep_unsat_clauses = greedy_parallel_mis_join(indep_unsat_clauses);

                for (Clause<uint32_t>* clause : *max_indep_unsat_clauses) {
                    unsat_clauses->at(clause->t_id)->erase(clause);
                }

                for (int t = 0; t < n_threads; t++) {
                    unsat_clauses->at(t)->resize(0);
                }

                statistics->avg_mis_size += max_indep_unsat_clauses->size(); // Update statistics

                // In parallel and dynamically, re-sample the variables appearing in the clauses of the MIS constructed
                #pragma omp parallel for schedule(dynamic) default(none) shared(max_indep_unsat_clauses, rbg_ensemble, statistics)
                for(ull c = 0; c < max_indep_unsat_clauses->size(); c++){
                    int t_id = omp_get_thread_num(); // Get thread identifier

                    // Re--sample every variable in the clause using the RBG associated with the thread t_id
                    for(auto& l : *(max_indep_unsat_clauses->at(c))->literals){
                        (var_arr->vars)[l >> 1] = (rbg_ensemble[t_id])->sample();
                    }

                    // Update statistics
                    statistics->n_thread_resamples.at(t_id) += (max_indep_unsat_clauses->at(c))->literals->size();
                }

                delete max_indep_unsat_clauses; // Memory management
            }

            for(int t = 0; t < n_threads; t++){
                statistics->n_resamples += (statistics->n_thread_resamples).at(t); // Update statistics
            }
            statistics->avg_mis_size = (ull) (statistics->avg_mis_size / statistics->n_iterations); // Update statistics

            return statistics;
        }

        /* Serial Algorithmic Lovasz Local Lemma of Moser and Tardos (2010), with 'shuffling' of the clauses vector to
         * satisfy the 'there exists an unsatisfied clause' condition stochastically.
         */
//        Statistics* sequential_solve() const{
//            auto statistics = new Statistics;
//
//            Clause<tV>* unsat_clause;
//
//            // Get the system default random generator with a random seed
//            auto engine = new default_random_engine(std::random_device{}());
//
//            // Initialise an instance of a random boolean generator...
//            auto rbg = new RBG<default_random_engine>(*engine);
//
//            ull iterator = P_9223372036854775783 % n_clauses; // Starting seed for an LCG based index to the this.clauses vector,
//            // which pseudo-randomly visits all the clauses exactly once with
//            // a period of n_clauses, allowing for efficient 'shuffling' of the
//            // clauses vector;
//
//            /* Note that by shuffling through the clauses array, we are effectively breaking any 'local minima' in the search
//             * space; in between iterations, we do not get stuck successively in the same clause until it becomes satisfied, but
//             * rather we visit other unsatisfied clauses and satisfy those (which in turn may influence out original unsatisfied
//             * clause).
//             */
//
//            bool solved = false;
//            while(!solved){ // While there exists a clause c which is not satisfied...
//                bool unsat_exists = false;
//                statistics->n_iterations += 1; // Update statistics
//
//                for(ull i = 0; i < n_clauses; i++){ // For each clause...
//                    // Fetch the clause specified at the (class variable) index C, and check if it is satisfied...
//
//                    // If it is not satisfied, return a ptr to the Clause instance
//                    if((clauses->at(iterator))->is_not_satisfied(var_arr->vars)){
//                        unsat_clause = clauses->at(iterator);
//                        unsat_exists = true;
//                        break;
//                    }
//                    else{ /* Otherwise, calculate the next index C by means of the LCG; note that the LCG has a period of
//                   * n_clauses and hence we iterate over each clause exactly once, in a non-linear fashion; this is
//                   * equivalent to a pseudo-random shuffling of the clauses vector, and is required for optimal
//                   * convergence of the Algorithmic Lovasz Local Lemma method.
//                   */
//                        iterator = (iterator + P_9223372036854775783) % n_clauses;
//                    }
//                }
//
//                if(unsat_exists){
//                    // For every variable in the clause (obtained by left shifting by 1 the literal encoding), randomly
//                    // re-sample
//                    for(auto& l : *(unsat_clause->literals)){
//                        (var_arr->vars)[l >> 1] = rbg->sample();
//                    }
//
//                    statistics->n_resamples += unsat_clause->literals->size(); // Update statistics
//                }
//                else{
//                    solved = true;
//                }
//            }
//
//            return statistics;
//        }

        // =============================================================================================================
        // ------------------------------------------------- UTILITIES -------------------------------------------------
        // =============================================================================================================

        // Given two Clause instances, establishes whether the two are dependent or not. The criteria of dependency is
        // having one or more variables in common between the literals of the two clauses.
        bool dependent_clauses(Clause<tV>* c1, Clause<tV>* c2, bool use_cache, ClauseCache* cached_clause_dependencies){
            ClausePair pair_cIDs;

            if(use_cache){
                pair_cIDs = minmax(c1, c2);
                if(auto d = cached_clause_dependencies->find(pair_cIDs); d != cached_clause_dependencies->end()){
                    return false;
                }
            }

            // For every pair (x, y) of literals between the two clauses...
            bool dependent = false;
            for(auto& l1 : *c1->literals){
                for(auto& l2 : *c2->literals){
                    // If the variable associated with x is equal to the variable associated with y, then the clauses are
                    // dependent and hence we return true (no need to check further - we require AT LEAST one).
                    if((l1 >> 1) == (l2 >> 1)){
                        dependent = true;
                        break;
                    }
                }

                if(dependent){
                    break;
                }
            }

            if(use_cache && !dependent){
                if(cache_depth < cached_clause_dependencies->size()){
                    cached_clause_dependencies->erase(cached_clause_dependencies->begin());
                }

                cached_clause_dependencies->insert(pair_cIDs);
            }

            // Otherwise, we return false (independent if FOR ALL literal pairs, the underlying variables are distinct)
            return dependent;
        }

        // Utility function for constructing a maximally independent set (MIS) from two other such sets
        ClausesArray* greedy_mis_join(ClausesArray* set1, ClausesArray* set2){
            auto result = new ClausesArray; // The resulting MIS from joining set1 and set2

            /* We begin by iterating across all the elements in set1 and checking it against every element in set2. If
             * for some element x in set1 there is an element y in set2 such that x and y are independent, then we delete
             * y from set2 (on the fly). Hence, the resulting MIS will always contain set1 as a subset (where set1 is
             * assumed to be a MIS).
             */
            auto idx1 = set1->begin();
            while(idx1 != set1->end()){ // Consider an element x in set1...
                auto idx2 = set2->begin();
                while(idx2 != set2->end()){ // for every element y in set2
                    if(dependent_clauses(*idx1, *idx2, false, nullptr)){ // if y and x are dependent, remove y from set2
                        idx2 = set2->erase(idx2);
                    }
                    else{ // else maintain y in set2 and check the next element y' (if any) in set2
                        ++idx2;
                    }
                }

                result->push_back(*idx1); // In any case, add x from set1 to the resulting MIS
                idx1 = set1->erase(idx1); // Erase x from set1 (simply for memory management)
            }

            /* As a result, all the remaining elements in set2 are independent from all the elements in set1, i.e. from
             * all the elements in the resulting MIS thus far. Hence we simply take the union of the two sets (the current
             * resulting set and the augmented set2).
             */
            auto idx2 = set2->begin();
            while(idx2 != set2->end()){ // Add every element in the augmented set2 to the resulting MIS
                result->push_back(*idx2);
                idx2 = set2->erase(idx2); // memory management
            }

            return result;
        }

        // Utility function for recursively and in parallel joining k disjoint maximally independent sets into a single one
        ClausesArray* greedy_parallel_mis_join(vector<ClausesArray*>* sets){
            // BASE CASE
            if(sets->size() == 1){ // If only a single MIS is passed as input, then simply return it
                return sets->at(0);
            }
            else{ // RECURSIVE CASE
                /* The join works by pairing the independent sets passed as input, and (in parallel) each pair is joined
                 * into a single MIS using the greedy_mis_join() algorithm. Note that in the case of an odd number of
                 * sets, one is not paired.
                 *
                 * The greedy_parallel_mis_join is then called once again on the resulting sets, until the base case is
                 * satisfied.
                 */

                auto joined_sets = new vector<ClausesArray*>; // Resulting MIS from pair--wise joins
                int offset = 0; // Index offset from which to start populating joined_sets (if odd then index 0 is
                // populated with the largest input MIS, and the remaining are paired off and joined_sets is populated
                // with the joined sets from index 1 onwards ie. offset is set to 1.

                if(sets->size() % 2){ // If odd number of initial sets...
                    offset = 1; // Set offset to 1

                    // Find the input MIS with the largest size
                    auto max_idx = sets->begin();
                    for(auto idx = sets->begin(); idx != sets->end(); ++idx){
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

                // In parallel join each pair of input MISs using the greedy_mis_join() algorithm (one for each thread;
                // number of pairs does not exceed n_threads)
                #pragma omp parallel for schedule(static, 1) default(none) shared(n_pairs, offset, sets, joined_sets)
                for(int t = 0; t < n_pairs; t++){
                    // Join pairs greedily into a single independent set...
                    if(sets->at((2 * t) + 1)->size() < sets->at(2 * t)->size()) {
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
};


#endif //ALLLSATISFIABILITYSOLVER_SATINSTANCE_H