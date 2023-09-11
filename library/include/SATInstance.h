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

#include "RandomBoolGenerator.h"
#include "ClauseGenerator.h"
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
template <typename T>
class SATInstance{

    public:
        using ClauseArray = typename Clause<T>::ClauseArray;

        T n_vars;
        ull n_clauses = 0;

        VariablesArray<T>* var_arr;

        // Constructor for a SATInstance object
        SATInstance(VariablesArray<T>* var_arr, int n_threads) {
            this->var_arr = var_arr;                        // Encoded variables array
            this->n_threads = n_threads;                    // Number of threads for solving (if using parallel solver)

            n_vars = var_arr->n_vars;                       // Number of variables in SAT instance
        }

        // SAT solver based on the Algorithmic Lovasz Local Lemma of Moser and Tardos (2010); wrapper to either sequential
        // or parallel solver depending on whether 1 or more threads specified during instantiation.
        Statistics* solve (vector<ClauseArray*>* clauses) {
            for (auto c : *clauses) {    // Number of clauses in SAT instance
                n_clauses += c->size();
            }

            return parallel_solve(clauses, false);
        }

        // Dynamic-generation SAT solver based on the Algorithmic Lovasz Local Lemma of Moser and Tardos (2010);
        // wrapper to either sequential or parallel solver depending on whether 1 or more threads specified during instantiation.
        Statistics* solve (Clause<T>* (*getEnumeratedClause)(T, unsigned short int), ull n_clauses, T batch_size) {
            auto statistics = new Statistics;
            
            this->n_clauses = n_clauses;
            T t_n_clauses = (T) n_clauses / n_threads;

            auto clauses = new vector<ClauseArray*>;
            auto generators = new vector<ClauseGenerator<T>*>;
            for (int t = 0; t < this->n_threads; t++) {
                T offset = (T) t * t_n_clauses;
                if (t == n_threads - 1) {
                    t_n_clauses = n_clauses - offset;
                }

                generators->push_back(new ClauseGenerator<T>(getEnumeratedClause, t, t_n_clauses, offset, batch_size));
            }
            
            bool solved = false;
            while (!solved) {
                solved = true;
                statistics->n_iterations += 1;

                bool finishedYielding = false;
                
                while (!finishedYielding) {
                    #pragma omp parallel for schedule(static, 1) default(none) shared(clauses, generators, n_threads)
                    for (int t = 0; t < n_threads; t++) {
                        clauses->push_back(generators->at(t)->yieldRandomClauseBatch());
                    }

                    finishedYielding = true;
                    for (int t = 0; t < n_threads; t++) {
                        if (!generators->at(t)->has_finished_yielding()) {
                            finishedYielding = false;
                            break;
                        }
                    }

                    auto result = parallel_solve(clauses, true);

                    statistics->avg_mis_size += result->avg_mis_size;
                    statistics->n_resamples += result->n_resamples;
                    for (int t = 0; t < n_threads; t++) {
                        statistics->n_thread_resamples.at(t) += result->n_thread_resamples.at(t);
                    }

                    #pragma omp parallel for schedule(static, 1) default(none) shared(clauses, n_threads)
                    for (int t = 0; t < n_threads; t++) {
                        // Memory management
                        for (auto c: *clauses->at(t)) {
                            c->literals->clear();
                            delete c->literals;
                            delete c;
                        }

                        clauses->at(t)->clear();
                    }

                    clauses->clear();
                }
                
                finishedYielding = false;
                while (!finishedYielding) {
                    #pragma omp parallel for schedule(static, 1) default(none) shared(clauses, generators, n_threads)
                    for (int t = 0; t < n_threads; t++) {
                        clauses->push_back(generators->at(t)->yieldOrderedClauseBatch());
                    }

                    finishedYielding = true;
                    for (int t = 0; t < n_threads; t++) {
                        if (!generators->at(t)->has_finished_yielding()) {
                            finishedYielding = false;
                            break;
                        }
                    }

                    solved = solved && verify_validity(clauses);

                    #pragma omp parallel for schedule(static, 1) default(none) shared(clauses, n_threads)
                    for (int t = 0; t < n_threads; t++) {
                        // Memory management
                        for (auto c: *clauses->at(t)) {
                            c->literals->clear();
                            delete c->literals;
                            delete c;
                        }

                        clauses->at(t)->clear();
                    }
                }
            }

            statistics->avg_mis_size /= statistics->n_iterations;

            return statistics;
        }

        // Convenience function for checking whether the assignments in var_arr represent a valid solution or not.
        bool verify_validity (vector<ClauseArray*>* clauses) const {
            volatile bool valid = true;

            #pragma omp parallel for schedule(static, 1) default(none) shared(clauses, valid)
            for (int t = 0; t < clauses->size(); t++) {
                for (auto clause = clauses->at(t)->begin(); clause != clauses->at(t)->end(); clause++) {
                    if(!valid) {
                        continue;
                    }

                    if((*clause)->is_not_satisfied(var_arr->vars)) {
                        valid = false;
                    }
                }
            }

            return valid;
        }

    private:
        // Prime number for LCG over the clauses array (which should be reasonably large enough...we hope...)
        int n_threads{};

        // =============================================================================================================
        // ----------------------------------------------- ALLL SOLVERS ------------------------------------------------
        // =============================================================================================================

        /* Parallel Algorithmic Lovasz Local Lemma of Moser and Tardos (2010), with parallel unsatisfied clause checking
         * and parallel divide-and-conquer generation of non-trivial independent sets of unsatisfied clauses; the sets
         * of independent unsatisfied clauses chosen by each thread is done in a pseudo-random manner.
         */
        Statistics* parallel_solve (vector<ClauseArray*>* clauses, bool stream) {
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

            vector<RBG<default_random_engine>*> rbg_ensemble;

            // Initialise ensembles of RBGs and clause_iterators...
            for (int i = 0; i < n_threads; i++) {
                // Get the system default random generator with a random seed
                auto engine = new default_random_engine(std::random_device{}());

                // Initialise an instance of a random boolean generator
                rbg_ensemble.push_back(new RBG<default_random_engine>(*engine));

                // Initialise statistics
                statistics->n_thread_resamples.push_back(0);
            }

            auto unsat_clauses = new vector<ClauseArray*>;

            // Initialise an empty vector U_t for the unsatisfied clauses found by each thread
            for (int t = 0; t < n_threads; t++) {
                unsat_clauses->push_back(new ClauseArray());
            }

            omp_set_num_threads(n_threads);
            while (true) { // While there exists a clause c which is not satisfied...
                statistics->n_iterations += 1; // Update statistics

                // Split clauses in n_threads batches {C_j}, in parallel through dynamic allocation
                #pragma omp parallel for schedule(static, 1) default(none) shared(clauses, unsat_clauses)
                for (int t = 0; t < n_threads; t++) {
                    for (auto clause = clauses->at(t)->begin(); clause != clauses->at(t)->end(); ++clause) {
                        if((*clause)->is_not_satisfied(var_arr->vars)) {
                            unsat_clauses->at(t)->push_back(*clause);
                        }
                    }
                }

                bool solved = true;
                for (auto unsatClauses : *unsat_clauses) {
                    if(!unsatClauses->empty()) {
                        solved = false;
                        break;
                    }
                }

                if (solved) {
                    // Memory management
                    unsat_clauses->clear();
                    delete unsat_clauses;

                    break;
                }

                // Initialise an empty vector I_t for the MIS of unsatisfied clauses associated with U_t of each thread
                auto indep_unsat_clauses = new vector<ClauseArray*>;
                for (int t = 0; t < n_threads; t++) {
                    indep_unsat_clauses->push_back(new ClauseArray);
                }

                // Find a subset MIS I_t in each set of unsatisfied clauses U_t, in parallel (statically, one for each thread);
                // As described earlier, we visit the clauses in each U_t in a pseudo-random order (effectively shuffling)
                #pragma omp parallel for schedule(static, 1) default(none) shared(unsat_clauses, indep_unsat_clauses)
                for (int t = 0; t < n_threads; t++) {

                    for (auto clause = unsat_clauses->at(t)->begin();
                        clause != unsat_clauses->at(t)->end(); ++clause) { 
                        // Visit each clause in U_t; note that the iterator i is not used as an index; 
                        // rather the LCG corresponding to the thread is. Since the LCG has period chunk_size,
                        // each time we update it with every iteration i, we ensure that no clause is visited twice. Hence,
                        // since i has a range from 0 to chunk_size - 1, we visit all the clauses in U_t in a non-linear
                        // fashion.

                        bool independent = true; // If current clause in U_t is dependent on some clause in the current
                        // state of the independent set I_t, then change flag to false.
                        for (auto& c: *(indep_unsat_clauses->at(t))) { // Check if dependent on any clause in I_t...
                            if (dependent_clauses((*clause), c)) {

                                // If dependent, update independent flag to false and break (no need to continue checking -
                                // the current clause in U_t will not be added to the independent set)
                                independent = false;

                                break;
                            }
                        }

                        if(independent) { // If independent, then add to independent set I_t
                            (*indep_unsat_clauses->at(t)).push_back(*clause);
                        }
                    } // At the end, I_t will be a MIS since all clauses in U_t have been exhausted
                }

                // Join the n_thread MISs {I_t} into a single MIS (through the greedy_parallel_mis_join() function)
                auto max_indep_unsat_clauses = greedy_parallel_mis_join(indep_unsat_clauses);

                for (int t = 0; t < n_threads; t++) {
                    unsat_clauses->at(t)->clear();
                }

                statistics->avg_mis_size += max_indep_unsat_clauses->size(); // Update statistics

                // In parallel and dynamically, re-sample the variables appearing in the clauses of the MIS constructed
                #pragma omp parallel for schedule(dynamic) default(none) shared(max_indep_unsat_clauses, rbg_ensemble, statistics)
                for (ull c = 0; c < max_indep_unsat_clauses->size(); c++) {
                    int t_id = omp_get_thread_num(); // Get thread identifier

                    // Re--sample every variable in the clause using the RBG associated with the thread t_id
                    for (auto& l : *(max_indep_unsat_clauses->at(c))->literals) {
                        (var_arr->vars)[l >> 1] = (rbg_ensemble[t_id])->sample();
                    }

                    // Update statistics
                    statistics->n_thread_resamples.at(t_id) += (max_indep_unsat_clauses->at(c))->literals->size();
                }

                max_indep_unsat_clauses->clear();
                delete max_indep_unsat_clauses; // Memory management

                if (stream) {
                    unsat_clauses->clear();
                    delete unsat_clauses;

                    break;
                }
            }

            for (int t = 0; t < n_threads; t++) {
                statistics->n_resamples += (statistics->n_thread_resamples).at(t); // Update statistics
            }
            statistics->avg_mis_size = (ull) (statistics->avg_mis_size / statistics->n_iterations); // Update statistics

            return statistics;
        }

        // =============================================================================================================
        // ------------------------------------------------- UTILITIES -------------------------------------------------
        // =============================================================================================================

        // Given two Clause instances, establishes whether the two are dependent or not. The criteria of dependency is
        // having one or more variables in common between the literals of the two clauses.
        bool dependent_clauses (Clause<T>* c1, Clause<T>* c2) {
            // For every pair (x, y) of literals between the two clauses...
            bool dependent = false;
            for (auto& l1 : *c1->literals) {
                for (auto& l2 : *c2->literals) {
                    // If the variable associated with x is equal to the variable associated with y, then the clauses are
                    // dependent and hence we return true (no need to check further - we require AT LEAST one).
                    if ((l1 >> 1) == (l2 >> 1)) {
                        dependent = true;
                        break;
                    }
                }

                if(dependent) {
                    break;
                }
            }

            // Otherwise, we return false (independent if FOR ALL literal pairs, the underlying variables are distinct)
            return dependent;
        }

        // Utility function for constructing a maximally independent set (MIS) from two other such sets
        ClauseArray* greedy_mis_join(ClauseArray* set1, ClauseArray* set2) {
            auto result = new ClauseArray; // The resulting MIS from joining set1 and set2

            /* We begin by iterating across all the elements in set1 and checking it against every element in set2. If
             * for some element x in set1 there is an element y in set2 such that x and y are independent, then we delete
             * y from set2 (on the fly). Hence, the resulting MIS will always contain set1 as a subset (where set1 is
             * assumed to be a MIS).
             */
            auto idx1 = set1->begin();
            while (idx1 != set1->end()) { // Consider an element x in set1...
                auto idx2 = set2->begin();
                while (idx2 != set2->end()) { // for every element y in set2
                    if (dependent_clauses(*idx1, *idx2)) { // if y and x are dependent, remove y from set2
                        idx2 = set2->erase(idx2);
                    } else { // else maintain y in set2 and check the next element y' (if any) in set2
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
            while (idx2 != set2->end()) { // Add every element in the augmented set2 to the resulting MIS
                result->push_back(*idx2);
                idx2 = set2->erase(idx2); // memory management
            }

            return result;
        }

        // Utility function for recursively and in parallel joining k disjoint maximally independent sets into a single one
        ClauseArray* greedy_parallel_mis_join(vector<ClauseArray*>* sets) {
            // BASE CASE
            if (sets->size() == 1) { // If only a single MIS is passed as input, then simply return it
                auto ret_set = sets->at(0);
                delete sets;

                return ret_set;
            } else { // RECURSIVE CASE
                /* The join works by pairing the independent sets passed as input, and (in parallel) each pair is joined
                 * into a single MIS using the greedy_mis_join() algorithm. Note that in the case of an odd number of
                 * sets, one is not paired.
                 *
                 * The greedy_parallel_mis_join is then called once again on the resulting sets, until the base case is
                 * satisfied.
                 */

                auto joined_sets = new vector<ClauseArray*>; // Resulting MIS from pair--wise joins
                int offset = 0; // Index n_yielded_clauses from which to start populating joined_sets (if odd then index 0 is
                // populated with the largest input MIS, and the remaining are paired off and joined_sets is populated
                // with the joined sets from index 1 onwards ie. n_yielded_clauses is set to 1.

                if (sets->size() % 2) { // If odd number of initial sets...
                    offset = 1; // Set n_yielded_clauses to 1

                    // Find the input MIS with the largest size
                    auto max_idx = sets->begin();
                    for (auto idx = sets->begin(); idx != sets->end(); ++idx) {
                        if ((*max_idx)->size() < (*idx)->size()) {
                            max_idx = idx;
                        }
                    }

                    joined_sets->push_back(*max_idx); // Append largest MIS to joined_sets
                    sets->erase(max_idx); // Memory management
                } // The remaining input MISs are paired off and joined into a single MIS

                auto n_pairs = (int) (sets->size() / 2); // Number of sets resulting from pairing
                for (int t = 0; t < n_pairs; t++) { // Initialisation
                    joined_sets->push_back(nullptr);
                }

                // In parallel join each pair of input MISs using the greedy_mis_join() algorithm (one for each thread;
                // number of pairs does not exceed n_threads)
                #pragma omp parallel for schedule(static, 1) default(none) shared(n_pairs, offset, sets, joined_sets)
                for (int t = 0; t < n_pairs; t++) {
                    // Join pairs greedily into a single independent set...
                    if (sets->at((2 * t) + 1)->size() < sets->at(2 * t)->size()) {
                        joined_sets->at(offset + t) = greedy_mis_join(sets->at(2 * t), sets->at((2 * t) + 1));
                    } else {
                        joined_sets->at(offset + t) = greedy_mis_join(sets->at((2 * t) + 1), sets->at(2 * t));
                    }

                    // Memory management (we only require the resulting joined set from two initial sets)
                    sets->at(2*t)->clear();
                    delete sets->at(2*t);

                    sets->at((2*t) + 1)->clear();
                    delete sets->at((2*t) + 1);
                }

                sets->clear();
                delete sets; // Memory management

                return greedy_parallel_mis_join(joined_sets); // Recursive join of the resulting MISs
            }
        }
};


#endif //ALLLSATISFIABILITYSOLVER_SATINSTANCE_H