//
// Created by Xandru Mifsud on 28/09/2021.
//

#ifndef ALLLSATISFIABILITYSOLVER_SATINSTANCE_H
#define ALLLSATISFIABILITYSOLVER_SATINSTANCE_H

#include <unordered_map>
#include <iostream>
#include <fstream>
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

            return parallel_solve(clauses, new ClauseArray(),false, true);
        }

        // Dynamic-generation SAT solver based on the Algorithmic Lovasz Local Lemma of Moser and Tardos (2010);
        // wrapper to either sequential or parallel solver depending on whether 1 or more threads specified during instantiation.
        Statistics* solve (Clause<T>* (*getEnumeratedClause)(T, unsigned short int), ull n_clauses, T batch_size) {
            auto statistics = new Statistics;
            
            this->n_clauses = n_clauses;
            T t_n_clauses = (T) n_clauses / n_threads;

            auto generators = new vector<ClauseGenerator<T>*>;
            for (int t = 0; t < this->n_threads; t++) {
                statistics->n_thread_resamples.push_back(0);

                T offset = (T) t * t_n_clauses;
                if (t == n_threads - 1) {
                    t_n_clauses = n_clauses - offset;
                }

                generators->push_back(new ClauseGenerator<T>(getEnumeratedClause, t, t_n_clauses, offset, batch_size));
            }

            volatile bool solved = false;
            while (!solved) {
                solved = true;
                statistics->n_iterations += 1;

                auto mis = new ClauseArray();

                bool finishedYielding = false;
                cout << "New solve iteration..." << endl;

                while (!finishedYielding) {
                    auto clauses = new vector<ClauseArray*>;

                    for (int t = 0; t < n_threads; t++) {
                        clauses->push_back(nullptr);
                    }

                    #pragma omp parallel for schedule(static, 1) default(none) shared(clauses, generators, var_arr, n_threads)
                    for (int t = 0; t < n_threads; t++) {
                        clauses->at(t) = generators->at(t)->yieldRandomUNSATClauseBatch(var_arr->vars);
                    }

                    finishedYielding = true;
                    for (int t = 0; t < n_threads; t++) {
                        if (!generators->at(t)->has_finished_yielding()) {
                            finishedYielding = false;
                            break;
                        }
                    }

                    auto result = parallel_solve(clauses, mis, true, finishedYielding);

                    statistics->avg_mis_size += result->avg_mis_size;
                    statistics->n_resamples += result->n_resamples;
                    for (int t = 0; t < n_threads; t++) {
                        statistics->n_thread_resamples.at(t) += result->n_thread_resamples.at(t);
                    }
                }

                delete mis;

                #pragma omp parallel for schedule(static, 1) default(none) shared(generators, solved, n_threads, var_arr)
                for (int t = 0; t < n_threads; t++) {
                    for (int k = 0; k < generators->at(t)->n_clauses; k++) {
                        if (!solved) {
                            continue;
                        }

                        if(generators->at(t)->yieldNextClause()->is_not_satisfied(var_arr->vars)) {
                            solved = false;
                        }
                    }
                }
            }

            statistics->avg_mis_size /= statistics->n_iterations;

            return statistics;
        }

        // Convenience function for checking whether the assignments in var_arr represent a valid solution or not.
        bool verify_validity (vector<ClauseArray*>* clauses) const {
            volatile bool valid = true;

            #pragma omp parallel for schedule(static, 1) default(none) shared(clauses, valid, var_arr)
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

        void writeDIMACS(Clause<T>* (*getEnumeratedClause)(T, unsigned short int), ofstream* out_f) {
            auto generator = new ClauseGenerator<T>(getEnumeratedClause, 0, n_clauses, 0, n_clauses);

            *out_f << "p cnf " << n_vars << " " << n_clauses << endl;

            for (ull i = 0; i < n_clauses; i++){
                auto clause = generator->yieldNextClause();

                for(auto& l : *(clause->literals)) {
                    if (l & 1) {
                        *out_f << " " << to_string(-((intmax_t) (l >> 1)) - 1);
                    } else {
                        *out_f << " " << to_string((l >> 1) + 1);
                    }
                }

                *out_f << " 0" << endl;

                clause->literals->clear();
                delete clause->literals;
                delete clause;
            }
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
        Statistics* parallel_solve (vector<ClauseArray*>* clauses, ClauseArray* mis, bool stream, bool resample) {
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

            // Initialise statistics...
            for (int i = 0; i < n_threads; i++) {
                statistics->n_thread_resamples.push_back(0);
            }

            omp_set_num_threads(n_threads);
            while (true) { // While there exists a clause c which is not satisfied...
                statistics->n_iterations += 1; // Update statistics
                vector<ClauseArray*>* unsat_clauses;

                if (!stream) {
                    unsat_clauses = new vector<ClauseArray *>;

                    // Initialise an empty vector U_t for the unsatisfied clauses found by each thread
                    for (int t = 0; t < n_threads; t++) {
                        unsat_clauses->push_back(new ClauseArray());
                    }

                    // Split clauses in n_threads batches {C_j}, in parallel through dynamic allocation
                    #pragma omp parallel for schedule(static, 1) default(none) shared(clauses, unsat_clauses, n_threads)
                    for (int t = 0; t < n_threads; t++) {
                        for (auto clause = clauses->at(t)->begin(); clause != clauses->at(t)->end(); ++clause) {
                            if((*clause)->is_not_satisfied(var_arr->vars)) {
                                unsat_clauses->at(t)->push_back(*clause);
                            }
                        }
                    }
                } else {
                    unsat_clauses = clauses; // Oracle in streaming scenario only supplies unsat clauses
                }

                if (!stream && check_if_noUNSAT(unsat_clauses)) {
                    break;
                }

                // Join the n_thread MISs {I_t} into a single MIS (through the greedy_parallel_mis_join() function)
                get_mis_parallel(unsat_clauses, mis, stream);
                statistics->avg_mis_size += mis->size(); // Update statistics

                if (resample) {
                    resample_clauses(mis, statistics);

                    if (stream) { // Memory management - can discard clauses at this point
                        #pragma omp parallel for schedule(dynamic) default(none) shared(mis)
                        for (int t = 0; t < mis->size(); t++) {
                            mis->at(t)->literals->clear();
                            delete mis->at(t)->literals;
                            delete mis->at(t);
                        }
                    }

                    mis->clear();
                }

                if (stream) {
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

        bool check_if_noUNSAT(vector<ClauseArray*>* unsat_clauses) {
            for (auto unsatClauses : *unsat_clauses) {
                if(!unsatClauses->empty()) {
                    return false;
                }
            }

            // Memory management
            unsat_clauses->clear();
            delete unsat_clauses;

            return true;
        }

        void resample_clauses(ClauseArray* mis, Statistics* statistics) {
            vector<RBG<default_random_engine>*> rbg_ensemble;

            // Initialise ensembles of RBGs and clause_iterators...
            for (int i = 0; i < n_threads; i++) {
                // Get the system default random generator with a random seed
                auto engine = new default_random_engine(std::random_device{}());

                // Initialise an instance of a random boolean generator
                rbg_ensemble.push_back(new RBG<default_random_engine>(*engine));
            }

            // In parallel and dynamically, re-sample the variables appearing in the clauses of the MIS constructed
            #pragma omp parallel for schedule(dynamic) default(none) shared(mis, rbg_ensemble, statistics)
            for (ull c = 0; c < mis->size(); c++) {
                int t_id = omp_get_thread_num(); // Get thread identifier

                // Re--sample every variable in the clause using the RBG associated with the thread t_id
                for (auto& l : *(mis->at(c))->literals) {
                    (var_arr->vars)[l >> 1] = (rbg_ensemble[t_id])->sample();
                }

                // Update statistics
                statistics->n_thread_resamples.at(t_id) += (mis->at(c))->literals->size();
            }
        }

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

        void get_mis_parallel(vector<ClauseArray*>* sets, ClauseArray* mis, bool memory_manage) {
            if (!mis->empty()) { // Filter out clauses which are dependent on current MIS
                for (auto clause : *mis) {
                    #pragma omp parallel for schedule(static, 1) default(none) shared(sets, clause, memory_manage)
                    for (int k = 0; k < sets->size(); k++) {
                        auto idx = sets->at(k)->begin();
                        while (idx != sets->at(k)->end()) { // for every element y in set
                            if (dependent_clauses(clause, *idx)) { // if y and x are dependent, remove y from set2
                                if (memory_manage) { // Can completely discard clause
                                    (*idx)->literals->clear();
                                    delete (*idx)->literals;
                                    delete (*idx);
                                }

                                idx = sets->at(k)->erase(idx);
                            } else { // else maintain y in set and check the next element y' (if any) in set2
                                ++idx;
                            }
                        }
                    }
                }
            }

            int t = 0;
            while (!sets->empty()) { // Extend input MIS from batches until all batches are exhausted
                t = (t + 1) % sets->size();
                if (sets->at(t)->empty()) {
                    sets->at(t)->clear();
                    delete sets->at(t);

                    sets->erase(sets->begin() + t);
                    continue;
                }

                auto clause = sets->at(t)->at(0);
                sets->at(t)->erase(sets->at(t)->begin());

                mis->push_back(clause);

                #pragma omp parallel for schedule(static, 1) default(none) shared(sets, clause, memory_manage)
                for (int k = 0; k < sets->size(); k++) {
                    auto idx = sets->at(k)->begin();
                    while (idx != sets->at(k)->end()) { // for every element y in set
                        if (dependent_clauses(clause, *idx)) { // if y and x are dependent, remove y from set2
                            if (memory_manage) { // Can completely discard clause
                                (*idx)->literals->clear();
                                delete (*idx)->literals;
                                delete (*idx);
                            }

                            idx = sets->at(k)->erase(idx);
                        } else { // else maintain y in set and check the next element y' (if any) in set2
                            ++idx;
                        }
                    }
                }
            }

            sets->clear();
            delete sets;

            return mis;
        }
};


#endif //ALLLSATISFIABILITYSOLVER_SATINSTANCE_H