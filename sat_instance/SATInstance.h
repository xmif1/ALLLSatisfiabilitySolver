//
// Created by Xandru Mifsud on 28/09/2021.
//

#ifndef ALLLSATISFIABILITYSOLVER_SATINSTANCE_H
#define ALLLSATISFIABILITYSOLVER_SATINSTANCE_H

#include <iostream>
#include <utility>
#include <vector>
#include <random>
#include <set>
#include <map>

#include <omp.h>

#include "../cnf_io/cnf_io.h"

#include "../core/RandomBoolGenerator.h"
#include "../core/VariablesArray.h"
#include "../core/Clause.h"

using namespace std;

typedef vector<Clause*> ClausesArray;

typedef struct Statistics{
    ull n_iterations = 0;
    ull n_resamples = 0;
    ull avg_mis_size = 0;
    vector<ull> n_thread_resamples;
} Statistics;

/* Represents a CNF SAT instance loaded from a DIMACS-formatted file, providing a number of functions to:
 *   i. Generate the Laplacian of the (full) dependency graph between the clauses of the instance.
 *  ii. Solve, by means of the Algorithmic Lovasz Local Lemma (Moser & Tardos, 2010), the SAT instance.
 * iii. Check whether a variable assignment satisfies the SAT instance.
 */
class SATInstance{
    public:
        uint32_t n_vars;
        uint32_t n_literals;
        uint32_t n_clauses;

        VariablesArray* var_arr;
        ClausesArray* clauses = new vector<Clause*>;

        explicit SATInstance(const string& cnf_file_name, int n_threads = 0);

        Statistics* solve();
        bool verify_validity() const;

    private:
        // Prime number for LCG over the clauses array (which should be reasonably large enough...we hope...)
        const ull P_9223372036854775783 = 9223372036854775783;
        int n_threads;

        Statistics* parallel_solve();
        Statistics* sequential_solve() const;
        ClausesArray* parallel_k_partite_mis(vector<ClausesArray*>* sets);
        static ClausesArray* bipartite_mis(ClausesArray* set1, ClausesArray* set2);
};


#endif //ALLLSATISFIABILITYSOLVER_SATINSTANCE_H