//
// Created by Xandru Mifsud on 04/10/2021.
//

#ifndef ALLLSATISFIABILITYSOLVER_SUBSATINSTANCE_H
#define ALLLSATISFIABILITYSOLVER_SUBSATINSTANCE_H

#include <iostream>
#include <utility>
#include <vector>
#include <set>

#include <omp.h>

#include "../core/VariablesArray.h"
#include "../core/Clause.h"

using namespace std;

class SubSATInstance{
    public:
        vector<Clause*>* clauses;
        VariablesArray* var_arr;

        SubSATInstance(VariablesArray* variables, vector<Clause*>* clauses, int n_threads = 0);

        void solve();
        Clause* sequential_is_satisfied();
        bool is_ALLL_compatible() const;

    private:
        // Prime number for LCG over the clauses array (which should be reasonably large enough...we hope...)
        const ull P_9223372036854775783 = 9223372036854775783;

        vector<RBG<default_random_engine>*> rbg_ensemble;

        int n_threads;
        uint32_t n_clauses;
        ull C;

        vector<Clause*>* parallel_k_partite_mis(vector<vector<Clause*>*>* sets);
        static vector<Clause*>* bipartite_mis(vector<Clause*>* set1, vector<Clause*>* set2);
};


#endif //ALLLSATISFIABILITYSOLVER_SUBSATINSTANCE_H
