//
// Created by Xandru Mifsud on 04/10/2021.
//

#ifndef ALLLSATISFIABILITYSOLVER_SUBSATINSTANCE_H
#define ALLLSATISFIABILITYSOLVER_SUBSATINSTANCE_H

#include <iostream>
#include <utility>
#include <vector>
#include <set>

#include <Eigen/Dense>
#include <omp.h>

#include "../core/VariablesArray.h"
#include "../core/Clause.h"

using namespace Eigen;
using namespace std;

typedef vector<set<ull>*>* blocks;

class SubSATInstance{
    public:
        blocks clausePartition = nullptr;
        vector<Clause*>* clauses;
        VariablesArray* var_arr;

        SubSATInstance(VariablesArray* variables, vector<Clause*>* clauses, ull parallel_resample = 0);

        void solve();
        Clause* is_satisfied();
        bool is_ALLL_compatible() const;

        template<typename T>
        void partition(T t, blocks (*p)(T)){
            clausePartition = p(t);
        }

    private:
        // Prime number for LCG over the clauses array (which should be reasonably large enough...we hope...)
        const ull P_9223372036854775783 = 9223372036854775783;

        vector<RBG<default_random_engine>*> rbg_ensemble;

        ull parallel_resample;
        ull n_clauses;
        ull C;
};


#endif //ALLLSATISFIABILITYSOLVER_SUBSATINSTANCE_H
