//
// Created by Xandru Mifsud on 04/10/2021.
//

#ifndef ALLLSATISFIABILITYSOLVER_SUBSATINSTANCE_H
#define ALLLSATISFIABILITYSOLVER_SUBSATINSTANCE_H

#include <iostream>
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

        template<typename T>
        void partition(T t, blocks (*p)(T)){
            clausePartition = p(t);
        }

    private:
        // Largest prime number (2^64 - 59) that fits in a 64-bit register; note that this limits the number of clauses
        // that can be in the SAT instance to 2^64 - 59 (which should be reasonably large enough...we hope...)
        const ull P_2e64_m59 = 18446744073709551557;

        vector<RBG<default_random_engine>*> rbg_ensemble;

        ull parallel_resample;
        ull n_clauses;
        ull C;
};


#endif //ALLLSATISFIABILITYSOLVER_SUBSATINSTANCE_H
