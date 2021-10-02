//
// Created by Xandru Mifsud on 28/09/2021.
//

#ifndef ALLLSATISFIABILITYSOLVER_SATINSTANCE_H
#define ALLLSATISFIABILITYSOLVER_SATINSTANCE_H

#include <iostream>
#include <random>
#include <vector>
#include <set>

#include <Eigen/Dense>
#include <cnf_io.hpp>

#include "../core/VariablesArray.h"
#include "../core/Clause.h"
#include "../core/RandomBoolGenerator.h"

using namespace Eigen;
using namespace std;

/* Represents a CNF SAT instance loaded from a DIMACS-formatted file, providing a number of functions to:
 *   i. Generate the Laplacian of the (full) dependency graph between the clauses of the instance.
 *  ii. Solve, by means of the Algorithmic Lovasz Local Lemma (Moser & Tardos, 2010), the SAT instance.
 * iii. Check whether a variable assignment satisfies the SAT instance.
 */
class SATInstance{
    public:
        VariablesArray* sat = nullptr;
        vector<Clause*> clauses;
        ull n_vars;

        explicit SATInstance(const string& cnf_file_name);

        pair<MatrixXd*, vector<vector<ull>*>*> getDependencyGraph();
        VariablesArray* solve();
        Clause* is_satisfied(VariablesArray* var_arr);

    private:
        // Largest prime number (2^64 - 59) that fits in a 64-bit register; note that this limits the number of clauses
        // that can be in the SAT instance to 2^64 - 59 (which should be reasonably large enough...we hope...)
        const ull P_2e64_m59 = 18446744073709551557;

        ull n_clauses;
        ull C;

        static bool dependent_clauses(Clause* c1, Clause* c2);
};


#endif //ALLLSATISFIABILITYSOLVER_SATINSTANCE_H
