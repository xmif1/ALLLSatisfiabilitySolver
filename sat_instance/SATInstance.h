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

class SATInstance{
    public:
        VariablesArray* sat = nullptr;
        vector<Clause*> clauses;
        ull n_vars;

        explicit SATInstance(const string& cnf_file_name);

        MatrixXd* getDependancyGraphLaplacian();
        VariablesArray* solve();
        Clause* is_satisfied(VariablesArray* var_arr);

    private:
        const ull P_2e64_m59 = 18446744073709551557;
        ull n_clauses;
        ull C;

        bool dependent_clauses(Clause* c1, Clause* c2);
};


#endif //ALLLSATISFIABILITYSOLVER_SATINSTANCE_H
