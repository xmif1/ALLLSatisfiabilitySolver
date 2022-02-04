//
// Created by Xandru Mifsud on 28/09/2021.
//

#ifndef ALLLSATISFIABILITYSOLVER_SATINSTANCE_H
#define ALLLSATISFIABILITYSOLVER_SATINSTANCE_H

#include <random>
#include <map>
#include <omp.h>

#include "../cnf_io/cnf_io.h"

#include "SubSATInstance.h"
#include "../core/RandomBoolGenerator.h"

using namespace std;

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

        SATInstance(const string& cnf_file_name, vector<Clause*>* clauses);

        VariablesArray* solve(vector<SubSATInstance*>* subInstances, bool parallel = true) const;
        vector<SubSATInstance*>* createSubSATInstances(vector<vector<Clause*>*>* components, int parallel_resample = 0) const;
        bool verify_validity(vector<Clause*>* clauses) const;

        static vector<vector<Clause*>*>* getDependencyGraphComponents(vector<Clause*>* clauses);
        static bool dependent_clauses(Clause* c1, Clause* c2);
};


#endif //ALLLSATISFIABILITYSOLVER_SATINSTANCE_H