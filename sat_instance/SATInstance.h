//
// Created by Xandru Mifsud on 28/09/2021.
//

#ifndef ALLLSATISFIABILITYSOLVER_SATINSTANCE_H
#define ALLLSATISFIABILITYSOLVER_SATINSTANCE_H

#include <random>

#include "../cnf_io/cnf_io.h"

#include "SubSATInstance.h"
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
        ull n_vars;
        ull n_literals;
        ull n_clauses;
        VariablesArray* var_arr;

        SATInstance(const string& cnf_file_name, vector<Clause*>* clauses);

        VariablesArray* solve(vector<SubSATInstance*>* subInstances, bool parallel = true) const;
        vector<SubSATInstance*>* createSubSATInstances(vector<vector<Clause*>*>* components, ull parallel_resample = 0) const;
        static pair<vector<MatrixXd*>*, vector<vector<Clause*>*>*> getDependencyGraph(vector<Clause*>* clauses);

        template<typename T>
        void partition(vector<SubSATInstance*>* subInstances, vector<T>* ts, blocks (*p)(T)){
            for(ull i = 0; i < subInstances->size(); i++){
                (subInstances->at(i))->partition(ts->at(i), p);
            }
        }

    private:
        static bool dependent_clauses(Clause* c1, Clause* c2);
};


#endif //ALLLSATISFIABILITYSOLVER_SATINSTANCE_H
