#include <iostream>
#include <random>
#include <vector>
#include <set>

#include <Eigen/Dense>
#include <cnf_io.hpp>

#include "core/VariablesArray.h"
#include "core/Clause.h"
#include "core/RandomBoolGenerator.h"

using namespace Eigen;
using namespace std;


Clause* is_satisfied(const vector<Clause*>& clauses, VariablesArray* var_arr);
VariablesArray* solve(const vector<Clause*>& clauses, ull n_vars);
MatrixXd* getDependancyGraphLaplacian(vector<Clause*> clauses);
bool dependent_clauses(Clause* c1, Clause* c2);

int main(){
    string cnf_file_name = "/Users/xandrumifsud/PLLL/ALLLSatisfiabilitySolver/aim-50-1_6-yes1-4.cnf";
    int l_num; int* l_c_num; int* l_val;
    int c_num;
    int v_num;

    if(cnf_header_read(cnf_file_name, &v_num, &c_num, &l_num)){
        cout << endl;
        cout << "The header information could not be read." << endl;
        return 1;
    }

    cout << endl;
    cout << "V_NUM  = " << v_num << endl;
    cout << "C_NUM  = " << c_num << endl;
    cout << "L_NUM  = " << l_num << endl;

    l_c_num = new int[c_num];
    l_val = new int[l_num];

    cnf_data_read(cnf_file_name, v_num, c_num, l_num, l_c_num, l_val);

    int c;
    int l;
    int l_c;

    vector<Clause*> clauses;

    l = 0;
    for(c = 0; c < c_num; c++){
        // l_c_num[c] = # of signed literals in clause c
        lli* literals; literals = new lli[l_c_num[c]];

        for(l_c = 0; l_c < l_c_num[c]; l_c++){
            literals[l_c] = l_val[l];

            l += 1;
        }

        clauses.push_back(new Clause(literals, l_c_num[c]));
    }

    MatrixXd* laplacian = getDependancyGraphLaplacian(clauses);
    solve(clauses, v_num);

    return 0;
}

MatrixXd* getDependancyGraphLaplacian(vector<Clause*> clauses){
    ull N = clauses.size();
    auto laplacian = new MatrixXd(N, N);
    laplacian->setZero();

    for(ull i = 0; i < N; i++){
        for(ull j = i+1; j < N; j++){
            if(dependent_clauses(clauses.at(i), clauses.at(j))){
                (*laplacian)(i, j) = -1;
                (*laplacian)(j, i) = -1;
                (*laplacian)(i, i) += 1;
                (*laplacian)(j, j) += 1;
            }
        }
    }

    return laplacian;
}

bool dependent_clauses(Clause* c1, Clause* c2){
    for(int i = 0; i < c1->n_literals; i++){
        for(int j = 0; j < c2->n_literals; j++){
            if((c1->literals)[i] == (c2->literals)[i] || (c1->literals)[i] == -((c2->literals)[i])){
                return true;
            }
        }
    }

    return false;
}

VariablesArray* solve(const vector<Clause*>& clauses, ull n_vars){
    default_random_engine engine(std::random_device{}());
    RBG<default_random_engine> rbg(engine);

    auto var_arr = new VariablesArray(n_vars);

    Clause* c = is_satisfied(clauses, var_arr);
    while(c){
        for(lli i = 0; i < c->n_literals; i++){
            (*var_arr)[(c->literals)[i]] = rbg.sample();
        }

        c = is_satisfied(clauses, var_arr);
    }

    return var_arr;
}

Clause* is_satisfied(const vector<Clause*>& clauses, VariablesArray* var_arr){
    for(auto c: clauses){
        if(c->is_not_satisfied(var_arr)){
            return c;
        }
    }

    return nullptr;
}
