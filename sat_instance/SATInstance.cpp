//
// Created by Xandru Mifsud on 28/09/2021.
//

#include "SATInstance.h"

SATInstance::SATInstance(const string& cnf_file_name){
    int l_num, c_num, v_num;
    int* l_c_num;
    int* l_val;

    if(cnf_header_read(cnf_file_name, &v_num, &c_num, &l_num)){
        cout << "The header information could not be read. Exiting..." << endl;
        exit(1);
    }

    cout << "V_NUM = " << v_num << endl;
    cout << "C_NUM = " << c_num << endl;
    cout << "L_NUM = " << l_num << endl;

    l_c_num = new int[c_num];
    l_val = new int[l_num];

    cnf_data_read(cnf_file_name, v_num, c_num, l_num, l_c_num, l_val);

    int c, l, l_c;

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

    n_vars = v_num;
    n_clauses = c_num;
    C = P_2e64_m59 % n_clauses;
}

bool SATInstance::dependent_clauses(Clause* c1, Clause* c2){
    for(int i = 0; i < c1->n_literals; i++){
        for(int j = 0; j < c2->n_literals; j++){
            if((c1->literals)[i] == (c2->literals)[i] || (c1->literals)[i] == -((c2->literals)[i])){
                return true;
            }
        }
    }

    return false;
}

MatrixXd* SATInstance::getDependancyGraphLaplacian(){
    auto laplacian = new MatrixXd(n_clauses, n_clauses);
    laplacian->setZero();

    for(ull i = 0; i < n_clauses; i++){
        for(ull j = i+1; j < n_clauses; j++){
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

Clause* SATInstance::is_satisfied(VariablesArray* var_arr){
    for(ull i = 0; i < n_clauses; i++){
        if((clauses.at(C))->is_not_satisfied(var_arr)){
            return clauses.at(C);
        }
        else{
            C = (C + P_2e64_m59) % n_clauses;
        }
    }

    return nullptr;
}

VariablesArray* SATInstance::solve(){
    default_random_engine engine(std::random_device{}());
    RBG<default_random_engine> rbg(engine);

    auto var_arr = new VariablesArray(n_vars);

    Clause* c = is_satisfied(var_arr);
    while(c){
        for(lli i = 0; i < c->n_literals; i++){
            (*var_arr)[(c->literals)[i]] = rbg.sample();
        }

        c = is_satisfied(var_arr);
    }

    sat = var_arr;

    return var_arr;
}