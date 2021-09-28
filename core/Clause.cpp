//
// Created by Xandru Mifsud on 27/09/2021.
//

#include "Clause.h"

Clause::Clause(lli* literals, ull n_literals){
    this->literals = literals;
    this->n_literals = n_literals;
}

bool Clause::is_not_satisfied(VariablesArray* var_arr) const{
    for(ull i = 0; i < n_literals; i++){
        if((literals[i] < 0) ? !((var_arr->vars)[-literals[i]]) : (var_arr->vars)[literals[i]]){
            return false;
        }
    }

    return true;
}