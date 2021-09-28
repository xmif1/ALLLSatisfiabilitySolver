//
// Created by Xandru Mifsud on 27/09/2021.
//

#include "Clause.h"

Clause::Clause(ull* literals, ull n_literals){
    this->literals = literals;
    this->n_literals = n_literals;
}

bool Clause::is_not_satisfied(VariablesArray* var_arr) const{
    for(ull i = 0; i < n_literals; i++){
        if((literals[i] & 1) ? !((var_arr->vars)[literals[i] >> 1]) : (var_arr->vars)[literals[i] >> 1]){
            return false;
        }
    }

    return true;
}