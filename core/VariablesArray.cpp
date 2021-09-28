//
// Created by Xandru Mifsud on 27/09/2021.
//

#include "VariablesArray.h"

VariablesArray::VariablesArray(ull n_vars){
    default_random_engine engine(std::random_device{}());
    RBG<default_random_engine> rbg(engine);

    this->n_vars = n_vars;
    vars = new bool[n_vars + 1];

    for(ull i = 1; i <= n_vars; i++){
        vars[i] = rbg.sample();
    }
}

bool& VariablesArray::operator[](lli index) const{
    return (index < 0) ? vars[-index] : vars[index];
}
