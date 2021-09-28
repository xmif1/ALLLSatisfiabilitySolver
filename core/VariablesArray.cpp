//
// Created by Xandru Mifsud on 27/09/2021.
//

#include "VariablesArray.h"

VariablesArray::VariablesArray(ull n_vars){
    default_random_engine engine(std::random_device{}());
    RBG<default_random_engine> rbg(engine);

    this->n_vars = n_vars;
    vars = new bool[n_vars];

    for(ull i = 0; i < n_vars; i++){
        vars[i] = rbg.sample();
    }
}
