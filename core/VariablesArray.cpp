//
// Created by Xandru Mifsud on 27/09/2021.
//

#include "VariablesArray.h"

// Constructor for a VariablesArray instance, which initialises a boolean array of size n_vars with randomly generated
// boolean values.
VariablesArray::VariablesArray(uint32_t n_vars){
    default_random_engine engine(std::random_device{}()); // Get the system default random generator with a random seed
    RBG<default_random_engine> rbg(engine); // Initialise an instance of a random boolean generator...

    this->n_vars = n_vars;
    vars = new bool[n_vars];

    // Populate each entry in the vars array with a randomly generated boolean value...
    for(uint32_t i = 0; i < n_vars; i++){
        vars[i] = rbg.sample();
    }
}
