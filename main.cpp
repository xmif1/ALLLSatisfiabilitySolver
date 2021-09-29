#include "sat_instance/SATInstance.h"

int main(int argc, char *argv[]){
    string cnf_fpath;

    // option checking...
    if(argc <= 1){
        throw std::runtime_error("The filepath to the DIMACS formatted CNF SAT instance was not specified...exiting...");
    }
    else if(argc > 2){
        throw std::runtime_error("Too many parameters specified...exiting...");
    }
    else{
        cnf_fpath = argv[1];
    }

    auto satInstance = new SATInstance(cnf_fpath);
    VariablesArray* sat = satInstance->solve();

    for(ull i = 0; i < satInstance->n_vars; i++){
        cout << "Var" << i + 1 << " = " << (sat->vars)[i] << endl;
    }

    return 0;
}
