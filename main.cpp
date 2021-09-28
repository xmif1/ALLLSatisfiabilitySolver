#include "sat_instance/SATInstance.h"

int main(){
    auto satInstance = new SATInstance("/Users/xandrumifsud/PLLL/ALLLSatisfiabilitySolver/aim-50-1_6-yes1-4.cnf");
    VariablesArray* sat = satInstance->solve();

    for(ull i = 1; i <= satInstance->n_vars; i++){
        cout << "Var" << i << " = " << (sat->vars)[i] << endl;
    }

    return 0;
}
