//
// Created by Xandru Mifsud on 28/09/2021.
//

#include "SATInstance.h"

/* Constructor for a SATInstance object, with the file path to a DIMACS formatted .cnf file accepted as input.
 * This constructor is responsible for:
 *  i. Loading meta data (number of variables, clauses etc) as well as reading the clauses from the specified .cnf file,
 *     using the CNF_IO package (available at https://people.sc.fsu.edu/~jburkardt/cpp_src/cnf_io/cnf_io.html)
 * ii. Encoding literals such that a variable x, where x is a non-negative integer, is mapped to 2x while its negation
 *     is mapped to 2x + 1. For this encoding, we can obtain the variable associated with a literal by a single left
 *     shift i.e. the variable v associated with a literal l is: v = l >> 1.
 */
SATInstance::SATInstance(const string& cnf_file_name, vector<Clause*>* clauses){
    int l_num, c_num, v_num;
    int* l_c_num;
    int* l_val;

    // Read the meta data from the .cnf file
    if(cnf_header_read(cnf_file_name, &v_num, &c_num, &l_num)){
        cout << "The header information could not be read. Exiting..." << endl;
        exit(1);
    }

    l_c_num = new int[c_num];
    l_val = new int[l_num];

    // Read the clauses data from the .cnf file
    cnf_data_read(cnf_file_name, v_num, c_num, l_num, l_c_num, l_val);

    int c, l, l_c;

    l = 0;
    for(c = 0; c < c_num; c++){ // For each read clauses, create a new Clause instance with its literals encoded in the
                                // aforementioned manner
        // l_c_num[c] = # of signed literals in clause c
        ull* literals; literals = new ull[l_c_num[c]];

        for(l_c = 0; l_c < l_c_num[c]; l_c++){ // For every literal in the clause...
            /* Note that the DIMACS format uses 0 as a special character, hence the variables are labelled as positive
             * (and hence non-zero) integers; Since we wish to use the variable x as an index to an array (where array
             * index starts from zero), we wish to shift the DIMACS variable label x to x - 1. Hence the DIMACS variable
             * x is encoded as 2x - 2. The negation of a variable x in DIMACS is the negative integer -x, hence applying
             * the shift and encoding, it is mapped to -2x - 1.
             */
            literals[l_c] = (0 < l_val[l]) ? (2 * l_val[l]) - 2 : ((-2) * l_val[l]) - 1;

            l += 1;
        }

        clauses->push_back(new Clause(literals, l_c_num[c]));
    }

    // Initialise internal state variables...
    n_vars = v_num; // Number of variables in SAT instance
    n_clauses = c_num; // Number of clauses in SAT instance
    n_literals = l_num; // Number of literals in SAT instance
    var_arr = new VariablesArray(n_vars);
}

// Given two Clause instances, establishes whether they are dependent or not. The criteria of dependency is having one
// or more variables in common between the literals of the two clauses.
bool SATInstance::dependent_clauses(Clause* c1, Clause* c2){
    // For every pair (x, y) of literals between the two clauses...
    for(int i = 0; i < c1->n_literals; i++){
        for(int j = 0; j < c2->n_literals; j++){
            // If the variable associated with x is equal to the variable associated with y, then the clauses are
            // dependent and hence we return true (no need to check further - we require AT LEAST one).
            if(((c1->literals)[i] >> 1) == ((c2->literals)[j] >> 1)){
                return true;
            }
        }
    }

    // Otherwise, we return false (independent if FOR ALL literal pairs, the underlying variables are distinct)
    return false;
}

vector<vector<Clause*>*>* SATInstance::getDependencyGraphComponents(vector<Clause*>* clauses){
    ull n_clauses = clauses->size();

    auto component = new vector<Clause*>;
    auto component_clauses = new vector<vector<Clause*>*>;

    vector<ull> neighbours;
    vector<ull> neighbours_queue;
    vector<ull> remaining_clauses;
    for(ull u = 0; u < n_clauses; u++){ remaining_clauses.push_back(u);}

    while(!remaining_clauses.empty()){
        if(neighbours.empty()){
            neighbours.push_back(remaining_clauses.front());
            component->push_back(clauses->at(remaining_clauses.front()));

            remaining_clauses.erase(remaining_clauses.begin());
        }

        for(auto n: neighbours){
            auto idx = remaining_clauses.begin();
            while(idx != remaining_clauses.end()){
                if(dependent_clauses(clauses->at(n), clauses->at(*idx))){
                    (clauses->at(n))->degree += 1;
                    (clauses->at(*idx))->degree += 1;

                    neighbours_queue.push_back(*idx);
                    component->push_back(clauses->at(*idx));

                    idx = remaining_clauses.erase(idx);
                }
                else{
                    idx++;
                }
            }
        }

        if(neighbours_queue.empty()){
            neighbours.clear();

            component_clauses->push_back(component);
            component = new vector<Clause*>;
        }
        else{
            neighbours = neighbours_queue;
            neighbours_queue.clear();
        }
    }

    if(!component->empty()){
        component_clauses->push_back(component);
    }
    return component_clauses;
}

// SAT solver based on the Algorithmic Lovasz Local Lemma of Moser and Tardos (2010)
VariablesArray* SATInstance::solve(vector<SubSATInstance*>* subInstances, bool parallel) const{
    #pragma omp parallel for if(parallel) schedule(dynamic) default(none) shared(subInstances)
    for(ull i = 0; i < subInstances->size(); i++){
        subInstances->at(i)->solve();
    }

    return var_arr;
}

vector<SubSATInstance*>* SATInstance::createSubSATInstances(vector<vector<Clause*>*>* components, ull parallel_resample) const{
    auto subInstance = new vector<SubSATInstance*>;
    for(auto c: *components){
        subInstance->push_back(new SubSATInstance(var_arr, c, parallel_resample));
    }

    return subInstance;
}

bool SATInstance::verify_validity(vector<Clause*>* clauses) const{
    for(auto c: *clauses){
        if(c->is_not_satisfied(var_arr)){
            return false;
        }
    }

    return true;
}