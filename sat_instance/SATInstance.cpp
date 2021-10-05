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

    // Display the meta data for monitoring purposes
    cout << "V_NUM = " << v_num << endl;
    cout << "C_NUM = " << c_num << endl;
    cout << "L_NUM = " << l_num << endl;

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
            if(((c1->literals)[i] >> 1) == ((c2->literals)[i] >> 1)){
                return true;
            }
        }
    }

    // Otherwise, we return false (independent if FOR ALL literal pairs, the underlying variables are distinct)
    return false;
}

// Utility function which checks the dependencies between the clauses of the SAT instance and constructs the associated
// Laplacian matrix L = D - A for each component, as an Eigen::MatrixXd instance, while also returning the vertex sets
// of the components of the dependency graph.
pair<vector<MatrixXd*>*, vector<vector<Clause*>*>*> SATInstance::getDependencyGraph(vector<Clause*>* clauses){
    ull n_clauses = clauses->size();

    // Initialise a MatrixXd instance and initialise all the entries to zero.
    auto laplacian = new MatrixXd(n_clauses, n_clauses);
    laplacian->setZero();

    auto component_laplacians = new vector<MatrixXd*>;
    auto component_clauses = new vector<vector<Clause*>*>;

    // Initialise a vector with the vertices (clauses) of the dependency graph the remain (initially all n_clauses)
    vector<ull> remaining_vertices;
    for(ull u = 0; u < n_clauses; u++){ remaining_vertices.push_back(u); }
    auto v = remaining_vertices.begin(); // Iterator over the remaining_vertices vector...

    // Lambda expression which recursively constructs the Laplacian and componets of the graph
    function<void(vector<ull>*)> _get_component;
    _get_component = [&](vector<ull>* component){
        ull i = *v; // Current vertex (clause) pointed to by iterator v
        component->push_back(i); // Add to current component
        remaining_vertices.erase(v); // Remove from remaining_vertices vector

        // Iterating over all the remaining vertices, find all the neighbours of v
        for(auto u = remaining_vertices.begin(); u != remaining_vertices.end();){
            ull j = *u; // Current vertex (clause) pointed to by iterator u

            // If the clauses i and j are dependent, then:
            if(dependent_clauses(clauses->at(i), clauses->at(j))){
                (*laplacian)(i, j) = -1; // The entry (i, j) of the Laplacian is -1
                (*laplacian)(j, i) = -1; // The entry (j, i) of the Laplacian is -1 by symmetry
                (*laplacian)(i, i) += 1; // The degree of i, i.e. the entry (i, i) of the Laplacian, increases by 1
                (*laplacian)(j, j) += 1; // The degree of j, i.e. the entry (j, j) of the Laplacian, increases by 1
            } // Otherwise if independent, the entries (i, j) and (j, i) of the Laplacian are 0

            ++u; // Increment iterator
        }

        // Recursively call _get_component() on the neighbours of the vertex pointed to by the iterator v
        for(; v != remaining_vertices.end();){
            ull j = *v;

            if((*laplacian)(i, j) == -1){
                _get_component(component);
            }
            else{
                ++v;
            }
        }

        v = remaining_vertices.begin(); // Reset to beginning (otherwise some vertices will not be visited)
    };

    // Until remaining_vertices is empty, add a new component and call the _get_component() lambda function to find the
    // next component and its adjacencies
    while(v != remaining_vertices.end()){
        auto component = new vector<ull>;
        _get_component(component);

        auto curr_laplacian = new MatrixXd(component->size(), component->size());
        *curr_laplacian = (*laplacian)(*component, *component);
        component_laplacians->push_back(curr_laplacian);

        auto curr_clauses = new vector<Clause*>;
        for(auto c: *component){
            curr_clauses->push_back(clauses->at(c));
        }

        component_clauses->push_back(curr_clauses);
        delete component;
    }

    delete laplacian;

    return {component_laplacians, component_clauses};
}

// SAT solver based on the Algorithmic Lovasz Local Lemma of Moser and Tardos (2010)
VariablesArray* SATInstance::solve(vector<SubSATInstance*>* subInstances) const{
    for(auto sat: *subInstances){
        sat->solve();
    }

    return var_arr;
}

vector<SubSATInstance*>* SATInstance::createSubSATInstances(vector<vector<Clause*>*>* components) const{
    auto subInstance = new vector<SubSATInstance*>;
    for(auto c: *components){
        subInstance->push_back(new SubSATInstance(var_arr, c));
    }

    return subInstance;
}