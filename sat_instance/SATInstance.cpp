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
SATInstance::SATInstance(const string& cnf_file_name){
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

        clauses.push_back(new Clause(literals, l_c_num[c]));
    }

    // Initialise internal state variables...
    n_vars = v_num; // Number of variables in SAT instance
    n_clauses = c_num; // Number of clauses in SAT instance
    C = P_2e64_m59 % n_clauses; // Starting seed for an LCG based index to the this.clauses vector, which pseudo-randomly
                                // visits all the clauses exactly once with a period of n_clauses, allowing for efficient
                                // 'shuffling' of the clauses vector.
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
// Laplacian matrix L = D - A, as an Eigen::MatrixXd instance of dimension n_clauses by n_clauses, while also returning
// the vertex sets of the components of the dependency graph.
pair<MatrixXd*, vector<vector<ull>*>*> SATInstance::getDependencyGraph(){
    // Initialise a MatrixXd instance and initialise all the entries to zero.
    auto laplacian = new MatrixXd(n_clauses, n_clauses);
    laplacian->setZero();

    // Initialise a vector of ull vectors (each such vector being a vertex set of some component)
    auto components = new vector<vector<ull>*>;

    // Initialise a vector with the vertices (clauses) of the dependency graph the remain (initially all n_clauses)
    vector<ull> remaining_vertices;
    for(ull u = 0; u < n_clauses; u++){ remaining_vertices.push_back(u); }
    auto v = remaining_vertices.begin(); // Iterator over the remaining_vertices vector...

    // Lambda expression which recursively constructs the Laplacian and componets of the graph
    function<void()> _get_component;
    _get_component = [&](){
        ull i = *v; // Current vertex (clause) pointed to by iterator v
        (components->back())->push_back(i); // Add to current component
        remaining_vertices.erase(v); // Remove from remaining_vertices vector

        // Iterating over all the remaining vertices, find all the neighbours of v
        for(auto u = remaining_vertices.begin(); u != remaining_vertices.end();){
            ull j = *u; // Current vertex (clause) pointed to by iterator u

            // If the clauses i and j are dependent, then:
            if(dependent_clauses(clauses.at(i), clauses.at(j))){
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
                _get_component();
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
        components->push_back(new vector<ull>);
        _get_component();
    }

    return {laplacian, components};
}

/* Checks whether the SAT instance is satisfied by a given variable assignment (specified as a VariablesArray instance),
 * by checking if each Clause instance in the clauses vector is satisfied. We iterate over the clauses vector in a pseudo
 * random manner, by means of an LCG. Note that we maintain the state of the LCG as a class variable.
 *
 * In the case that a clause is not satisfied, we return a pointer to the Clause instance. Otherwise, if every clause is
 * satisfied i.e. the SAT instance is satisfied, we return a nullptr.
 */
Clause* SATInstance::is_satisfied(VariablesArray* var_arr){
    for(ull i = 0; i < n_clauses; i++){ // For each clause...
        // Fetch the clause specified at the (class variable) index C, and check if it is satisfied...
        if((clauses.at(C))->is_not_satisfied(var_arr)){ // If it is not satisfied, return a ptr to the Clause instance
            return clauses.at(C);
        }
        else{ /* Otherwise, calculate the next index C by means of the LCG; note that the LCG has a period of n_clauses
               * and hence we iterate over each clause exactly once, in a non-linear fashion; this is equivalent to a
               * pseudo-random shuffling of the clauses vector, and is required for optimal convergence of the Algorithmic
               * Lovasz Local Lemma method.
               */
            C = (C + P_2e64_m59) % n_clauses;
        }
    }

    return nullptr;
}

// SAT solver based on the Algorithmic Lovasz Local Lemma of Moser and Tardos (2010)
VariablesArray* SATInstance::solve(){
    default_random_engine engine(std::random_device{}()); // Get the system default random generator with a random seed
    RBG<default_random_engine> rbg(engine); // Initialise an instance of a random boolean generator...

    auto var_arr = new VariablesArray(n_vars);

    Clause* c = is_satisfied(var_arr);
    while(c){ // While there exists a clause c which is not satisfied...
        // For every variable in the clause (obtained by left shifting by 1 the literal encoding), randomly re-sample
        for(ull i = 0; i < c->n_literals; i++){
            (var_arr->vars)[(c->literals)[i] >> 1] = rbg.sample();
        }

        c = is_satisfied(var_arr);
    }

    sat = var_arr;

    return var_arr;
}