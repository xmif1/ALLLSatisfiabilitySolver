#include <iostream>
#include <chrono>
#include <ctime>

#include "sat_instance/SATInstance.h"
#include <omp.h>

void output(const string& str, ofstream& out_f);

int main(int argc, char *argv[]){
    bool parallel = false;
    int n_threads = 0;

    string cnf_fpath;

    // Option checking...
    if(argc <= 1){
        throw std::runtime_error("The filepath to the DIMACS formatted CNF SAT instance was not specified...exiting...");
    }
    else if(argc > 4){
        throw std::runtime_error("Too many parameters specified...exiting...");
    }
    else if(argc == 3 && strcmp(argv[1], "-p") == 0){
        cnf_fpath = argv[2];
        parallel = true;
    }
    else if(argc == 4 && strcmp(argv[1], "-p") == 0){
        cnf_fpath = argv[3];
        parallel = true;
        n_threads = stoi(argv[2]);
        if(n_threads > omp_get_num_procs() || n_threads < 2){
            throw std::runtime_error("Invalid number of threads specified...exiting...");
        }
    }
    else if(argc == 2){
        cnf_fpath = argv[1];
    }
    else{
        throw std::runtime_error("Invalid options specified...exiting...");
    }

    string out_fpath = cnf_fpath; out_fpath.replace(out_fpath.size() - 4, 4, ".out");
    ofstream out_f(out_fpath);

    // Initialise new SATInstance from specified CNF file
    auto clauses = new vector<Clause*>;
    auto satInstance = new SATInstance(cnf_fpath, clauses);

    // Display the meta data for monitoring purposes
    output("V_NUM = " + to_string(satInstance->n_vars) + "\nC_NUM = " + to_string(satInstance->n_clauses) +
           "\nL_NUM = " + to_string(satInstance->n_literals) + "\n\n", out_f);

    auto components = new vector<vector<Clause*>*>;
    components->push_back(clauses);

    auto subSATInstances = satInstance->createSubSATInstances(components, n_threads);
    output("# of components = " + to_string(subSATInstances->size()) + "\n\n", out_f);

    for(ull i = 0; i < subSATInstances->size(); i++){
        string alll_str;
        if(parallel && (subSATInstances->at(i))->is_ALLL_compatible()){
            alll_str = " (ALLL Compatible)";
        }

        output("Component " + to_string(i + 1) + alll_str + ": # of clauses = " +
               to_string((subSATInstances->at(i))->clauses->size()) + "\n", out_f);
    }

    output("\n", out_f);

    // Solve the SAT instance (using the Algorithmic Lovasz Local Lemma)
    // Logging...
    auto timestart = chrono::system_clock::to_time_t(chrono::system_clock::now());
    string log_start = "Log "; output(log_start.append(ctime(&timestart)) + "\tStarting solve...\n", out_f);
    auto start = chrono::high_resolution_clock::now();

    VariablesArray* sat = satInstance->solve(subSATInstances, parallel);

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    auto timeend = chrono::system_clock::to_time_t(chrono::system_clock::now());
    string log_end = "Log "; output(log_end.append(ctime(&timeend)) + "\tCompleted solve...Duration: " +
                                    to_string(duration.count()) + "\n\n", out_f);

    if(satInstance->verify_validity(clauses)){
        // Print the variable assignment for the solution
        output("SATISFIABLE\n", out_f);

        out_f.close();
        return 0;
    }
    else{
        output("ERROR: Solver converged to an invalid solution!\n", out_f);

        out_f.close();
        return 1;
    }

}

void output(const string& str, ofstream& out_f){
    cout << str;
    out_f << str;
}