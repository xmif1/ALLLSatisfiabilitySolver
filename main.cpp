#include <iostream>
#include <fstream>
#include <chrono>
#include <ctime>

#include <omp.h>

#include "sat_instance/SATInstance.h"

void output(const string& str, ofstream& out_f);

int main(int argc, char *argv[]){
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
        n_threads = omp_get_num_procs();
        if(n_threads < 2){
            n_threads = 0; // Do not parallelise
        }
    }
    else if(argc == 4 && strcmp(argv[1], "-p") == 0){
        cnf_fpath = argv[3];
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

    // Reading SAT instance from .cnf file...
    // Logging...
    auto timestart = chrono::system_clock::to_time_t(chrono::system_clock::now());
    auto time_str = string(ctime(&timestart)); time_str.pop_back();
    string log_start = "Log "; output(log_start.append(time_str) + ": Reading CNF file\n", out_f);
    auto start = chrono::high_resolution_clock::now();

    // Initialise new SATInstance from specified CNF file
    auto satInstance = new SATInstance(cnf_fpath, n_threads);

    auto stop = chrono::high_resolution_clock::now();
    auto read_duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    auto timeend = chrono::system_clock::to_time_t(chrono::system_clock::now());
    time_str = string(ctime(&timeend)); time_str.pop_back();
    string log_end = "Log "; output(log_end.append(time_str) + ": Read complete; Duration: " +
                                    to_string(read_duration.count() / 1000.0) + "s\n\n", out_f);

    // Display the meta data for monitoring purposes
    output("---------- INFORMATION ----------\n# Variables\t= " + to_string(satInstance->n_vars) +
    "\n# Clauses\t= " + to_string(satInstance->n_clauses) +
    "\n# Literals\t= " + to_string(satInstance->n_literals) +
    "\n---------------------------------\n\n", out_f);


    string solve_info = "Starting sequential solve (# Threads = 1)";
    if(n_threads){
        solve_info = "Starting parallel solve (# Threads = " + to_string(n_threads) + ")";
    }

    // Solve the SAT instance (using the Algorithmic Lovasz Local Lemma)
    // Logging...
    timestart = chrono::system_clock::to_time_t(chrono::system_clock::now());
    time_str = string(ctime(&timestart)); time_str.pop_back();
    log_start = "Log "; output(log_start.append(time_str) + ": " + solve_info + "\n", out_f);
    start = chrono::high_resolution_clock::now();

    VariablesArray* sat = satInstance->solve();

    stop = chrono::high_resolution_clock::now();
    auto solve_duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    timeend = chrono::system_clock::to_time_t(chrono::system_clock::now());
    time_str = string(ctime(&timeend)); time_str.pop_back();
    log_end = "Log "; output(log_end.append(time_str) + ": Completed solve; Duration: " +
                             to_string(solve_duration.count() / 1000.0) + "s\n\n", out_f);

    if(satInstance->verify_validity()){
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