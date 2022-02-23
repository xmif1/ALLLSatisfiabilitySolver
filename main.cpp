#include <iostream>
#include <fstream>
#include <cstring>
#include <chrono>
#include <ctime>

#include <omp.h>

#include "sat_instance/SATInstance.h"

void output(const string& str, ofstream* out_f, bool dump);

int main(int argc, char *argv[]){
    int n_threads = 0;
    bool dump = false;
    string cnf_fpath;

    // Option checking...
    if(argc <= 1){
        throw std::runtime_error("The filepath to the DIMACS formatted CNF SAT instance was not specified...exiting...");
    }
    else if(argc > 5){
        throw std::runtime_error("Too many parameters specified...exiting...");
    }
    else if(argc == 3 && strcmp(argv[1], "-o") == 0){
        cnf_fpath = argv[2];
        dump = true;
    }
    else if(argc == 3 && strcmp(argv[1], "-p") == 0){
        cnf_fpath = argv[2];
        n_threads = omp_get_num_procs();
        if(n_threads < 2){
            n_threads = 0; // Do not parallelise
        }
    }
    else if(argc == 4 && strcmp(argv[1], "-p") == 0 && strcmp(argv[2], "-o") != 0){
        cnf_fpath = argv[3];
        n_threads = stoi(argv[2]);
        if(n_threads > omp_get_num_procs() || n_threads < 2){
            throw std::runtime_error("Invalid number of threads specified...exiting...");
        }
    }
    else if(argc == 4 && ((strcmp(argv[1], "-o") == 0 && (strcmp(argv[2], "-p")) == 0) ||
                          (strcmp(argv[1], "-p") == 0 && (strcmp(argv[2], "-o")) == 0))){
        cnf_fpath = argv[3];
        dump = true;

        n_threads = omp_get_num_procs();
        if(n_threads < 2){
            n_threads = 0; // Do not parallelise
        }
    }
    else if(argc == 5 && (strcmp(argv[1], "-o") == 0 && (strcmp(argv[2], "-p")) == 0)){
        cnf_fpath = argv[4];
        dump = true;

        n_threads = stoi(argv[3]);
        if(n_threads > omp_get_num_procs() || n_threads < 2){
            throw std::runtime_error("Invalid number of threads specified...exiting...");
        }
    }
    else if(argc == 5 && (strcmp(argv[1], "-p") == 0 && (strcmp(argv[3], "-l")) == 0)){
        cnf_fpath = argv[4];
        dump = true;

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

    ofstream* out_f;
    ofstream* stat_f;

    if(dump){
        string out_fpath = cnf_fpath; out_fpath.replace(out_fpath.size() - 4, 4, ".out");
        out_f = new ofstream(out_fpath);

        string stat_fpath = cnf_fpath; stat_fpath.replace(stat_fpath.size() - 4, 4, ".csv");
        stat_f = new ofstream(stat_fpath);
    }
    else{
        out_f = nullptr;
        stat_f = nullptr;
    }

    // Reading SAT instance from .cnf file...
    // Logging...
    auto timestart = chrono::system_clock::to_time_t(chrono::system_clock::now());
    auto time_str = string(ctime(&timestart)); time_str.pop_back();
    string log_start = "Log "; output(log_start.append(time_str) + ": Reading CNF file\n", out_f, dump);
    auto start = chrono::high_resolution_clock::now();

    // Initialise new SATInstance from specified CNF file
    auto satInstance = new SATInstance(cnf_fpath, n_threads);

    auto stop = chrono::high_resolution_clock::now();
    auto read_duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    auto timeend = chrono::system_clock::to_time_t(chrono::system_clock::now());
    time_str = string(ctime(&timeend)); time_str.pop_back();
    string log_end = "Log "; output(log_end.append(time_str) + ": Read complete; Duration: " +
                                    to_string(read_duration.count() / 1000.0) + "s\n\n", out_f, dump);

    // Display the meta data for monitoring purposes
    output("------------ INFORMATION ------------\n# Variables\t= " + to_string(satInstance->n_vars) +
    "\n# Clauses\t= " + to_string(satInstance->n_clauses) +
    "\n# Literals\t= " + to_string(satInstance->n_literals) +
    "\n-------------------------------------\n\n", out_f, dump);

    if(dump){
        *stat_f << to_string(read_duration.count() / 1000.0) + ",";
        *stat_f << to_string(satInstance->n_vars) + ",";
        *stat_f << to_string(satInstance->n_clauses) + ",";
        *stat_f << to_string(satInstance->n_literals) + ",";
    }


    string solve_info = "Starting sequential solve (# Threads = 1)";
    if(n_threads){
        solve_info = "Starting parallel solve (# Threads = " + to_string(n_threads) + ")";
    }

    // Solve the SAT instance (using the Algorithmic Lovasz Local Lemma)
    // Logging...
    timestart = chrono::system_clock::to_time_t(chrono::system_clock::now());
    time_str = string(ctime(&timestart)); time_str.pop_back();
    log_start = "Log "; output(log_start.append(time_str) + ": " + solve_info + "\n", out_f, dump);
    start = chrono::high_resolution_clock::now();

    Statistics* statistics = satInstance->solve();

    stop = chrono::high_resolution_clock::now();
    auto solve_duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    timeend = chrono::system_clock::to_time_t(chrono::system_clock::now());
    time_str = string(ctime(&timeend)); time_str.pop_back();
    log_end = "Log "; output(log_end.append(time_str) + ": Completed solve; Duration: " +
                             to_string(solve_duration.count() / 1000.0) + "s\n\n", out_f, dump);

    if(dump){
        *stat_f << to_string(solve_duration.count() / 1000.0) + ",";
    }

    if(n_threads){
        output("------------ STATISTICS -------------\n# Iterations\t= " + to_string(statistics->n_iterations) +
        "\n# Resamples\t= " + to_string(statistics->n_resamples), out_f, dump);

        for(int t = 0; t < n_threads; t++){
            output("\n\tThread " + to_string(t+1) + ": " +
            to_string((statistics->n_thread_resamples).at(t)), out_f, dump);
        }

        output("\n\nAvg. UNSAT MIS Size = " + to_string(statistics->avg_mis_size) +
        "\n-------------------------------------\n\n", out_f, dump);

        if(dump){
            *stat_f << to_string(n_threads) + ",";
            *stat_f << to_string(statistics->n_iterations) + "\n";
        }
    }
    else{
        output("------------ STATISTICS -------------\n# Iterations\t= " + to_string(statistics->n_iterations) +
               "\n# Resamples\t= " + to_string(statistics->n_resamples) +
               "\n-------------------------------------\n\n", out_f, dump);

        if(dump){
            *stat_f << "1,";
            *stat_f << to_string(statistics->n_iterations) + "\n";
        }
    }

    if(satInstance->verify_validity()){
        // Print the variable assignment for the solution
        output("SATISFIABLE\n", out_f, dump);

        if(dump){
            for(ull i = 0; i < satInstance->n_vars; i++){
                *out_f << "\nVariable " + to_string(i + 1) + " = " + to_string((satInstance->var_arr->vars)[i]);
            }

            stat_f->close();
            out_f->close();
        }

        return 0;
    }
    else{
        output("ERROR: Solver converged to an invalid solution!\n", out_f, dump);

        if(dump){
            stat_f->close();
            out_f->close();
        }

        return 1;
    }

}

void output(const string& str, ofstream* out_f, bool dump){
    cout << str;
    if(dump){ *out_f << str;}
}