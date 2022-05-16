#include <iostream>
#include <fstream>
#include <cstring>
#include <chrono>
#include <ctime>

#include <omp.h>

#include "sat_instance/SATInstance.h"

void output(const string& str, ofstream* out_f, bool dump);

// =====================================================================================================================
// ---------------------------------------------------- INFORMATION ----------------------------------------------------
// =====================================================================================================================

/* OPTIONS: [-o] [-p [n_threads]] <cnf_file_path>
 *
 * INPUT: SAT CNF instance in DIMACS format at <cnf_file_path>
 *
 * OUTPUT:
 *      Terminal: Instance meta-data and time-stamped solver logging; if solution found, 'SATISFIABLE' is printed along
 *                with some solve statistics, including running time, no. of iterations and no. of variable resamples.
 *      If output flag -o is specified:
 *          .out file: Copy of terminal output; if satisfiable, the variable assignment is also outputted to this file.
 *          .csv file: Solve statistics structured as follows:  [0] instance read/load time
 *                                                              [1] number of variables
 *                                                              [2] number of clauses
 *                                                              [3] number of literals
 *                                                              [4] instance solve time
 *                                                              [5] number of threads
 *                                                              [6] number of solve iterations
 */

int main(int argc, char *argv[]){
    int n_threads = 0; // number of threads (will be 0 if sequential, > 1 if parallel)
    bool dump = false; // flag signally whether solver meta-data and statistics is to be dumped to a text file
    string cnf_fpath;  // path to CNF instance

    // =================================================================================================================
    // ------------------------------------------------ OPTION PARSING -------------------------------------------------
    // =================================================================================================================

    /* OPTIONS: [-o] [-p [n_threads]] <cnf_file_path>
     * where: i. -o is an argument signalling output of solver meta-data and statistics to text files
     *       ii. -p is an argument signalling parallel solving; if the optional int value 'n_threads' is not passed, the
     *           solver will use all available (physical) processor threads; otherwise, 'n_threads' are used.
     */
    if(argc <= 1){ // if file path to a CNF instance not given, throw runtime error
        throw std::runtime_error("The filepath to the DIMACS formatted CNF SAT instance was not specified...exiting...");
    }
    else if(argc > 5){ // else if too many options given, throw runtime error
        throw std::runtime_error("Too many parameters specified...exiting...");
    }
    else if(argc == 3 && strcmp(argv[1], "-o") == 0){ // if output flag -o passed
        cnf_fpath = argv[2];
        dump = true;
    }
    else if(argc == 3 && strcmp(argv[1], "-p") == 0){ // if parallelisation flag -p passed (but n_threads not specified)
        cnf_fpath = argv[2];
        n_threads = omp_get_num_procs(); // get number of (physical) threads available
        if(n_threads < 2){ // if number of physical threads is 1
            n_threads = 0; // then do not parallelise and use sequential solver
        }
    }
    // if parallelisation flag -p passed and n_threads specified
    else if(argc == 4 && strcmp(argv[1], "-p") == 0 && strcmp(argv[2], "-o") != 0){
        cnf_fpath = argv[3];
        n_threads = stoi(argv[2]);
        if(n_threads > omp_get_num_procs() || n_threads < 2){ // check that n_threads is a valid value, else throw error
            throw std::runtime_error("Invalid number of threads specified...exiting...");
        }
    }
    // if parallisation flag -p passed (but n_threads not specified) and output flag -o passed (in any order)
    else if(argc == 4 && ((strcmp(argv[1], "-o") == 0 && (strcmp(argv[2], "-p")) == 0) ||
                          (strcmp(argv[1], "-p") == 0 && (strcmp(argv[2], "-o")) == 0))){
        cnf_fpath = argv[3];
        dump = true;

        n_threads = omp_get_num_procs(); // get number of (physical) threads available
        if(n_threads < 2){ // if number of physical threads is 1
            n_threads = 0; // then do not parallelise and use sequential solver
        }
    }
    // if the output flag -o is passed, followed by parallisation flag -p passed and n_threads specified
    else if(argc == 5 && (strcmp(argv[1], "-o") == 0 && (strcmp(argv[2], "-p")) == 0)){
        cnf_fpath = argv[4];
        dump = true;

        n_threads = stoi(argv[3]);
        if(n_threads > omp_get_num_procs() || n_threads < 2){ // check that n_threads is a valid value, else throw error
            throw std::runtime_error("Invalid number of threads specified...exiting...");
        }
    }
    // if parallisation flag -p passed and n_threads specified, followed by the output flag -o
    else if(argc == 5 && (strcmp(argv[1], "-p") == 0 && (strcmp(argv[3], "-o")) == 0)){
        cnf_fpath = argv[4];
        dump = true;

        n_threads = stoi(argv[2]);
        if(n_threads > omp_get_num_procs() || n_threads < 2){ // check that n_threads is a valid value, else throw error
            throw std::runtime_error("Invalid number of threads specified...exiting...");
        }
    }
    else if(argc == 2){ // if no options passed, only path to CNF instance
        cnf_fpath = argv[1];
    }
    else{ // otherwise throw a runtime error if any other argument combination is given
        throw std::runtime_error("Invalid options specified...exiting...");
    }

    // =================================================================================================================
    // -------------------------------------------- LOGGING FILE STRUCTURES --------------------------------------------
    // =================================================================================================================

    ofstream* out_f;  // for solver logging etc as well as final solution variable assignment
    ofstream* stat_f; // for solve statistics

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

    // =================================================================================================================
    // --------------------------------------------- LOADING SAT INSTANCE ----------------------------------------------
    // =================================================================================================================

    // Reading SAT instance from .cnf file
    // Logging read start...
    auto timestart = chrono::system_clock::to_time_t(chrono::system_clock::now());
    auto time_str = string(ctime(&timestart)); time_str.pop_back();
    string log_start = "Log "; output(log_start.append(time_str) + ": Reading CNF file\n", out_f, dump);
    auto start = chrono::high_resolution_clock::now();

    // Initialise new SATInstance from specified CNF file
    auto satInstance = new SATInstance(cnf_fpath, n_threads);

    // Logging read complete...
    auto stop = chrono::high_resolution_clock::now();
    auto read_duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    auto timeend = chrono::system_clock::to_time_t(chrono::system_clock::now());
    time_str = string(ctime(&timeend)); time_str.pop_back();
    string log_end = "Log "; output(log_end.append(time_str) + ": Read complete; Duration: " +
                                    to_string(read_duration.count() / 1000.0) + "s\n\n", out_f, dump);

    // Display the SAT instance meta-data for monitoring purposes
    output("------------ INFORMATION ------------\n# Variables\t= " + to_string(satInstance->n_vars) +
    "\n# Clauses\t= " + to_string(satInstance->n_clauses) +
    "\n# Literals\t= " + to_string(satInstance->n_literals) +
    "\n-------------------------------------\n\n", out_f, dump);

    if(dump){ // Save statistics to CSV file
        *stat_f << to_string(read_duration.count() / 1000.0) + ",";
        *stat_f << to_string(satInstance->n_vars) + ",";
        *stat_f << to_string(satInstance->n_clauses) + ",";
        *stat_f << to_string(satInstance->n_literals) + ",";
    }

    // =================================================================================================================
    // --------------------------------------------- SOLVING SAT INSTANCE ----------------------------------------------
    // =================================================================================================================

    string solve_info = "Starting sequential solve (# Threads = 1)";
    if(n_threads){
        solve_info = "Starting parallel solve (# Threads = " + to_string(n_threads) + ")";
    }

    // Solve the SAT instance (using the Algorithmic Lovasz Local Lemma)
    // Logging solve start...
    timestart = chrono::system_clock::to_time_t(chrono::system_clock::now());
    time_str = string(ctime(&timestart)); time_str.pop_back();
    log_start = "Log "; output(log_start.append(time_str) + ": " + solve_info + "\n", out_f, dump);
    start = chrono::high_resolution_clock::now();

    Statistics* statistics = satInstance->solve(); // Solving SAT instance and maintaining solve statistics

    // Logging solve completion...
    stop = chrono::high_resolution_clock::now();
    auto solve_duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    timeend = chrono::system_clock::to_time_t(chrono::system_clock::now());
    time_str = string(ctime(&timeend)); time_str.pop_back();
    log_end = "Log "; output(log_end.append(time_str) + ": Completed solve; Duration: " +
                             to_string(solve_duration.count() / 1000.0) + "s\n\n", out_f, dump);

    if(dump){ // Add solve time (in milliseconds) to statistics CSV
        *stat_f << to_string(solve_duration.count()) + ",";
    }

    // =================================================================================================================
    // ------------------------------------------- LOGGING SOLVE STATISTICS --------------------------------------------
    // =================================================================================================================

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

    // =================================================================================================================
    // ---------------------------------------------- VERIFYING SOLUTION -----------------------------------------------
    // =================================================================================================================

    if(satInstance->verify_validity()){ // If solution is valid
        // Print that instance is SATISFIABLE
        output("SATISFIABLE\n", out_f, dump);

        if(dump){
            // Output solution variable assignment to .out file (and NOT to console -- handy esp. for large instances)
            for(ull i = 0; i < satInstance->n_vars; i++){
                *out_f << "\nVariable " + to_string(i + 1) + " = " + to_string((satInstance->var_arr->vars)[i]);
            }

            // Close file pointers
            stat_f->close();
            out_f->close();
        }

        return 0;
    }
    else{ // If converged to an invalid solution, report so and terminate erroneously
        output("ERROR: Solver converged to an invalid solution!\n", out_f, dump);

        if(dump){
            // Close file pointers
            stat_f->close();
            out_f->close();
        }

        return 1;
    }

}

// Convenience function for duplicating the console output to file (if desired)
void output(const string& str, ofstream* out_f, bool dump){
    cout << str;
    if(dump){ *out_f << str;}
}