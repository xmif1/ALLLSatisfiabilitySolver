#include <iostream>
#include <chrono>
#include <ctime>

#include <boost/program_options.hpp>

#include "SATInstance.h"
#include "cnf_io/cnf_io.h"

using namespace boost::program_options;

typedef uint32_t UINT_T;
typedef SATInstance<UINT_T>::ClauseArray ClauseArray;

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
    int n_threads = 1;    // number of threads (will be 0 if sequential, > 1 if parallel)
    bool dump = false;    // flag signalling whether solver meta-data and statistics is to be dumped to a text file
    string cnf_fpath;     // path to CNF instance

    // =================================================================================================================
    // ------------------------------------------------ OPTION PARSING -------------------------------------------------
    // =================================================================================================================

    /* OPTIONS: [-h] [-o] [-p n_threads] --sat <cnf_file_path>
     * where: i. -o is an argument signalling output of solver meta-data and statistics to text files
     *       ii. -p is an argument signalling parallel solving; if the int value 'n_threads' is < 0, the
     *           solver will use all available (physical) processor threads; otherwise, 'n_threads' are used.
     *           Note: the number of threads used is min(n_threads, omp_get_num_procs())
     */

    try {
        options_description desc{"Options"};
        desc.add_options()
            ("help,h", "Help")
            ("output,o", "Output meta-data and statistics to separate files")
            ("parallel,p", value<int>()->default_value(0), "Use parallel solver")
            ("sat", value<string>()->required(), "Path to SAT instance in DIMACS-CNF format");

        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);

        if(vm.count("help")) {
            std::cout << desc << '\n';
            exit(0);
        } else {
            notify(vm);

            if(vm.count("output")){
                dump = true;
            }

            if(vm.count("parallel")){
                if (vm["parallel"].as<int>() < 0 || vm["parallel"].as<int>() > omp_get_num_procs()) {
                    n_threads = omp_get_num_procs();
                } else if (vm["parallel"].as<int>() > 0) {
                    n_threads = vm["parallel"].as<int>();
                } else {
                    n_threads = 1;
                }
            }

            if (vm.count("sat")) {
                cnf_fpath = vm["sat"].as<string>();
            }
        }
    } catch (const error &e) {
        cerr << e.what() << endl;
        exit(1);
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

    /* Loading a SAT instance is a two-step operation consisting of the following:
     *  i. Loading meta-data (number of variables, clauses etc) as well as reading the clauses from the specified .cnf file,
     *     using the CNF_IO package (available at https://people.sc.fsu.edu/~jburkardt/cpp_src/cnf_io/cnf_io.html)
     * ii. Encoding literals such that a variable x, where x is a non-negative integer, is mapped to 2x while its negation
     *     is mapped to 2x + 1. For this encoding, we can obtain the variable associated with a literal by a single left
     *     shift i.e. the variable v associated with a literal l is: v = l >> 1.
     */

    int c_num, v_num, l_num;
    int* l_c_num;
    int* l_val;

    // Read the meta-data from the .cnf file
    if(cnf_header_read(cnf_fpath, &v_num, &c_num, &l_num)){
        cout << "The header information could not be read. Exiting..." << endl;
        exit(1);
    }

    l_c_num = new int[c_num];
    l_val = new int[l_num];

    // Read the clause data from the .cnf file
    cnf_data_read(cnf_fpath, v_num, c_num, l_num, l_c_num, l_val);

    int chunk_size = ceil(c_num / (double) n_threads);

    auto clauses = new vector<ClauseArray*>();
    for(int t = 0; t < n_threads; t++){
        clauses->push_back(new ClauseArray());
    }

    int c, l, l_c; l = 0; unsigned short int t = 0;
    for(c = 0; c < c_num; c++){ // For each read clauses, create a new Clause instance with its literals encoded in the
        // aforementioned manner
        // l_c_num[c] = # of signed literals in clause c
        auto literals = new vector<UINT_T>;
        for(l_c = 0; l_c < l_c_num[c]; l_c++){ // For every literal in the clause...
            /* Note that the DIMACS format uses 0 as a special character, hence the variables are labelled as positive
             * (and hence non-zero) integers; Since we wish to use the variable x as an index to an array (where array
             * index starts from zero), we wish to shift the DIMACS variable label x to x - 1. Hence the DIMACS variable
             * x is encoded as 2x - 2. The negation of a variable x in DIMACS is the negative integer -x, hence applying
             * the shift and encoding, it is mapped to 2x - 1.
             */
            literals->push_back((0 < l_val[l]) ? (2 * l_val[l]) - 2 : ((-2) * l_val[l]) - 1);

            l += 1;
        }

        if(c > (t + 1) * chunk_size){
            t += 1;
        }

        clauses->at(t)->push_back(new Clause<UINT_T>(literals, t)); // Populate clauses array
    }

    // Initialise new SATInstance from specified CNF file
    auto satInstance = new SATInstance<UINT_T>(new VariablesArray<UINT_T>(v_num), n_threads);

    // Logging read complete...
    auto stop = chrono::high_resolution_clock::now();
    auto read_duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    auto timeend = chrono::system_clock::to_time_t(chrono::system_clock::now());
    time_str = string(ctime(&timeend)); time_str.pop_back();
    string log_end = "Log "; output(log_end.append(time_str) + ": Read complete; Duration: " +
                                    to_string(read_duration.count() / 1000.0) + "s\n\n", out_f, dump);

    // Display the SAT instance meta-data for monitoring purposes
    output("------------ INFORMATION ------------\n\t\t\t# Variables\t= " + to_string(satInstance->n_vars) +
    "\n\t\t\t# Clauses\t= " + to_string(satInstance->n_clauses) +
    "\n-------------------------------------\n\n", out_f, dump);

    if(dump){ // Save statistics to CSV file
        *stat_f << to_string(read_duration.count() / 1000.0) + ",";
        *stat_f << to_string(satInstance->n_vars) + ",";
        *stat_f << to_string(satInstance->n_clauses) + ",";
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

    Statistics* statistics = satInstance->solve(clauses); // Solving SAT instance and maintaining solve statistics

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

    if(satInstance->verify_validity(clauses)){ // If solution is valid
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