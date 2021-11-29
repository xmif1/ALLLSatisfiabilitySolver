#include <iostream>
#include <chrono>
#include <ctime>

#include "core/LaplacianLambdaCDP.h"
#include "sat_instance/SATInstance.h"
#include <omp.h>

void MatrixXd_to_CSV(MatrixXd* matrix, const string& fp);
void output(const string& str, ofstream& out_f);

int main(int argc, char *argv[]){
    bool parallel = false;
    bool partition = false;
    ull min_parallel_clauses = 0;

    string cnf_fpath;

    // Option checking...
    if(argc <= 1){
        throw std::runtime_error("The filepath to the DIMACS formatted CNF SAT instance was not specified...exiting...");
    }
    else if(argc > 5){
        throw std::runtime_error("Too many parameters specified...exiting...");
    }
    else if(argc == 3 && strcmp(argv[1], "-l") == 0){
        cnf_fpath = argv[2];
        partition = true;
    }
    else if(argc == 3 && strcmp(argv[1], "-p") == 0){
        cnf_fpath = argv[2];
        parallel = true;
    }
    else if(argc == 4 && strcmp(argv[1], "-p") == 0 && strcmp(argv[2], "-l") != 0){
        cnf_fpath = argv[3];
        parallel = true;
        min_parallel_clauses = stoll(argv[2]);
    }
    else if(argc == 4 && ((strcmp(argv[1], "-l") == 0 && (strcmp(argv[2], "-p")) == 0) ||
                          (strcmp(argv[1], "-p") == 0 && (strcmp(argv[2], "-l")) == 0))){
        cnf_fpath = argv[3];
        parallel = true;
        partition = true;
    }
    else if(argc == 5 && (strcmp(argv[1], "-l") == 0 && (strcmp(argv[2], "-p")) == 0)){
        cnf_fpath = argv[4];
        partition = true;
        parallel = true;
        min_parallel_clauses = stoll(argv[3]);
    }
    else if(argc == 5 && (strcmp(argv[1], "-p") == 0 && (strcmp(argv[3], "-l")) == 0)){
        cnf_fpath = argv[4];
        partition = true;
        parallel = true;
        min_parallel_clauses = stoll(argv[2]);
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

    // Then get the Laplacian describing the dependency graph of the SAT instance, along the the vertex sets of the
    // components of the dependency graph
    auto components = SATInstance::getDependencyGraphComponents(clauses);

    auto subSATInstances = satInstance->createSubSATInstances(components, min_parallel_clauses);
    output("# of components = " + to_string(subSATInstances->size()) + "\n\n", out_f);

    // Get the Laplacian Lambda Core Distance Partition for the dependency graph; note that we get the 'optimal' CDP for
    // every component; for logging purposes, we print the CDP blocks for each component
    if(partition){
        for(ull i = 0; i < subSATInstances->size(); i++){
            MatrixXd* laplacian = getLaplacian((subSATInstances->at(i))->clauses);
            string csv_fpath = cnf_fpath; csv_fpath.replace(csv_fpath.size() - 4, 4, "_" + to_string(i) + ".csv");
            //MatrixXd_to_CSV(laplacian, csv_fpath);

            ((subSATInstances->at(i))->partition<MatrixXd*>)(laplacian, getLaplacianLambdaCDP);
            delete laplacian;
        }

        for(ull i = 0; i < subSATInstances->size(); i++){
            string alll_str;
            if((subSATInstances->at(i))->is_ALLL_compatible()){
                alll_str = " (ALLL Compatible)";
            }

            output("Component " + to_string(i + 1) + alll_str + ": # of clauses = " +
                   to_string((subSATInstances->at(i))->clauses->size()) + ", # of partitions = " +
                   to_string((subSATInstances->at(i))->clausePartition->size()) + "\n", out_f);

            for(ull j = 0; j < (subSATInstances->at(i))->clausePartition->size(); j++){
                output("\tPartition " + to_string(j + 1) + ": ", out_f);

                for(auto c: *((subSATInstances->at(i))->clausePartition->at(j))){
                    output(to_string(c) + " ", out_f);
                }

                output("\n\n", out_f);
            }
        }
    }
    else{
        for(ull i = 0; i < subSATInstances->size(); i++){
            string alll_str;
            if((subSATInstances->at(i))->is_ALLL_compatible()){
                alll_str = " (ALLL Compatible)";
            }

            output("Component " + to_string(i + 1) + alll_str + ": # of clauses = " +
                   to_string((subSATInstances->at(i))->clauses->size()) + "\n", out_f);
        }

        output("\n", out_f);
    }

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
        output("SATISFIABLE", out_f);
//        for(ull i = 0; i < satInstance->n_vars; i++){
//            output("\nVariable " + to_string(i + 1) + " = " + to_string((sat->vars)[i]), out_f);
//        }

        out_f.close();
        return 0;
    }
    else{
        output("ERROR: Solver converged to an invalid solution!\n", out_f);

        out_f.close();
        return 1;
    }

}

// Utility function for exporting an Eigen MatrixXd instance to a .csv file
void MatrixXd_to_CSV(MatrixXd* matrix, const string& fp){
    const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ",", "\n");

    ofstream file(fp);
    if(file.is_open()){
        file << matrix->format(CSVFormat);
        file.close();
    }
}

void output(const string& str, ofstream& out_f){
    cout << str;
    out_f << str;
}