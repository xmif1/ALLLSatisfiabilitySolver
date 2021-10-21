#include <iostream>
#include <fstream>

#include "core/LaplacianLambdaCDP.h"
#include "sat_instance/SATInstance.h"
#include <omp.h>

#define PARALLEL_RESAMPLE_LOWERBOUND 5

void MatrixXd_to_CSV(MatrixXd* matrix, const string& fp);

int main(int argc, char *argv[]){
    bool parallel = false;
    bool partition = false;

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
    else if(argc == 3 && strcmp(argv[1], "-l") == 0){
        cnf_fpath = argv[2];
        partition = true;
    }
    else if(argc == 4 && ((strcmp(argv[1], "-l") == 0 && (strcmp(argv[2], "-p")) == 0) ||
                          (strcmp(argv[1], "-p") == 0 && (strcmp(argv[2], "-l")) == 0))){

        cnf_fpath = argv[3];
        parallel = true;
        partition = true;
    }
    else if(argc == 2){
        cnf_fpath = argv[1];
    }
    else{
        throw std::runtime_error("Invalid options specified...exiting...");
    }

    // Initialise new SATInstance from specified CNF file
    auto clauses = new vector<Clause*>;
    auto satInstance = new SATInstance(cnf_fpath, clauses);

    // Then get the Laplacian describing the dependency graph of the SAT instance, along the the vertex sets of the
    // components of the dependency graph
    auto graph = SATInstance::getDependencyGraph(clauses);
    ull avg_literals_per_clause = (ull) (satInstance->n_literals / satInstance->n_clauses);
    ull parallel_resample = (avg_literals_per_clause < PARALLEL_RESAMPLE_LOWERBOUND || !parallel) ? 0 : avg_literals_per_clause;
    auto subSATInstances = satInstance->createSubSATInstances(graph.second, parallel_resample);

    // Export the Laplacian into a CSV file, for debugging and correctness checking purposes
    for(ull i = 0; i < (graph.first)->size(); i++){
        string csv_fpath = cnf_fpath; csv_fpath.replace(csv_fpath.size() - 4, 4, "_" + to_string(i) + ".csv");
        MatrixXd_to_CSV((graph.first)->at(i), csv_fpath);
    }

    // Get the Laplacian Lambda Core Distance Partition for the dependency graph; note that we get the 'optimal' CDP for
    // every component; for logging purposes, we print the CDP blocks for each component
    if(partition){
        satInstance->partition<MatrixXd*>(subSATInstances, graph.first, getLaplacianLambdaCDP, parallel);
        cout << "# of components = " << subSATInstances->size() << endl;
        for(ull i = 0; i < subSATInstances->size(); i++){
            cout << "Component " << i + 1 << ": # of partitions = " << (subSATInstances->at(i))->clausePartition->size() << endl;
            for(ull j = 0; j < (subSATInstances->at(i))->clausePartition->size(); j++){
                cout << "\tPartition " << j + 1 << ": ";

                for(auto c: *((subSATInstances->at(i))->clausePartition->at(j))){
                    cout << c << " ";
                }

                cout << endl;
            }
        }
    }

    // Solve the SAT instance (using the Algorithmic Lovasz Local Lemma)
    VariablesArray* sat = satInstance->solve(subSATInstances, parallel);

    // Print the variable assignment for the solution
//    for(ull i = 0; i < satInstance->n_vars; i++){
//        cout << "Var" << i + 1 << " = " << (sat->vars)[i] << endl;
//    }

    return 0;
}

// Utility function for exporting an Eigen MatrixXd instance to a .csv file
void MatrixXd_to_CSV(MatrixXd* matrix, const string& fp){
    const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ",", "\n");

    ofstream file(fp);
    if (file.is_open()){
        file << matrix->format(CSVFormat);
        file.close();
    }
}
