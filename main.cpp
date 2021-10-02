#include <iostream>
#include <fstream>

#include "core/LaplacianLambdaCDP.h"
#include "sat_instance/SATInstance.h"

void MatrixXd_to_CSV(MatrixXd* matrix, const string& fp);

int main(int argc, char *argv[]){
    string cnf_fpath;

    // Option checking...
    if(argc <= 1){
        throw std::runtime_error("The filepath to the DIMACS formatted CNF SAT instance was not specified...exiting...");
    }
    else if(argc > 2){
        throw std::runtime_error("Too many parameters specified...exiting...");
    }
    else{
        cnf_fpath = argv[1];
    }

    // Initialise new SATInstance from specified CNF file
    auto satInstance = new SATInstance(cnf_fpath);

    // Then get the Laplacian describing the dependency graph of the SAT instance, along the the vertex sets of the
    // components of the dependency graph
    pair<MatrixXd*, vector<vector<ull>*>*> graph = satInstance->getDependencyGraph();

    // Export the Laplacian into a CSV file, for debugging and correctness checking purposes
    string csv_fpath = cnf_fpath; csv_fpath.replace(csv_fpath.size() - 3, 3, "csv");
    MatrixXd_to_CSV(graph.first, csv_fpath);

    // Get the Laplacian Lambda Core Distance Partition for the dependency graph; note that we get the 'optimal' CDP for
    // every component; for logging purposes, we print the CDP blocks for each component
    vector<cdp>* satCDPs = getLaplacianLambdaCDP(graph.first, graph.second);
    cout << "# of components = " << satCDPs->size() << endl;
    for(ull i = 0; i < satCDPs->size(); i++){
        cout << "Component " << i + 1 << ": # of partitions = " << (satCDPs->at(i))->size() << endl;
        for(ull j = 0; j < (satCDPs->at(i))->size(); j++){
            cout << "\tPartition " << j + 1 << ": ";

            for(auto c: *(*(satCDPs->at(i))).at(j)){
                cout << c << " ";
            }

            cout << endl;
        }
    }

    // Solve the SAT instance (using the Algorithmic Lovasz Local Lemma)
    VariablesArray* sat = satInstance->solve();

    // Print the variable assignment for the solution
    for(ull i = 0; i < satInstance->n_vars; i++){
        cout << "Var" << i + 1 << " = " << (sat->vars)[i] << endl;
    }

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
