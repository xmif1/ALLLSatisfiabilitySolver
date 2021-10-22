//
// Created by Xandru Mifsud on 09/09/2021.
// Please cite as: Mifsud X., "\[Lambda]--Core distance partitions", Linear Algebra Appl. (2021), https://doi.org/10.1016/j.laa.2020.12.012
//

#ifndef LAPLACIANLAMBDACDP_H
#define LAPLACIANLAMBDACDP_H

// Macros for determining if two values are approximately equal to each other for a given epsilon
#define EPSILON 1.0e-15
#define approx_eig(l1, l2, N) ((l1 < (l2 + (EPSILON*N))) && (l1 > (l2 - (EPSILON*N))))
#define approx_zero(l, N) ((l < (EPSILON*N)) && (l > (-1*EPSILON*N)))

#include <iostream>
#include <vector>
#include <set>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "../sat_instance/SATInstance.h"

using namespace Eigen;
using namespace std;

typedef unsigned long long ull;
typedef vector<set<ull>*>* cdp; // CDPs are given as a pointer to a vectors of pointers to sets, each set being a block

MatrixXd* getLaplacian(vector<Clause*>* clauses);
cdp getLaplacianLambdaCDP(MatrixXd* laplacian);
cdp nestedNeighbourhoods(MatrixXd* laplacian, set<ull>* cores);
double LambdaCDPIndex(cdp lambdaCDP, ull N);

MatrixXd* getLaplacian(vector<Clause*>* clauses){
    ull n_clauses = clauses->size();
//    ull sum_of_degrees = 0;
//    ull avg_clause_literals = 0;
//    for(ull i = 0; i < n_clauses; i++){
//        sum_of_degrees += (clauses->at(i))->degree;
//        avg_clause_literals += (clauses->at(i))->n_literals;
//    }
//
//    avg_clause_literals = (ull) avg_clause_literals / n_clauses;

    // Initialise a MatrixXd instance and initialise all the entries to zero.
    auto laplacian = new MatrixXd(n_clauses, n_clauses);

    for(ull i = 0; i < n_clauses; i++){
        for(ull j = i + 1; j < n_clauses; j++){
            // If the clauses i and j are dependent, then:
            if(SATInstance::dependent_clauses(clauses->at(i), clauses->at(j))) {
                (*laplacian)(i, j) = -1; // The entry (i, j) of the Laplacian is -1
                (*laplacian)(j, i) = -1; // The entry (j, i) of the Laplacian is -1 by symmetry
            } // Otherwise if independent, the entries (i, j) and (j, i) of the Laplacian are 0
        }

        (*laplacian)(i, i) = (clauses->at(i))->degree;
    }

    return laplacian;
}

/* Computes the 'optimal' Lambda-CDP (LCDP) for the Laplacian of a given graph. The optimal LCDP is defined as the LCDP
 * with the largest LCDP Index; this is generally the LCDP with the largest number of blocks, of approximately equal size
 * (or the closest possible from all generated LCDPs).
 *
 * Parameters:
 *    MatrixXd* laplacian : A pointer to a MatrixXd instance representing an Laplacian matrix of a graph. It is assumed
 *                          that laplacian is a valid N x N symmetric matrix satisfying Laplace's equation L = D - A.
 *
 * Outline:
 *  (1) Compute the eigen-decomposition of L, given by L = USU^T, there the columns of U are the eigenvectors of L and
 *      the diagonal entries of S are the corresponding eigenvalues.
 *  (2) Compute the set of core vertices for each eigenvalue of L via the eigen-decomposition.
 *  (3) Compute the LCDP associated with each eigenvalue via the corresponding set of core vertices.
 *  (4) Return the LCDP with the largest LCDP Index, while freeing any memory associated with the un-returned LCDPs.
 */
cdp getLaplacianLambdaCDP(MatrixXd* laplacian){
    ull N = laplacian->rows();

    /* Theorem: For a real-symmetric positive semi-definite matrix, its singular value decomposition (SVD) and eigen-
     *          decomposition correspond to each other.
     *
     * Proposition: The Laplacian of a graph is a real-symmetric positive semi-definite matrix.
     *
     * Consequently, we can use the divide-and-conquer SVD implementation of the Eigen package, which albeit slower then
     * the eigen-decomposition implementation for small dimensionality, the SVD implementation scales better for larger
     * dimensionality while retaining superior accuracy.
     */
    auto svd = laplacian->bdcSvd(ComputeThinU); // compute the SVD for the Laplacian L
    const auto& eigs = svd.singularValues(); // get the eigenvalues of L; note that these are sorted in decreasing order

    // Initialise a vector of set<ull> instances which will correspond to the possible cores for all eigenvalues of L.
    // Note that 0 is always an eigenvalue of the Laplacian, and hence we are guaranteed to have at least one core vertex
    // set (and hence LCDP).
    auto core_sets = new vector<set<ull>*>;

    // Compute the set of core vertices for each eigenvalues of L
    for(ull eig_min = 0; eig_min < N; eig_min++){
        // Find the maximum range of indices (eig_min to eig_max) of the eigs vector such that all entries are (approx)
        // the same eigenvalue. We use an approximate, up to some epsilon, to account for numerical errors.
        ull eig_max = eig_min + 1;
        for(; eig_max < N; eig_max++){
            if(!approx_eig(eigs(eig_max), eigs(eig_min), N)){
                break;
            }
        }
        eig_max--;

        /* Hence the columns of U in the SVD of L with indices between eig_min and eig_max for an eigen-basis for the
         * same eigenvalue. From these columns we can then construct the core vertex set of this eigenvalue for every
         * component of G; a vertex v of a component is in the core vertex set of that component if there is one such
         * column vector with the v^th entry (approx.) not zero. We use an approximate, up to some epsilon, to account
         * for numerical errors.
         */
        auto cores = new set<ull>;

        for(ull v = 0; v < N; v++){ // for each vertex v
            for(ull e = eig_min; e <= eig_max; e++){ // for each column vector of U for this eigenvalue
                // if the v^th entry is non-zero, v is a core vertex
                if(!approx_zero((svd.matrixU())(v, e), N)){
                    cores->insert(v); // hence we insert it into the core vertex set
                    break; // and there is no need to check the v^th entry of the remaining columns
                }
            }
        }

        // If the graph has at least 1 core vertex for the current eigenvalue, push the set of core vertices to the
        // collection of core vertex sets.
        if(!cores->empty()){
            core_sets->push_back(cores);
        }

        eig_min = eig_max;
    }

    /* For each set of core vertices associated with a distinct eigenvalue, we find the full LCDP and its LCDP Index.
     * Note that we maintain the largest (optimal) LCDP Index and the associated LCDP; if the LCDP associated with a
     * core vertex set has and LCDP Index less than the current maximum, then it is not 'optimal' and hence we discard
     * it (in particular we carry out memory management on the go, especially useful for large numbers of vertices).
     */
    cdp ret_cdp, curr_cdp;
    double max_index = 0;
    for(auto& cv: *core_sets){
        curr_cdp = nestedNeighbourhoods(laplacian, cv); // calculate the full LCDP
        double index = LambdaCDPIndex(curr_cdp, N); // and find the corresponding index
        if(index >= max_index){ // if index exceed the current maximum, this LCDP is more optimal
            max_index = index;
            ret_cdp = curr_cdp;
        }
        else{ // otherwise we discard it by freeing any associated memory
            for(auto p: *curr_cdp){ delete p;}
            delete curr_cdp;
        }
    }

    delete core_sets;

    return ret_cdp; // lastly, we return the optimally computed LCDP
}

/* Given a Laplacian matrix of a graph and a subset of the vertex set, this function computes the nested neighbourhoods
 * of this subset, which is a partition of the vertex set when the graph is connected. For a subset S of the vertex set,
 * we define its neighbourhood N(S) as the set of vertices not in S that are adjacent to some vertex in S.
 */
cdp nestedNeighbourhoods(MatrixXd* laplacian, set<ull>* cores){
    // Vector structure to hold the nested neighbourhood set instances, in order of distance.
    auto nestedN = new vector<set<ull>*>;
    nestedN->push_back(cores);

    // The set 'remaining_vertices' initially has all the vertices except those in the initial subset passed. Termination
    // follows when this set is empty. Termination is guaranteed if the graph is connected.
    set<ull> remaining_vertices;
    for(ull v = 0; v < laplacian->rows(); v++){ // for each vertex in the graph
        if(cores->count(v) == 0){ // insert it into the set 'remaining_vertices' if it is not in the set 'cores'
            remaining_vertices.insert(v);
        }
    }

    // Construct nested neighbourhoods until we have visited all the vertices in the graph...
    while(!remaining_vertices.empty()){
        auto curr_neighbourhood = new set<ull>; // new set instance for the next nested neighbourhood

        for(auto v = remaining_vertices.begin(); v != remaining_vertices.end();){ // Iterate through the unvisited vertices
            bool v_deleted = false;
            for(auto& u: *(nestedN->back())){ // for every vertex in the last neighbourhood...
                if((*laplacian)(*v, u) == -1){ // if adjacent to some vertex in the last neighbourhood
                    curr_neighbourhood->insert(*v); // then insert into new neighbourhood
                    v = remaining_vertices.erase(v); // and remove from the list of remaining vertices
                    v_deleted = true;

                    break;
                }
            }

            // If not deleted, i.e. not a neighbour to any vertex in the last neighbourhood, increment iterator.
            if(!v_deleted){ ++v;}
        }

        nestedN->push_back(curr_neighbourhood);
    }

    return nestedN;
}

/* Computes the LCDP Index for a given LCDP for a graph on N vertices. The formula is given by:
 *                                (-1/N)*Sum(|p| * Log10(|p| / N), p in LCDP)
 */
double LambdaCDPIndex(cdp lambdaCDP, ull N){
    double sum = 0;

    for(auto p: *lambdaCDP){ // Sum over p in LCDP
        sum += ((p->size()) * log10((double) (p->size())/N)); // |p| * Log10(|p| / N)
    }

    return (-1.0) * (sum / N);
}

#endif //LAPLACIANLAMBDACDP_H
