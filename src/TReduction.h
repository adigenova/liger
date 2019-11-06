//
// Created by Alex Digenova on 7/5/18.
//

#ifndef LIGER_TREDUCTION_H
#define LIGER_TREDUCTION_H

#include <unordered_map>
#include <vector>
#include <memory>
#include <algorithm>
//external libraries
#include <lemon/list_graph.h>
#include <lemon/connectivity.h>

//local libraries
#include "GraphS.h"
#include "PathS.h"

using namespace std;
using namespace lemon;

//to declare in functions in a shorter way
typedef unordered_map<int, unordered_map<int,bool> > hashhash;
typedef unordered_map<int, vector<int> > hash2vec;
typedef unordered_map<int, int> simplehash;
typedef unordered_map<uint64_t, pair<int,int>> hash2pair64;
//class that perform a Transitive reduction of the graph
class TReduction {
private:
    //this store the reduced edges and the distance selected
    hash2pair64 reduced_edges;

public:
    GraphS* transitive_reduction(int lme, GraphS* g);

    int look_for_a_path(int lme,GraphS* g, ListGraph &bp, EdgeS* e, int bi, hashhash &bcc2nodes, simplehash &gToLemon,simplehash &lemonTog);
    //print the reduced edges ID, path length, path variance?, score?
    void print_reduced_edges();
    //get a pointer to the reduced edges to use them in the validation step
    hash2pair64 & get_reduced_edges(){return reduced_edges;};
};


#endif //LIGER_TREDUCTION_H
