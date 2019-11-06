//
// Created by Alex Digenova on 9/19/18.
//

#ifndef LIGER_VALIDATEBACKBONE_H
#define LIGER_VALIDATEBACKBONE_H


//#include <cstdint>

#include <lemon/list_graph.h>
#include <lemon/connectivity.h>
#include <lemon/bfs.h>

//#include <cstdlib>

//this class validate the backbone of the assembly using the reduced edges, paths and contigs

#include "GraphS.h"
#include "Utils.h"
#include "SLines.h"

using namespace std;
using namespace lemon;

//to declare in functions in a shorter way
typedef unordered_map<int, unordered_map<int,bool> > hashhash;
typedef unordered_map<int, vector<int> > hash2vec;
typedef unordered_map<int, int> simplehash;
typedef unordered_map<uint64_t, int> simplehash64;
typedef unordered_map<uint64_t, pair<int,int>> hash2pair64;




class ValidateBackbone {
    //scaffolds as lines in the graph
    vector<SLines*> lines;
    //the nodes2lines simplehash
    simplehash nodes2lines;

public:
    //check_backbone
    GraphS* Check_Backbone(int min_frag_len,int min_long_edge,GraphS *cover,hash2pair64 &reduced_edges,simplehash &circular_nodes);
    //create the scaffolds paths
    void create_backbone_lines(GraphS *g,simplehash &circular_nodes);
    int validate_backbone(int max_scaff_size,int min_long_edge,GraphS *g);

};


#endif //LIGER_VALIDATEBACKBONE_H
