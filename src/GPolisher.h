//
// Created by Alex Digenova on 2/20/19.
//

//graph based Polisher, uses the reduced and the scaffolds lines to polish the filled gaps with repetitive and short contigs
//The steps are to build
#ifndef LIGER_GPOLISHER_H
#define LIGER_GPOLISHER_H

//external libs
//EDlib library
#include "edlib.h"
//lemon lib
#include <lemon/list_graph.h>
#include <lemon/connectivity.h>
#include <lemon/bfs.h>

#include "GraphS.h"
#include "PathS.h"
#include "CJoiner.h"

using namespace std;
using namespace lemon;

typedef unordered_map<int, int> simplehash;

class GPolisher {

private:
    GraphS *reduced;
    GraphS *lines;
    ListGraph bp;
    simplehash g2l;//graph to lemon
    simplehash l2g;//lemon to graph
    simplehash ctgsinlines;//store the contigs that are already in lines
    simplehash nodesinlines;
    vector<uint64_t> edges2polish;//edges to polish
    vector<PathS*> path2polish;//store the path that could be polished
    vector<uint64_t> path2edges;//store a relation between the paths and the edges
    simplehash ctgused2polish;//hash that store which contigs were used in the polishing step
    int cpu;//number of cores to use
    void _create_lemon_graph(GraphS *g, ListGraph &bp, simplehash &gToLemon,simplehash &lemonTog);
    int get_alternative_path(GraphS* g, ListGraph &bp,EdgeS* e, simplehash &gToLemon,simplehash &lemonTog);
public:
    GPolisher(GraphS *reduced, GraphS* lines,int ncpu);
    void polish_paths(EdgeS* e, PathS* p);
    void polish();
    //return the contigs used in polishing graph step
    simplehash & get_ctg_used_in_polishing(){return ctgused2polish;};

    //simplehash & ctg_in_scaffolds(){return ctgsinlines;};

};


#endif //LIGER_GPOLISHER_H
