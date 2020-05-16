//
// Created by Alex Digenova on 7/10/18.
//

#ifndef LIGER_MATCHINGS_H
#define LIGER_MATCHINGS_H

#include <lemon/list_graph.h>
#include <lemon/connectivity.h>
#include <lemon/matching.h>

#include "GraphS.h"

using namespace std;
using namespace lemon;

//to declare in functions in a shorter way
typedef unordered_map<int, unordered_map<int,bool> > hashhash;
typedef unordered_map<int, vector<int> > hash2vec;
typedef unordered_map<int, int> simplehash;

class MatchingS {
    //track which nodes are cirlular in case of breaking edges
    simplehash circular;
public:
    GraphS* matching_cover(GraphS* g);
    //return a hash containining the nodes being circular
    simplehash & get_circular_nodes(){return circular;};
    GraphS* matching_cover(GraphS* g, float max_repeat_factor, int short_ctg_length);


};


#endif //LIGER_MATCHINGS_H
