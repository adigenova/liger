//
// Created by Alex Digenova on 7/10/18.
//

#ifndef LIGER_GREEDYS_H
#define LIGER_GREEDYS_H

#include <algorithm>
#include <lemon/list_graph.h>
#include <lemon/connectivity.h>

#include "GraphS.h"

using namespace std;
using namespace lemon;


//to declare in functions in a shorter way
typedef unordered_map<int, unordered_map<int,bool> > hashhash;
typedef unordered_map<int, vector<int> > hash2vec;
typedef unordered_map<int, int> simplehash;


class GreedyS {
    simplehash circular;
public:
    GraphS* greedy_cover(GraphS* g);
    simplehash & get_circular_nodes(){return circular;};
};


#endif //LIGER_GREEDYS_H
