//
// Created by Alex Digenova on 7/12/18.
//

#ifndef LIGER_G2SEQ_H
#define LIGER_G2SEQ_H

#include <lemon/list_graph.h>
#include <lemon/connectivity.h>
#include <lemon/bfs.h>


#include "GraphS.h"
#include "Utils.h"


using namespace std;
using namespace lemon;

//to declare in functions in a shorter way
typedef unordered_map<int, unordered_map<int,bool> > hashhash;
typedef unordered_map<int, vector<int> > hash2vec;
typedef unordered_map<int, int> simplehash;

//class for output  the graph to sequence
class G2Seq {

public:
    void graph2agp(GraphS* g,string file);
    //in the future
    void graph2seq(GraphS* g, string prefix, bool singletons);
};


#endif //LIGER_G2SEQ_H
