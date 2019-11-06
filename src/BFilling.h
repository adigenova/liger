//
// Created by Alex Digenova on 10/2/18.
//

#ifndef LIGER_BFILLING_H
#define LIGER_BFILLING_H

#include<algorithm>

//lemon lib
#include <lemon/list_graph.h>
#include <lemon/connectivity.h>
#include <lemon/bfs.h>


//external library for consensus generation
#include "spoa/spoa.hpp"
//EDlib library
#include "edlib.h"

//std class
#include <vector>
#include <algorithm>

//own clases
#include "BSeqDB.h"
#include "GraphS.h"
#include "SLines.h"
#include "WConsensus.h"
#include "CJoiner.h"


using namespace std;
using namespace lemon;

//to declare in functions in a shorter way
typedef unordered_map<int, unordered_map<int,bool> > hashhash;
typedef unordered_map<int, vector<int> > hash2vec;
typedef unordered_map<int, int> simplehash;
typedef unordered_map<uint64_t, int> simplehash64;
typedef unordered_map<uint64_t, pair<int,int>> hash2pair64;



//class that use the reduced scaffodling graph and the long reads to fill the gaps of the backbone
class BFilling {

private:
    //GraphS* reducedG;//the reduced graph
    GraphS* linesG; //the scaffolds lines
    unordered_map<int,bool> nodesinlines; //nodes used in the lines
    vector<SLines*> linesP;//linas as paths
    vector<uint64_t> eids;//edge order
    string filename; //name of the long read file
    string prefix;//prefix of the project for ouput read information
    BSeqDB* longreadDB;//pointer to object that store the long reads in a binary file
    simplehash nodes2lines;//simple hash that store wich node to wich line
    simplehash nodes2orientations;//simple hash that store the orientation of nodes in the line
    unordered_map<uint64_t,bool> edge_dir;
    ofstream bflog;//log file

    WConsensus *wcns;//pointer to a WConsensus class
    //method that build the scaffold paths
    void build_lines();
public:
    BFilling(GraphS* lines, string fastalongreads,string prefix);
    //funtion that create the edge consensus
    void edge2consensus2(uint64_t eid,vector<string> & lseqs, vector<linking> & lin , int &lid, string &lseq);
    void fillseqs2(uint64_t eid, int n_reads,vector<string> & lseqs,  vector<linking> & lin , int &lid, string &lseq);
    //threaded version of the create edge consensus code
    void create_edge_consensus(int number_threads);
    void get_edgeseq_ends(EdgeS* e,int d,string &seq_end_source, string &seq_end_target);
    void print_cns_info(uint64_t eid);
    //unparallel version of create edges
    void create_edge_consensus();
};


#endif //LIGER_BFILLING_H
