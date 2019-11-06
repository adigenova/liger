//
// Created by Alex Digenova on 2/15/19.
//

#ifndef LIGER_CJOINER_H
#define LIGER_CJOINER_H

#include<algorithm>
#include <vector>
#include "edlib.h"

#include <string>

#include "EdgeS.h"

using namespace std;
//structure to save the alignmemt of a query and target
struct cehits{
    int  qstart=-1;
    int  qstop=-1;
    int  tstart=-1;
    int  tstop=-1;
    int  identity=0;
    int  edit_distance=0;
    int  qcov=0;
    int  tcov=0;
    bool overlap=false;
    int  overlapb=-1;
    bool strand=true;
    int w=-1;
    int qlen=0;
    int tlen=0;
    bool alive=false;
    bool partial=false;
    float iden=0;
};


class CJoiner {

public:

    cehits compute_optimal_aligment(string &q, string & t, bool log);
    //chop the target sequence after each search,stop when the alignment of the next sequence is half the score of the best aligment
    vector<cehits> compute_n_optimal_aligment(string &query, string & target, int n);
    //funtion that determine the best fit alignment to the expected distance
    void contigs_best_fit(string &source_end, string &target_end, string &consensus,string &coverage_s, string &longread,  EdgeS* e, bool revcns, bool leadrev);
    int get_obs_distance_hits(cehits *s2c, cehits *t2c);
    //print alignment
    void printAlignment(const char* query, const char* target, const unsigned char* alignment, const int alignmentLength, const int position, const EdlibAlignMode modeCode) const;
    //funtion that fill the variables of the edge
    void fill_edge_variables(int distance_obs,bool revcns,EdgeS *e, string &consensus, string &coverage_s, cehits sh, cehits th);

    cehits TrimALN(EdlibAlignResult  &result, int targetL, int queryL, int wsize, float tiden);

    cehits compute_optimal_aligment_trim(string &q, string &t, bool log, int cns_cov);
};


#endif //LIGER_CJOINER_H
