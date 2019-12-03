//
// Created by Alex Digenova on 2/1/19.
//

#ifndef LIGER_WCONSENSUS_H
#define LIGER_WCONSENSUS_H


//external library for consensus generation
#include "spoa/spoa.hpp"
//EDlib library for aligments
#include "edlib.h"
//#include "IntervalTree.h"
#include "MPLibrary.h"
#include "Utils.h"


struct cwindow{
    int w;//id window
    int qstart;//start query
    int qlen;//len query
    int tstart;//start target
    int tlen;//len target
    bool  partial;//indicate if w is partial
    bool alive;//save state of window
    float iden;//idetity of window
    int qindex;//index of the query sequence
};


class WConsensus {

private :
        //matrix that store the windows to process
        vector< vector<cwindow> > windows;
        //window size
        int windowsize=200;//store the windowsize
        vector<uint32_t> coverage;//store the base coverage of the consensus sequence
        //base consensus quality is Q8=~15% error, //we start from Q33+QI
public:
    WConsensus(vector<linking> & lin, vector<string> & seqs, string &leadseq, uint32_t lid, int wsize);
    ~WConsensus();
    //create the consensus from the lead sequence in windows of size wize
    string consensusfromlead(vector<linking> & lin, vector<string> & seqs, string &leadseq, uint32_t lid);
    //create the windows
    int CheckAlignment3(EdlibAlignResult  &result, int targetL, int queryL, string  &leadseq, string &queryseq, int qindex);
    //create consensus using a backbone
    void printAlignment(const char* query, const char* target, const unsigned char* alignment,  int alignmentLength, const int correctedaln,const int position, const EdlibAlignMode modeCode) const;
    //funtion that build each consensus sequence of each window
    string generate_consensus2(vector<cwindow> &window, vector<string> &seqs, string &leadseq);
    vector<uint32_t> & get_coverage_cns() {return this->coverage;};
    string get_coverage_cns_string();
};

#endif //LIGER_WCONSENSUS_H
