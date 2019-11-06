//
// Created by Alex Digenova on 3/5/19.
//

#ifndef LIGER_SPOLISHER_H
#define LIGER_SPOLISHER_H


#include "GraphS.h"
#include "Kmer.h"
#include "CJoiner.h"

#include <algorithm>
#include <vector>
#include <cmath>

//Sequence based polisher, it align the skyped contigs to the edge consensus sequences
typedef unordered_map<uint32_t , uint32_t > simplehash2;
typedef unordered_map<kmer_type , uint32_t > simplehash3;//for 64-bits
typedef struct mmindex {
    uint64_t mh;//minimizer hash
    uint32_t seid;//edge id
    uint32_t pos;//j position
    bool strand;//save from what position the minimizer was taken
} mmindex;

//structure that store the hsps create from minimizers
typedef struct mmhsp{
    uint32_t seid;//edge id
    uint32_t cid;//contig id
    uint32_t pos;// edge position
    uint32_t cpos;//contig position
    //kmer_type m;//save the minimizer
    bool strand;
}mmhsp;

//structure that store the hits between contigs and edges
typedef struct mmhit{
    uint32_t tid;//target id ()
    uint32_t qid;//query id ()
    bool strand;//save the strand of the query with respect to the target
    uint32_t cnt=0;            // number of minimizers; if on the reverse strand
    uint32_t qs=0, qe=0, rs=0, re=0; // query start and end; reference start and end
    uint32_t mlen=0, blen=0;     // seeded exact match length; seeded alignment block length //used for stimated the identity of each hit
}mmhit;




class SPolisher {
    //simplehash2 kmers2counts;//save the kmers to count
    vector<mmindex> mpositions;
    //simplehash2 kmers2counts;
    simplehash3 minimizers2counts;
    simplehash2 ctginedges;
    int kmersize=15;
    int wsize=5;
    int maxfreq=100;//maximum minimizer freq for storing hits in the mpositions;
    GraphS *g;
    //two temporal vectors for testing only
    vector<ctg> targets;
    vector<ctg> queries;
    simplehash2 t2pos;//store the position in the array of ctgs
    simplehash2 q2pos;//store the position in the arrya of ctgs

public:
    SPolisher(GraphS *vmcg);
    SPolisher(int k, int w, int mf, GraphS *vmcg);
    SPolisher(int k, int w, int mf,string fctg, string fref);
    void build_edge_index();
    //funtion that polish the contigs
    void polish(unordered_map<int, int> & map);

    void _map_contig2edges(int cid, vector<mmhit> &hitsdb);
    void _map_contig2edges(ctg &c, vector<mmhit> &hitsdb);
    //find the largest collinear block in O(nlog(n)) time
    void find_largest_collinear_block(vector<int> &a, vector<int> &b, int offset);

    void hsps2hits(vector<mmhsp> &hsps,vector<mmhit> & localhits);

    void best_collinear2hits(vector<mmhsp> &hsps, vector<mmhit> &localhits, vector<int> &lis);

    void _hits2polish_edge(vector<mmhit> &hitsdb);

    void _polish_edge(vector<mmhit> & edgehits);

    void _get_best_ctg_for_edge(vector<mmhit> &edgehits, vector<mmhit> &bestctgs);


    bool _compute_overlap(mmhit &h, vector<mmhit> &bestctgs);

    uint32_t _compute_edge_length();

    //funtions for testing only
    void read_fasta(string filename, vector<ctg> & seqs);
    void build_index();

    void map_contigs();

    void _polish_edge2(vector<mmhit> &edgehits);

    void _hits2polish_edge2(vector<mmhit> &hitsdb);
};


#endif //LIGER_SPOLISHER_H
