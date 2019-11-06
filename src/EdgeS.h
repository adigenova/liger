//
// Created by Alex Digenova on 7/4/18.
//

#ifndef LIGER_EDGES_H
#define LIGER_EDGES_H

#include <vector>
#include <unordered_map>
#include <algorithm>

#include "MPCollection.h"
#include "Utils.h"

#include "DNA.h"

using namespace std;

class EdgeS {
private:
   // uint32_t source;
   // uint32_t target;
    uint64_t edge;//key to edge
    int edge_id=0;
    bool deleted=0;
    bool dummy=0;
    //bundling distance
    int bd=0;
    int bd_std=0;//bundling std
    int bw=0; //bundling w
    int rw=0; //original w
    //I should move this to unit32 instead
    unordered_map<int,int> long_reads;//save the long reads used

    vector<pair<int,int>> edgelinks; //save the links used to build the edge
    //pointer to the MPCollection
    MPCollection* libs;
    void _blunding_edge();//private method to perform the edge bundling

    //todo: shold we move this code to another class?
    //variables related to the consensus gap sequence
    DnaBitset* edge_gapseq;


private:
// a pointer to edge sequence it is stored in the direction of the source->target
    bool edge_has_seq=0;//flag to know the seq
    int32_t edge_overlap=0;//set the overlap amongs contig ends
    int32_t edge_avg_cns=0;//set the consensus average coverage
    int32_t idens=0;//avg identity of source contig end aligned to the edge consensus
    int32_t covs=0;//avg coverage of source contig end aligned to the edge consensus
    int32_t ident=0;//avg identity of target contig end aligned to the edge consensus
    int32_t covt=0;//avg coverage of target contig end aligned to the edge consensus
    bool polished=0; //save if an edge has been polished
    bool edge_low_identity=0;//flag to know if the identity of contigs ends is low in any of the both ends
    bool edge_proper_gap=0;//flag that indicates if the gap was filled properly

public:
    //build the edge and
    EdgeS(uint64_t edge, MPCollection* libs,vector<pair<int,int>> elinks);

    explicit EdgeS(uint64_t edge);


    uint32_t getSource() const;

    //void setSource(uint32_t source);

    uint32_t getTarget() const;

    //void setTarget(uint32_t target);

    uint64_t getEdge() const;

    void setEdge(uint64_t edge);

    int getBd() const;

    void setBd(int bd);

    int getBd_std() const;

    void setBd_std(int bd_std);

    int getBw() const;

    void setBw(int bw);

    int getRw() const;

    void setRw(int rw);

    int getEdge_id() const;

    void setEdge_id(int edge_id);

    bool isDeleted() const;

    void setDeleted(bool deleted);

    bool is_dummy() const;

    const unordered_map<int, int> &getLong_reads() const;

    void get_longread_links(uint maxlr, vector<linking> &llinks);
    int long_reads_pass(const int lr);
    //get links larger than the specified distance
    pair<int,int> get_link_larger_than(int d);

    //gap edge methods
    bool isEdge_has_seq() const;

    void setEdge_has_seq(bool edge_has_seq);

    bool isEdge_low_identity() const;

    void setEdge_low_identity(bool edge_low_identity);

    int32_t getEdge_overlap() const;

    void setEdge_overlap(int32_t edge_overlap);

    int32_t getEdge_avg_cns() const;

    void setEdge_avg_cns(int32_t edge_avg_cns);

    bool isEdge_proper_gap() const;

    void setEdge_proper_gap(bool edge_proper_gap);

    void set_edge_gapseq(const string cns);

    int32_t getIdens() const;

    void setIdens(int32_t idens);

    int32_t getCovs() const;

    void setCovs(int32_t covs);

    int32_t getIdent() const;

    void setIdent(int32_t ident);

    int32_t getCovt() const;

    void setCovt(int32_t covt);

    string getEdge_gapseq() const;

    int get_leadseq_id(int rank);

    bool isPolished() const;

    void setPolished(bool polished);

    void get_bestn_long_reads(int n, vector<int> &long_ids);
};


#endif //LIGER_EDGES_H
