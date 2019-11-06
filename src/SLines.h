//
// Created by Alex Digenova on 9/20/18.
//

#ifndef LIGER_SLINES_H
#define LIGER_SLINES_H


#include <algorithm>
#include <cstdlib>
#include <cmath>



#include "PathS.h"
//class for represent an scaffold sequence as a line, is an extension of the class PathS

typedef unordered_map<int, int> simplehash;
//we create a fragment to add to the lines


//Low Quality Interval LQI
typedef struct LQI{
    int64_t edgeid=-1;
    int  start=-1;
    int  stop=-1;
    int type=-1;//0 header, 1=interior,2=tail, always is type 1, beceause we are cutting edges
} LowQualityInterval;


//structure for building an intervalTREE use 80 bites per record
typedef  struct fragmentMP{
    uint64_t edgeid=0;
    int pos=0;
    int d=0;
    bool contig=false;
    bool edge=false;
    bool path=false;

} fragmentMP;



class SLines : public PathS{

private:
    //save if the line is a circular or not
    bool is_circular=false;
    int d_circular=0;//we save the length of the deleted edge
    //variables particulars of SLines class
    simplehash node2posinline;
    //vector containing the fragmentMP used to validate physically the scaffolds
    //the fragmentns comes from Reduced edges + selected edges + contigs lengths
    vector<fragmentMP> frags;//arrays of fragments
    vector<fragmentMP> dfrags;//arrays of tentative wrong edges



public:
    //this call the PathS constructor
    SLines(int node_id);
    //we declare the destructor of the class Sline
    //this method assign the coordinates of each node wihin the SLine
    void compute_pos_in_line(GraphS* g);
    //compute the distance among two  nodes
    int compute_distance_in_line(int a, int b);
    int compute_distance_in_line_circular(int a, int b);
    int add_fragment_to_line(uint64_t edge, uint32_t d, uint32_t var);
    void add_edges_contigs_to_line(GraphS* g, int min_frag_len );
    //set the line as circular or not
    void set_as_circular(int d){this->is_circular=true; this->d_circular=d;};
    int get_circular_d(){return this->d_circular;};
    int LQI2break_lines(GraphS* g, int min_long_edge, vector<int>& cov );
    void save_edges_id(GraphS* g,vector<uint64_t> &edgesid);
    //functions to implement the LQI search in another place
    const int get_line_size() {return this->node2posinline[this->get_last_in_path()];} ;
    const vector<fragmentMP> & get_fragments()  {return this->frags;};
    const vector<fragmentMP> & get_fragments_dangereus()  {return this->dfrags;};
    void sort_fragments();
};


#endif //LIGER_SLINES_H
