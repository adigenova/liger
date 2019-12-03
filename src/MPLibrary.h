//
// Created by Alex Digenova on 2/20/18.
//

#ifndef HSCAFF_DEV_MPLIBRARY_H
#define HSCAFF_DEV_MPLIBRARY_H

#include <vector>
#include <map>
#include <math.h>
#include<algorithm>


#include "Contig.h"
//#include "DNA.h"
#include "SAMR.h"
#include "Utils.h"

using namespace std;

//struct for links amongs contigs
typedef  struct linking{
    int lid=0; // link id
    //short int lib_id=0;
    /*int ctg1=-1;
    int ctg2=-1;*/

    uint32_t c1=0;
    uint32_t c2=0;
    //link positions
    int p1=-1;
    int p2=-1;
    //inferred contig orientation
    bool cor1=0;
    bool cor2=0;
    //actual mate orientation from alignment
    bool ror1=0;
    bool ror2=0;
    //distance between contigs
    int dist=0;
    //long read information
    int lonread_id=0;
    int gaps=0;
    int gape=0;
    //encoded gap sequence
   // DnaBitset* seq;
    //edge for graph construction
    //string edge="";
    //some useful flags
    bool switched=0;
   //bool hasseq=0;

    linking(int idl, uint32_t c1, uint32_t c2,int p_1, int p_2, bool co1, bool co2, bool ro1, bool ro2, int d,bool sw, int lr, int gaps, int gape )
            :lid(idl),c1(c1),c2(c2),p1(p_1),p2(p_2), cor1(co1),cor2(co2),ror1(ro1),ror2(ro2),dist(d),switched(sw),lonread_id(lr),gaps(gaps),gape(gape) {/*this->addseq(gapseq);*/};


} linkingp;


class MPLibrary {

private:
    int avg_insert_size;
    int std_insert_size;
    int rank;//libraries are ranked according to their average insert size
    int mobs;//store the maximal number of obs to infer average insert size
    float pout; //fractions of outliers to discards while computing the average insert sizes
    string file;//file saving the SAM file
    void _infer_insert_size();
    //save the long reds that pass througth the contig
    //map<string, uint16_t > long_reads;
    vector<linking> links;
    //map to save the long reads present in the contigs
    unordered_map<int, unordered_map<int,int> > ctg2lr;

public:
    //default constructor
    MPLibrary(string file);
    //another constructor which specified the outlayers
    MPLibrary(string file, int avg_ins, int std_ins);
    //void compute_lib_average(vector<int> obs_inserts);
    vector<int> get_orientation(shortread *a, shortread *b);

    //Getter & Setter methods
    int getAvg_insert_size() const {
        return avg_insert_size;
    }

    void setAvg_insert_size(int avg_insert_size) {
        MPLibrary::avg_insert_size = avg_insert_size;
    }

    int getStd_insert_size() const {
        return std_insert_size;
    }

    void setStd_insert_size(int std_insert_size) {
        MPLibrary::std_insert_size = std_insert_size;
    }

    int getRank() const {
        return rank;
    }

    void setRank(int rank) {
        MPLibrary::rank = rank;
    }

    const string &getFile() const {
        return file;
    }

    void setFile(const string &file) {
        MPLibrary::file = file;
    }

   // void compute_lib_std(vector<int> obs_inserts);
    //function that load the library and discard the pair
    void load_library(Contig *a);

    void print_links();

    vector<linking> getLinks();
    linking get_link(int linkid){ return this->links[linkid];};
    void load_reads2ctg(Contig *a);
    //clear links
    void clear_links(){this->links.erase(this->links.begin(),this->links.end());};
    void clear_links2(){this->links.clear();this->links.shrink_to_fit();};
};


#endif //LIGER_MPLIBRARY_H
