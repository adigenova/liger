//
// Created by Alex Digenova on 3/6/18.
//

#ifndef LIGER_CONTIG_H
#define LIGER_CONTIG_H

#include <vector>
#include <map>
#include <unordered_map>
#include <iterator>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>


#include "Utils.h"


using namespace std;

//struct for contigs
typedef  struct ctg{
    int ctgid=-1;
    string namectg="";
    int length=-1;
    string seq="";//contigs are
    float coverage=-1;
    int num_reads=0;
    bool repeat=0;//save if the contig is repeat
    unordered_map<int, int> long_reads;
    uint32_t header=0;//header of contig(+)
    uint32_t tail=0;//tail of contig(-)
    bool used_in_polishing=false;
} ctgs;

class Contig {
private:
    string fasta="";
    string fai="";
    vector<ctg> contigs;
    map<string,int> string2id;
    float avg_ctg_cov;
    float std_ctg_cov;

public:
    Contig(string filename);
    int ctg2id(string name);
    void print_contigs();
    void print_contigs_file(string file);
    ctg getcontig(string name);
    ctg getcontig(int id);
    vector<ctg> get_all_contigs();
    int getcontig_length(string name);
    int getcontig_length(int id);
    string getcontig_seq(string name);
    string getcontig_seq(int id);
    int get_number_contig(){return this->contigs.size();};
    void set_ctg_coverage(int id, float c){this->contigs[id].coverage=c;};
    void add_ctg_reads(int id, int c){this->contigs[id].num_reads+=c;};

    //funtions to add long reads and check the long reads stored in the map
    void add_long_reads(int id, const int lr){this->contigs[id].long_reads[lr]++;};
    void add_long_reads(int id, const int lr, const int counter){this->contigs[id].long_reads[lr]+=counter;};
    int check_long_read(int id, const int lr);//{return this->contigs[id].long_reads.count(lr) > 0;};

    int get_contig_id(string name);
    uint32_t get_contig_header(string name);
    uint32_t get_contig_tail(string name);

    void mark_repeats();
    //the coverage comes from illumina reads
    void mark_repeats(string filename);
    //getters and setters
    float getAvg_ctg_cov() const;

    void setAvg_ctg_cov(float avg_ctg_cov);

    float getStd_ctg_cov() const;

    void setStd_ctg_cov(float std_ctg_cov);

    string get_contig_name(int id){return this->contigs[id].namectg;};

    bool is_repeat(int id){return this->contigs[id].repeat;};
    void set_repeat(int id){ this->contigs[id].repeat=true;};
    void set_used_in_polishing(int id){this->contigs[id].used_in_polishing=true;};

    void mark_repeats_astat(string filename,float rf);
};



#endif //LIGER_CONTIG_H
