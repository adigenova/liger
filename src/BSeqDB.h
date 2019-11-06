//
// Created by Alex Digenova on 8/6/18.
//

#ifndef LIGER_BSEQDB_H
#define LIGER_BSEQDB_H

//std classes
#include <unistd.h>
#include <string>
#include <iostream>
#include <fstream>
#include <random>
#include <ctime>
#include <thread>
#include <iostream>     // std::ios, std::istream, std::cout
#include <fstream>
#include <unordered_map>
#include <cstring>
#include <queue>
#include <algorithm>
#include <math.h>

#include "DNA.h"

using namespace std;
//class for store in binary files a set of long read

/* useful constants */

typedef unordered_map<uint32_t , uint32_t > simplehashu;//unsinged integers
//typedef unordered_map<int , int > simplehash;
/*#define BASE_MASK 0x3 *//* binary: 11 *//*

enum
{
    BASE_A = 0x0, *//* binary: 00 *//*
    BASE_C = 0x1, *//* binary: 01 *//*
    BASE_G = 0x2, *//* binary: 10 *//*
    BASE_T = 0x3, *//* binary: 11 *//*
};*/


struct bseq{
    uint32_t id;
    string seq;
    //save the fileid from which the sequence was read
    ushort fid;
};

class BSeqDB {

private:
    //hash that save the position of the read in the file (tellp)
    unordered_map<uint32_t,uint64_t > seq2pos;

    unordered_map<uint32_t,uint32_t > seq2chunk;//save the chunk of the file where the seq reside
    vector<pair<uint64_t ,uint64_t >> chunks;//save the chunk memory where the seq reside
    unordered_map<uint32_t, DnaBitset*> lrcache;//structure that store in memory the lonreads that should be cache

    //hash that save a refererence to the file in which the sequence was stored
    //unordered_map<uint32_t,ushort> seq2file;
    //variables for the database
    uint max_file_size=209715200 * 5; //1Gb
    //uint max_file_size=509715*10; //50Mb
    uint totalseqs=0;
    uint numberfiles=0;
    //filebuf *rb;
    filebuf index;
    //istream *is;
    string prefix;
    uint chunk_size=uint(max_file_size/4);//size of the chunk of memory 256Mb
    char* bchunk;//buffer for reading  chunk of the file instead of a single seq.
    int loaded_chunk=-1;
    uint64_t max_chunk_size=0;

public:
    BSeqDB(string fasta, simplehashu &sorder, simplehashu &lcache);
    ~BSeqDB();
    //public interfaces of the storage
    void Storeseq(uint32_t val, string seq, uint32_t len, ostream &file);

    string ReadSeq(istream &file,uint32_t& id);

    char* encodeseq(string seq);
    char* decodeseq(char* seq, int len);

    void testRead(unordered_map<uint32_t,uint32_t> & lorder);
    void DBdump();
    void DumpCache();

    bseq getSeq(uint32_t seqid);
    void load_chunk(uint32_t cid);

    bseq ChunkGetSeq(uint32_t seqid);
    //Get variables from
    uint getMax_file_size() const;

    uint getTotalseqs() const;

    uint getNumberfiles() const;

    uint getNumberLRCahed() const;
};


#endif //LIGER_BSEQDB_H
