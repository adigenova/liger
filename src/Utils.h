//
// Created by Alex Digenova on 7/2/18.
//

#ifndef LIGER_UTILS_H
#define LIGER_UTILS_H

#include <math.h>
//STL classes
#include <string>
#include <vector>
#include <regex>


using namespace std;
//join in c++
template <typename T>
string join(const string& delim, const T& v ) {
    ostringstream s;
    for (const auto& i : v) {
        if (&i != &v[0]) {
            s << delim;
        }
        s << i;
    }
    return s.str();
}

//split in c++
template<typename T>
vector<string> splits(const string& delim, const T& text){
    regex ws(delim);
    sregex_token_iterator it(text.begin(), text.end(), ws, -1);
    sregex_token_iterator reg_end;
    vector<string> results;
    for (; it != reg_end; ++it) {
        results.push_back(it->str());
    }
    return results;
}

//format fasta seq for output
//join in c++
template <typename T>
string seqformat(const int width, const T& seq ) {
    ostringstream s;
    uint32_t j=0;
    for( j=0; j<seq.length(); j+=width){
        s << seq.substr(j,width)+"\n";
    }
    if(j < seq.length()){
        s << seq.substr(j,seq.length()-j)+"\n";
    }
    //we return the formated seq
    return s.str();
}


//global array for  to compute revcomp
char const comp_tab2[] = {
        0,   1,       2,       3,       4,   5,       6,       7,       8,   9,  10,  11,      12,  13,  14,  15,
        16,  17,  18,  19,      20,  21,  22,  23,      24,  25,  26,  27,      28,  29,  30,  31,
        32,  33,  34,  35,      36,  37,  38,  39,      40,  41,  42,  43,      44,  45,  46,  47,
        48,  49,  50,  51,      52,  53,  54,  55,      56,  57,  58,  59,      60,  61,  62,  63,
        64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
        'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',  91,      92,  93,  94,  95,
        64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
        'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};

//heng li revcomp this is for printing only
template <typename T>
string revcomp(T kmer) {
    int c0, c1;
    int klen=kmer.length();
    for (int i = 0; i < klen>>1; ++i) { // reverse complement sequence
        c0 = comp_tab2[(int)kmer[i]];
        c1 = comp_tab2[(int)kmer[klen - 1 - i]];
        kmer[i] = c1;
        kmer[klen - 1 - i] = c0;
    }
    if (klen & 1) // complement the remaining base
        kmer[klen>>1] = comp_tab2[(int)kmer[klen>>1]];
    return kmer;
}

//todo: include the avg, std or other easy funtions that are used in several clases
template <typename T>
float compute_average(const T& v) {
    float avg=0;
    for (auto i : v) {
        avg+=i;
    }
    return float(avg/v.size());
}

template <typename  T>
float compute_std(const T& v, const float avg) {

    double sqtotal=0;
    for(auto i : v){
        sqtotal+=pow((avg-i),2);
    }

    float std=float(sqrt(sqtotal/(v.size()-1)));

    if(std <1){
        return 1;
    }else{
        return std;
    }
}

//functions for fast hashing of edges
template <typename  T>
uint64_t pairf(const T&  xx, const T&  yy)
{
    uint64_t p=0;
    int i=0;
    uint32_t x=xx;
    uint32_t y=yy;

    while (x||y)
    {
        p|= ((uint64_t)(x&1)<<i);     x>>=1;
        p|= ((uint64_t)(y&1)<<(i+1)); y>>=1;
        i+=2;
    }

    return p;
}

template <typename  T>
void depairf(uint64_t p, T& x, T& y)
{
    x=0;
    y=0;
    int i=0;
    while (p)
    {
        x|=((uint32_t)(p&1)<<i); p>>=1;
        y|=((uint32_t)(p&1)<<i); p>>=1;
        i++;
    }
}



#endif //LIGER_UTILS_H
