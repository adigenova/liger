//
// Created by Alex Digenova on 2/21/18.
//

#ifndef HSCAFF_DEV_MPCOLLECTION_H
#define HSCAFF_DEV_MPCOLLECTION_H


#include <vector>
#include <algorithm>
#include <thread>


class MPLibrary;
#include "MPLibrary.h"
#include "Contig.h"

using namespace std;

class MPCollection {
private:
    vector<MPLibrary*> libs;
    string filelibs;
    bool is_sort;

public:
    MPCollection(string file);
    /*void add_lib(MPLibrary* l, Contig* g);
    void sort_libs();
    void set_rank();
    MPLibrary* get_lib_by_rank(int rank);
    int size();
    vector<MPLibrary*> get_all_libs();
    void print_info();*/
   // void compute_coverage();
   // void links2file();

    void sort_libs();
    void read_libs(Contig *a,int ncpu);
    void load_lib(int lid, Contig *a);
    void print_link_libs();
    void clear_links();
    //void compute_ctg_coverage(Contig *a);
    vector<MPLibrary *> get_all_libs();
    linking get_link_by_lib_id(int libid, int linkid);
    int get_std_by_lib(int libid){return this->libs[libid]->getStd_insert_size();};
};

#endif //HSCAFF_DEV_MPCOLLECTION_H
