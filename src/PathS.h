//
// Created by Alex Digenova on 7/9/18.
//

#ifndef LIGER_PATHS_H
#define LIGER_PATHS_H

#include <vector>
#include <unordered_map>
#include <iostream>

#include "GraphS.h"

using namespace std;

//class  paths in our implementation
class PathS {


protected:
    vector<int>  path;
private:
    unordered_map<int, bool>  inpath;
    int plen=0;//path lenght
    int clen=0;//ctg length
    int pvar=0;//path var
    int lastp=0;//last node added
    int size=0;//size of the path in number of nodes
    int pw=0; //save the number of time a long read is


protected:
    //protected methods
    bool _check_edge(int i, int i1);

public:
    //explicit PathS();
    //default constructor
    explicit PathS(int nid);
    //create a copy of the actual path
    explicit PathS(PathS *p);
    //class destructor
    ~PathS(){
        path.erase(path.begin(),path.end());
        //delete(path);
        inpath.erase(inpath.begin(),inpath.end());
        //delete(inpath);
    };

    void add_node(int nid);
    bool is_in_path(int nid);
    int get_last_in_path() const {return this->lastp;};
    int get_path_size() const {return this->size;};
    void print_path();
    bool is_valid(GraphS* g);
    //function that compute the lenght of a given path
    int compute_path_length(GraphS* g);
    //return the nodes in path
    vector<int> get_nodes_in_path();

    int check_long_reads(GraphS* g, EdgeS* e);
    void reduce_edge(GraphS* g, EdgeS* e);

    //Getters
    int getPlen() const;

    int getClen() const;

    int getPvar() const;

    int getSize() const;

    int get_path_hits() const;


};

#endif //LIGER_PATHS_H
