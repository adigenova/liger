//
// Created by Alex Digenova on 7/3/18.
//

#ifndef LIGER_GRAPHS_H
#define LIGER_GRAPHS_H

#include <vector>

#include "Contig.h"
#include "MPCollection.h"




#include "EdgeS.h"
#include "NodeS.h"


using namespace std;


class GraphS{

private:
      //internal MAPs for nodes and edges
      vector<NodeS*> Nodes;
      vector<EdgeS*> Edges;
      //internal map for fast access to edges by key
      unordered_map<uint64_t , int> edges2pos;
      map<int, int> eid2pos;//the eid never change is like string
      unordered_map<uint32_t ,int> nodes2pos;
      //edges that goes from contig+ to contig-
      //map<string, int> dummy_edges;
      //save all the edges using the nodes ids
      //vector<pair<int,int> >edges_ids;
      //Pointers to contigs and libs clases
      Contig* contigs;
      MPCollection* mplibs;

public:
    GraphS(Contig* a, MPCollection* libs,int short_contig_len);
    //lemon is hard to use, So I'll implement my own GraphS ...
    void print_graph();
    //printing with prefix
    void print_graph(string prefix);

    int get_number_nodes(){return this->Nodes.size();};
    NodeS* get_node(int id){return this->Nodes[id];};
    NodeS* get_node(uint32_t name){return this->Nodes[this->nodes2pos[name]];};
    //vector<pair<int,int>> get_all_edges(){ return this->edges_ids;};
    //bool is_repeat(int cid){return this->contigs->is_repeat(cid);};
    //functions for get edges
    //return the edge from the array by pos
    EdgeS* get_edge(int id){return this->Edges[id];};
    //get edge by id
    EdgeS* get_edge_by_id(int id){return this->Edges[this->eid2pos[id]];};
    //return the edge by id
    //EdgeS* get_edge(string &e_name);
    EdgeS* get_edge(uint64_t e_name);
    EdgeS* get_edge(NodeS* u, NodeS* v);
    int check_node_long_reads(int id_ctg,const int lr);


    GraphS(GraphS *g);

    //Getter /setters
    Contig *getContigs() const;
    MPCollection *getMplibs() const;

    int get_number_edges();
    //get contig seq based on node_id
    string get_seq_ctg_node(uint node_id){return contigs->getcontig_seq(this->get_node(node_id)->getCtg_id());};

};


#endif //LIGER_GRAPHS_H
