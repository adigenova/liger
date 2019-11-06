//
// Created by Alex Digenova on 7/3/18.
//

#include "GraphS.h"


GraphS::GraphS(Contig* a, MPCollection* libs, int short_contig_len) {
    //Step 1: for each contig we create two nodes a + and -
    int plus=0;
    int ecounter=0;//edge counter
    for(int i=0; i<a->get_number_contig();i++){
        auto ctg=a->getcontig(i);
        //int id, string n, int ctg
        bool cshort=0;

        //todo: this shold be a parameter of the method
        if(ctg.length < short_contig_len){
            cshort=1;
        }

        auto u=new NodeS(ctg.header,ctg.ctgid,ctg.repeat,cshort,ctg.length);
        this->nodes2pos[u->getId()]=plus;
        plus++;
        auto v=new NodeS(ctg.tail,ctg.ctgid,ctg.repeat,cshort,ctg.length);

        this->Nodes.push_back(u);
        this->Nodes.push_back(v);
        this->nodes2pos[v->getId()]=plus;
        plus++;
        auto edge=new EdgeS(pairf(u->getId(),v->getId()));
        edge->setEdge_id(ecounter);
        this->Edges.push_back(edge);
        this->edges2pos[edge->getEdge()]=ecounter;
        //todo remove eid2pos because we have unit64t to post
        this->eid2pos[edge->getEdge_id()]=ecounter;
        ecounter++;
    }

    //Step 2: we iterate all the links to create the initial edges
    unordered_map<uint64_t , vector< pair<int,int> > > edges2links;
    for(auto lib:libs->get_all_libs()){
        //cout << lib->getRank() <<" "<<lib->getAvg_insert_size() <<endl;
        for(auto link:lib->getLinks()){
            pair<int,int> tmp(lib->getRank(),link.lid);
            //we build the unique id
            edges2links[pairf(link.c1,link.c2)].push_back(tmp);
        }
    }

    //Step 3: we create the edges (bundling of the links)
    for(const auto &e : edges2links){
        //it do the bundling of the edge
        auto edge=new EdgeS(e.first,libs,e.second);
        //we ask if the weigth is larger than 5
        if(edge->getBw() >= 5){
            //we save our map
            edge->setEdge_id(ecounter);
            this->Edges.push_back(edge);
            this->edges2pos[edge->getEdge()]=ecounter;
            this->eid2pos[edge->getEdge_id()]=ecounter;
            ecounter++;
        }else{
            //todo: can we recover this low quality edges? for the moment we discard them
        }
    }
    //we save the reference to contigs and libs objects
    this->contigs=a;
    this->mplibs=libs;
    //we clear the edges2links matrix
    edges2links.erase(edges2links.begin(),edges2links.end());
}


void GraphS::print_graph() {
    //we print edges including the dummy ones
    for (auto e:this->Edges){
        cout << "GRAPH: "<< e->getSource() <<" "<<e->getTarget()<<" "<<e->getBw()<<" "<<e->getBd()<<" "<<e->getEdge()<<endl;
    }
}

void GraphS::print_graph(string prefix) {
    //we print edges including the dummy ones
    for (auto e:this->Edges){
        cout <<prefix<< " : "<< e->getSource() <<" "<<e->getTarget()<<" "<<e->getBw()<<" "<<e->getBd()<<" "<<e->getEdge_id()<<endl;
    }
}
//we obtain the edge by two nodes
EdgeS* GraphS::get_edge(NodeS *u, NodeS *v) {
    uint64_t e=pairf(u->getId(),v->getId());
    if(v->getId() < u->getId()){
          e=pairf(v->getId(),u->getId());
    }

    return this->get_edge(e);
}


EdgeS *GraphS::get_edge(uint64_t e_name) {
    if(this->edges2pos.count(e_name)>0){
        return this->Edges[this->edges2pos[e_name]];
    }else{
        cout << "Edge : "<<e_name<<" Not in graph"<<endl;
        exit(1);
    }
}


//we ask if the long read pass thought the contig
int GraphS::check_node_long_reads(int id_ctg,const int lr) {
    //return this->contigs->check_long_read(id_ctg,lr) > 0;
    return this->contigs->check_long_read(id_ctg,lr);
}

//constructor that copy a graphS object discarding deleted edges
GraphS::GraphS(GraphS *g) {
    //we copy the nodes of the graph and the node id is immutable
    //we copy because the adj list change
    //todo: check if is necesary a copy, because we are not using the adj_list for the moment
    for (int i = 0; i <g->get_number_nodes() ; ++i) {
        auto n=g->get_node(i);
        //auto n=new NodeS(g->get_node(i));
        this->Nodes.push_back(n);
        this->nodes2pos[n->getId()]=i;
    }
    //we save references for the non deleted edges of the graph
    int ecounter=0;
    for (int j = 0; j <g->get_number_edges() ; ++j) {
        auto e=g->get_edge(j);
        //the edge was reduced or deleted
        if(e->isDeleted()){
            continue;
        }
        //we store the edge in the new array
        this->Edges.push_back(e);
        this->edges2pos[e->getEdge()]=ecounter;
        this->eid2pos[e->getEdge_id()]=ecounter;
        ecounter++;
    }
    //we copy the graph g excluding deleted edges
    this->contigs=g->getContigs();
    this->mplibs=g->getMplibs();
}

Contig* GraphS::getContigs() const {
    return contigs;
}

MPCollection* GraphS::getMplibs() const {
    return mplibs;
}

int GraphS::get_number_edges() {
    return this->Edges.size();
}


