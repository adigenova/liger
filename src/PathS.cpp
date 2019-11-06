//
// Created by Alex Digenova on 7/9/18.
//

#include "PathS.h"


void PathS::add_node(int nid) {
    this->path.push_back(nid);
    this->inpath[nid]=true;
    this->size++;
    this->lastp=nid;
}

//default constructor
PathS::PathS(int nid) {
    this->plen=0;//path lenght
    this->clen=0;//ctg length
    this->pvar=0;//path var
    this->lastp=0;//last node added
    this->size=0;//size of the path in number of nodes
    this->path.push_back(nid);
    this->inpath[nid]=true;
    this->size++;
    this->lastp=nid;
}
//copy a path from a previous one
//default constructor
PathS::PathS(PathS* p) {
    this->plen=0;//path lenght
    this->clen=0;//ctg length
    this->pvar=0;//path var
    this->lastp=0;//last node added
    this->size=0;//size of the path in number of nodes

    //copy the path composition
    for(auto nid:p->get_nodes_in_path()){
        this->add_node(nid);
    }
}

//we ask if a node is already in path
bool PathS::is_in_path(int nid) {
    return this->inpath.count(nid) != 0;
}

//we print the path
void PathS::print_path() {
    cout << "CP : "<<this->plen<<" "<<this->pvar<<" "<<this->size<<" PW "<<this->pw<<" NODES: ";
    for(auto n:this->path){
        cout <<n<<" ";
    }
    cout <<endl;
}

//check if a path is valid or not
bool PathS::is_valid(GraphS *g) {
    //there are less than 2 edges for the moment is a valid path
    if(this->get_path_size() < 3){
        return true;
    }

    auto l=this->_check_edge(g->get_node(this->path[0])->getCtg_id(), g->get_node(this->path[1])->getCtg_id());

    for(int i=1; i<this->get_path_size() - 1; i++){
        bool n=this->_check_edge(g->get_node(this->path[i])->getCtg_id(),g->get_node(this->path[i+1])->getCtg_id());
        //means that is not a valid path, because there is no a transition between contigs and edges
        if(l == n){
            return false;
        }
        l=n;
    }
    return true;
}

bool PathS::_check_edge(int i, int i1){
    return i == i1;
}

int PathS::compute_path_length(GraphS *g) {

    //the path is short there are not sufficient edges to compute the length
    if(this->get_path_size() < 2 ){
        return 0;
    }
    //tmp var for path size
    int tclen=0;
    int tplen=0;
    float tvar=0;

    for(int i=0; i< this->get_path_size() - 1 ; i++){
        auto u=g->get_node(this->path[i]);
        auto v=g->get_node(this->path[i+1]);
        //is a contig
        if(this->_check_edge(u->getCtg_id(),v->getCtg_id())){
            //lengh of the path in contig bases
            tclen+=u->getnode_len();
            //lengh of the paths considering edges
            tplen+=u->getnode_len();
        }else{
          //is a edge
           auto e=g->get_edge(u,v);
            tplen+=e->getBd();
            tvar+=(e->getBd_std()*e->getBd_std());
        }
    }

    //we update the lenght of the paths
    this->clen=tclen;
    this->plen=tplen;
    this->pvar=int(sqrt(tvar));

    return this->plen;
}

vector<int> PathS::get_nodes_in_path() {
    return this->path;
}

int PathS::getPlen() const {
    return plen;
}

int PathS::getClen() const {
    return clen;
}

int PathS::getPvar() const {
    return pvar;
}

int PathS::getSize() const {
    return size;
}


//we check that the edge that we want reduce pass by the candidate path
int PathS::check_long_reads(GraphS *g, EdgeS* e) {
    //unordered_map<string,int > long_reads;
    if(this->get_path_size() < 2 ){
        return 3;
    }

    vector<int> long_reads;
    int nlr=0;
    int pwt=0;
    for(auto lr:e->getLong_reads()){
        //if the long read has at least 3 hits in the edge otherwise we think is an error of FAST-SG
        if(lr.second >2){
            nlr++;
            //long_reads[lr.first]=lr.second;
            long_reads.push_back(lr.first);
            //cout <<" LR in E "<<lr.first<<" "<<lr.second<<" "<<nlr<<endl;
        }
    }

    //now we check the composition of the PATH and the contigs forming the path
    int ec=0;//edge counter
    int nc=0;//contig counter
    int total_e=0;
    int total_n=0;
    //we iterate the for check the long_read composition
    for(int i=0; i< this->get_path_size() - 1 ; i++){
        auto u=g->get_node(this->path[i]);
        auto v=g->get_node(this->path[i+1]);
        //is a contig
        if(this->_check_edge(u->getCtg_id(),v->getCtg_id())){
            total_n++;
            //we have a contig and we must iterate over the contigs long reads
            int tnc=0;
            for(const auto &lr:long_reads){
                auto h=g->check_node_long_reads(u->getCtg_id(),lr);
                if(h){
                    pwt+=h;
                    //cout << lr.first << " in contig " << u->getCtg_id()<<" "<<u->getName()<<endl;
                    tnc++;
                }
            }
            if(tnc == nlr)
                nc++;
        }else{
            //is a edge
            total_e++;
            auto e2=g->get_edge(u,v);
            //we have to iterate to check if the long reads are presents
            int tec=0;
            for(const auto &lr:long_reads){
                auto h=e2->long_reads_pass(lr);
                if(h){
                    //cout << lr.first << " in edge " << e2->getEdge()<<endl;
                    pwt+=h;
                    tec++;
                }
            }
            if(tec == nlr)
                ec++;
        }
    }
    //cout << "COUNTERS " << "E "<< ec <<" N "<<nc<<" NLR "<<nlr<<" TE "<<total_e<<" TC "<<total_n<<endl;
    //the reduced edge pass by all the contigs and edges
    //we count the number of hits in the path
    this->pw=pwt;
    if(ec == total_e and nc == total_n)
        return 3;
    //the reduce edge pass by all the edges
    if(ec == total_e)
        return 2;
    //the other case is that the long reads pass by the contigs only and not the edges, which is weird but could happens
    if(nc == total_n)
        return 1;
    //there is no support in the long reads for the propossed path
    return 0;
}

void PathS::reduce_edge(GraphS *g, EdgeS *e) {
    int w=e->getBw();
    for(int i=0; i< this->get_path_size() - 1 ; i++) {
        auto u = g->get_node(this->path[i]);
        auto v = g->get_node(this->path[i + 1]);
        if(this->_check_edge(u->getCtg_id(),v->getCtg_id())){
            //is contig
        }else{
            auto e2=g->get_edge(u,v);
            int wc=e2->getBw();
            //we increase the weigth of edge path in relation to the weight of the reduced edge
            e2->setBw(wc+w);
        }
    }
    //we mark the reduced edge as deleted
    e->setDeleted(true);
}

int PathS::get_path_hits() const {
    return pw;
}



