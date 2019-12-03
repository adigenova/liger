//
// Created by Alex Digenova on 7/4/18.
//

#include "EdgeS.h"



EdgeS::EdgeS(uint64_t edge, MPCollection *libs, vector<pair<int, int>> elinks) {
    this->edgelinks=move(elinks);
    //we save a reference to the libraries
    this->libs=libs;
    //depairf(edge,source,target);
    //we split the edge name to obtain the source and target
    this->dummy=0;
    //we perform the edge blunding
    this->_blunding_edge();
    this->edge=edge;
    //cout <<this->edge<<" "<<this->source<<" "<<this->target<<endl;
}

//source and target are computed in the fly from the uint64_t
EdgeS::EdgeS(uint64_t edge) {
    //will save the node code for source and target
    //depairf(edge,source,target);
    this->dummy=true;
    this->edge=edge;
}


//we perform the edge bundling
void EdgeS::_blunding_edge() {
    //we compute the edge std and avg_distance and weigth
    vector<pair<int,int>> distance;
    int i=0;
    for (auto l:this->edgelinks) {
        auto  tmp=this->libs->get_link_by_lib_id(l.first,l.second);
        //cout << this->edge <<" "<<l.first<<" "<<l.second<<" "<<tmp.longread<<" "<<tmp.has_seq()<<endl;
        pair<int,int>b(i,tmp.dist);
        distance.push_back(b);
        i++;
    }
    //we sort the list of distances for determine a simple
    sort(distance.begin(),distance.end(),[](pair<int,int> a, pair<int,int> b){ return a.second < b.second;});
    //select median using a greedy strategy
    int  m=distance[int(distance.size()/2)].first;
    auto s=this->libs->get_link_by_lib_id(edgelinks[m].first,edgelinks[m].second);
    //we obtain the boundaries of the upper and lower
    int upper=s.dist+3*this->libs->get_std_by_lib(edgelinks[m].first);
    int lower=s.dist-3*this->libs->get_std_by_lib(edgelinks[m].first);
    //cout << "selected "<<m<<" "<<s.dist<<" "<<this->libs->get_std_by_lib(edgelinks[m].first)<<" "<<lower<<" "<<upper<<endl;

    float p=0;
    float q=0;
    int w=0;

    vector<pair<int,int>> blinks;
    //we compute the std avg and w for each edge
    for (auto l:this->edgelinks) {
        auto link = this->libs->get_link_by_lib_id(l.first, l.second);
        //we have a valid distance
        if(lower <= link.dist and upper >= link.dist ){
            int stdl=this->libs->get_std_by_lib(l.first);
            //cout << lower<<" "<<upper<<" "<<this->edge<<" "<<tmp.dist<<" "<<stdl<<endl;
            p+= (float) link.dist/(stdl*stdl);
            q+= (float) 1/(stdl*stdl);
            w++;
            //save how much times this long reads is supporting the edge
            this->long_reads[link.lonread_id]++;
            blinks.push_back(l);
        }
    }

    if(w > 0) {
        //bundled w
        this->bw = w;
        //bundled distance
        this->bd = int(p / q);
        //bundled std
        this->bd_std=int(1/sqrt(q));
        //raw w
        this->rw = int(distance.size());
        //if the selected links is smaller than the total number of links means that we discarded some edges
        if(blinks.size() < this->edgelinks.size()){
            this->edgelinks=blinks;
            //we sort the links by lib id from more shorter to more larger this is useful for obtain the long reads after
            sort(edgelinks.begin(),edgelinks.end(),[](const pair<int,int> a, const  pair<int,int> b){return a.first < b.first;});
        }
    }
     //cout << "BUNDLED: "<<this->edge<<" D "<<this->bd<<" STD "<<this->bd_std <<" BW "<<this->bw<<" RW "<<this->w<<endl;
    //cout <<"RAWE: "<<this->edge<<" "<<this->bw<<" "<<this->bd<<" "<<this->bd_std<<" "<<this->rw<<endl;
}


int EdgeS::getBd() const {
    return bd;
}

void EdgeS::setBd(int bd) {
    EdgeS::bd = bd;
}

int EdgeS::getBd_std() const {
    return bd_std;
}

void EdgeS::setBd_std(int bd_std) {
    EdgeS::bd_std = bd_std;
}

int EdgeS::getBw() const {
    return bw;
}

void EdgeS::setBw(int bw) {
    EdgeS::bw = bw;
}

int EdgeS::getRw() const {
    return rw;
}

void EdgeS::setRw(int rw) {
    EdgeS::rw = rw;
}

int EdgeS::getEdge_id() const {
    return edge_id;
}

void EdgeS::setEdge_id(int edge_id) {
    EdgeS::edge_id = edge_id;
}

bool EdgeS::isDeleted() const {
    return deleted;
}

void EdgeS::setDeleted(bool deleted) {
    EdgeS::deleted = deleted;
}

bool EdgeS::is_dummy() const {
    return dummy;
}

const unordered_map<int, int> &EdgeS::getLong_reads() const {
    return long_reads;
}

void EdgeS::get_bestn_long_reads(int n,vector<int> &long_ids){
    //we erase the vector
    long_ids.clear();
    //unordered_map<int ,int> usedl;//save whit long reads as been used already
    //cout << "Number of long reads:" << long_reads.size()<<endl;
    //we sort the long reads and we take the first 15
    vector<pair<int,int>> vlonreads;
    for(auto l:long_reads){
        pair<int,int> lc(l.first,l.second);
        vlonreads.push_back(lc);
    }

    //we sort by number of pairs suporting the long edge
    sort(vlonreads.begin(),vlonreads.end(),[](const pair<int,int> a,const pair<int,int> b){return a.second > b.second;});
    for(auto i=0; i < vlonreads.size(); ++i){
        //cout <<i<<" "<<vlonreads[i].first<<" "<<vlonreads[i].second<<endl;
        //usedl[vlonreads[i].first]=0;
        long_ids.push_back(vlonreads[i].first);
        if(i >= n){
            break;
        }
    }

}


//determine if the long reads is in the set of long_reads that traverse the edge
int EdgeS::long_reads_pass(const int lr) {
    //this return 0 or 1 in case when the long reads is in the map
    //return this->long_reads.count(lr) > 0;
    if(this->long_reads.count(lr)){
      return this->long_reads[lr];
    }else{
        return 0;
    }
}


//determine if the long reads is in the set of long_reads that traverse the edge
int EdgeS::get_leadseq_id(int rank) {
    int max=0;
    int id=-1;
    //we find the sequence having the most overlaps
    for (auto lr:long_reads){
            if(lr.second > max)
                id=lr.first;
    }

    return id;
}


uint32_t EdgeS::getSource() const {
    uint32_t source,target;
    depairf(edge,source,target);
    return source;
}

/*void EdgeS::setSource(uint32_t source) {
    EdgeS::source = source;
}*/

uint32_t EdgeS::getTarget() const {
    uint32_t source,target;
    depairf(edge,source,target);
    return target;
}

/*void EdgeS::setTarget(uint32_t target) {
    EdgeS::target = target;
}*/

uint64_t EdgeS::getEdge() const {
    return edge;
}

void EdgeS::setEdge(uint64_t edge) {
    EdgeS::edge = edge;
}

pair<int,int> EdgeS::get_link_larger_than(int d) {
    int counter=0;
    int lread=0;
    unordered_map<uint32_t, bool> hcount;
    for (auto l:this->edgelinks) {
        auto link = this->libs->get_link_by_lib_id(l.first, l.second);
        if(link.gape-link.gaps+400 >= d){
            //we print the link
            //cout << "LINK "<<l.first<<" "<<link.c1<<" "<<link.p1<<" "<<link.c2<<" "<<link.p2<<" "<<link.dist<<" "<<link.lonread_id<<" "<<link.gape-link.gaps+400<<endl;
            counter++;
            if(hcount.count(link.lonread_id) == 0){
                lread++;
                hcount[link.lonread_id]=true;
            }

        }

    }
    pair<int,int> res(counter,lread);
    //cout << "LINKS larger than "<<d<<" "<<res.first<<" long reads "<<res.second<<endl;
    //cout << "WT "<<this->getBw()<<endl;

    return res;
}




//funtion that select the best long reads for construction the consensus sequence
void EdgeS::get_longread_links(uint maxlr, vector<linking> &llinks) {
    //vector<linking> llinks;
    unordered_map<int ,int> usedl;//save whit long reads as been used already
    //cout << "Number of long reads:" << long_reads.size()<<endl;
    //we sort the long reads and we take the first 15
    vector<pair<int,int>> vlonreads;
    for(auto l:long_reads){
        pair<int,int> lc(l.first,l.second);
        vlonreads.push_back(lc);
    }

    //we sort by number of pairs suporting the long edge
    sort(vlonreads.begin(),vlonreads.end(),[](const pair<int,int> a,const pair<int,int> b){return a.second > b.second;});

    for(auto i=0; i < vlonreads.size(); ++i){
           //cout <<i<<" "<<vlonreads[i].first<<" "<<vlonreads[i].second<<endl;
           usedl[vlonreads[i].first]=0;
           if(i >= maxlr){
               break;
           }
    }
    //we select a template:
    //int best_backbone=vlonreads[0].first;
    // todo: should I select a template for building the consesnsus?,
    // the template is the longer sequence that span the edge
    int counter=0;
    //auto best=this->libs->get_link_by_lib_id(edgelinks[0].first, edgelinks[0].second);
    for (auto l:this->edgelinks) {
        auto link = this->libs->get_link_by_lib_id(l.first, l.second);
        //maybe I should remove this , but essentially avoid to use short fragments to build the consensus sequence, this increase
        //a little bit the time but we can build better consensus sequences for repetitives regions, essentially.
        /*if(l.first == 0){
            continue;
        }*/
        //cout << "links order in edges "<< l.first<<" "<<l.second<<" "<<link.lonread_id<<" "<<long_reads_pass(link.lonread_id)<<endl;
        if(usedl.count(link.lonread_id)){
            if(usedl[link.lonread_id]< 1){
            usedl[link.lonread_id]++;
            llinks.push_back(link);
              //  cout << "***links order in edges "<< l.first<<" "<<l.second<<" "<<link.lonread_id<<" "<<long_reads_pass(link.lonread_id)<<endl;
                counter++;
            }
        }

        if(counter >= maxlr){
            break;
        }
    }

    //we clean the local hash before leave the funtion
    usedl.erase(usedl.begin(),usedl.end());
    vlonreads.erase(vlonreads.begin(),vlonreads.end());
}





bool EdgeS::isEdge_has_seq() const {
    return edge_has_seq;
}

void EdgeS::setEdge_has_seq(bool edge_has_seq) {
    EdgeS::edge_has_seq = edge_has_seq;
}

bool EdgeS::isEdge_low_identity() const {
    return edge_low_identity;
}

void EdgeS::setEdge_low_identity(bool edge_low_identity) {
    EdgeS::edge_low_identity = edge_low_identity;
}

int32_t EdgeS::getEdge_overlap() const {
    return edge_overlap;
}

void EdgeS::setEdge_overlap(int32_t edge_overlap) {
    EdgeS::edge_overlap = edge_overlap;
}

int32_t EdgeS::getEdge_avg_cns() const {
    return edge_avg_cns;
}

void EdgeS::setEdge_avg_cns(int32_t edge_avg_cns) {
    EdgeS::edge_avg_cns = edge_avg_cns;
}

bool EdgeS::isEdge_proper_gap() const {
    return edge_proper_gap;
}

void EdgeS::setEdge_proper_gap(bool edge_proper_gap) {
    EdgeS::edge_proper_gap = edge_proper_gap;
}

void EdgeS::set_edge_gapseq(const string  cns) {
    if(cns.length() > 0){
        //we save the sequence only if it larger than 0
        this->edge_gapseq=new DnaBitset(cns.c_str(),cns.length());
        this->setEdge_has_seq(true);
    }
}

int32_t EdgeS::getIdens() const {
    return idens;
}

void EdgeS::setIdens(int32_t idens) {
    EdgeS::idens = idens;
}

int32_t EdgeS::getCovs() const {
    return covs;
}

void EdgeS::setCovs(int32_t covs) {
    EdgeS::covs = covs;
}

int32_t EdgeS::getIdent() const {
    return ident;
}

void EdgeS::setIdent(int32_t ident) {
    EdgeS::ident = ident;
}

int32_t EdgeS::getCovt() const {
    return covt;
}

void EdgeS::setCovt(int32_t covt) {
    EdgeS::covt = covt;
}

string EdgeS::getEdge_gapseq() const {
    if(this->isEdge_has_seq())
    return this->edge_gapseq->to_string();
    else
        return "";
}

bool EdgeS::isPolished() const {
    return polished;
}

void EdgeS::setPolished(bool polished) {
    EdgeS::polished = polished;
}
