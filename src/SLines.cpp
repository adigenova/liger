//
// Created by Alex Digenova on 9/20/18.
//

#include "SLines.h"

SLines::SLines(int node_id) : PathS(node_id) {
        //this call the corresponding constructor of PATHS
}

void SLines::compute_pos_in_line(GraphS *g) {

    //the is short there are not sufficient edges to compute the length
    if(this->get_path_size() < 2 ){
        return;
    }
    //tmp var for path size
    int tplen=0;
    //float tvar=0;
    int i=0;
    for(i=0; i< this->get_path_size() - 1 ; i++){
        auto u=g->get_node(this->path[i]);
        auto v=g->get_node(this->path[i+1]);
        node2posinline[this->path[i]]=tplen;
        //we print the coordiantes in the scf for each scaff
       // cout << this->path[i] <<" "<<tplen<<" "<<this->node2posinline[this->path[i]]<<endl;
        //is a contig
        if(this->_check_edge(u->getCtg_id(),v->getCtg_id())){
            //length of the paths considering edges
            tplen+=u->getnode_len();
        }else{
            //is a edge
            auto e=g->get_edge(u,v);
            tplen+=e->getBd();
            //tvar+=(e->getBd_std()*e->getBd_std());
        }
    }

    node2posinline[this->path[i]]=tplen;
    //cout << this->path[i] <<" "<<tplen<<" "<<this->node2posinline[this->path[i]]<<endl;
}

//we check if both nodes are in the Line
int SLines::compute_distance_in_line(int a, int b){
    if(node2posinline.count(a) and node2posinline.count(b)){
            //we return the absolute value, because a and b could be forward  or reverse to the line
            return abs(node2posinline[a]-node2posinline[b]);
    }else{
        //means that the nodes are not found in the current line
        return 0;
    }
}
//check the distance if is a circular line
int SLines::compute_distance_in_line_circular(int a, int b) {

    if(node2posinline.count(a) and node2posinline.count(b)){
        //means that a is at the end and b is at the begining
        if(node2posinline[a] > node2posinline[b]){
            return node2posinline[b]+(node2posinline[this->get_last_in_path()]-node2posinline[a])+get_circular_d();
        }else{
            return node2posinline[a]+(node2posinline[this->get_last_in_path()]-node2posinline[b])+get_circular_d();
        }
    }else{
        //means that the nodes are not found in the current line
        return 0;
    }
}
//we attempt to add a reduced edge to the fragments collection
int SLines::add_fragment_to_line(uint64_t edge, uint32_t d, uint32_t var) {
    uint32_t source,target;
    //this are the reduced edges
    depairf(edge,source,target);
    //we try in normal mode not circular
    auto dline=compute_distance_in_line(source,target);
    if(abs(int(d)-dline) <= 4*max(int(var),this->getPvar())){
        //we found a compatible fragment in the normal mode
        //we create and store the fragment
        fragmentMP f;
        f.edgeid=edge;
        f.pos= node2posinline[source];
        if(node2posinline[source] > node2posinline[target]){
            f.pos=node2posinline[target];
        }

        f.d=dline; //we save the observed distance

        f.path=true;
        //we save the fragment
        frags.push_back(f);

        return 1;
    }
    //if not we try in the circular fashion, if the Lines is a circle only. Information obtained from the scaffolder
    if(this->is_circular) {
        dline = compute_distance_in_line_circular(source, target);
        //we check is the distance is compatioble with:
        if (abs(int(d)-dline) <= 4 * max(int(var), this->getPvar())) {
            //we create two fragments one at the begining and the other at the end  of the line
            fragmentMP begin, end;
            if(node2posinline[source] < node2posinline[target]){
                    //first fragment
                    begin.edgeid=edge;
                    begin.pos=0;
                    begin.d=node2posinline[source];
                    //second fragment
                    end.edgeid=edge;
                    end.pos=node2posinline[target];
                    end.d=node2posinline[this->get_last_in_path()]-node2posinline[target];

            }else{
                //source is less than the target in the line
                begin.edgeid=edge;
                begin.pos=0;
                begin.d=node2posinline[target];
                //second fragment
                end.edgeid=edge;
                end.pos=node2posinline[source];
                end.d=node2posinline[this->get_last_in_path()]-node2posinline[source];
            }
            //we store both fragments
            begin.path=true;
            end.path=true;


            frags.push_back(begin);
            frags.push_back(end);
            return 1;
        }
    }
    //we edge was not added
    return 0;
}

//method that add the edges and contigs to the Line, edges are added only if they
void SLines::add_edges_contigs_to_line(GraphS *g, int min_frag_len) {
    //we iterate the Line

    //the is short there are not sufficient edges to compute the length
    if(this->get_path_size() < 2 ){
        return;
    }
    int i=0;
    //we create the contig fragments which extend only upto the border of the edge
    for(i=0; i< this->get_path_size() - 2 ; i++){
        auto u=g->get_node(this->path[i]);
        auto v=g->get_node(this->path[i+1]);
        //is a contig
        if(this->_check_edge(u->getCtg_id(),v->getCtg_id())){
            //we adjust the contig fragments accoding to the edges locations
            int start=node2posinline[this->path[i]],end=node2posinline[this->path[i+1]];
            //it overlap to the left with an edge
            if(i > 0 && node2posinline[this->path[i]] <= node2posinline[this->path[i-1]])
                start=node2posinline[this->path[i-1]];
            //it overlap to the rigth with an edge
            if(node2posinline[this->path[i+1]] >= node2posinline[this->path[i+2]])
                end=node2posinline[this->path[i+2]];

            fragmentMP f;
            f.edgeid=u->getCtg_id();
            f.pos= start;
            if(start > end){
                f.pos=end;
            }
            f.d=abs(start-end);

            //we save the fragment representing the contig in the collection
            f.contig=true;
            frags.push_back(f);

        }else{
           //is a edge
            auto e=g->get_edge(u,v);
            //we should look for a link that comes from a lib  that is larger than the edge
            //auto l=e->get_long_link()
            auto l=e->get_link_larger_than(min_frag_len);
            //cout << l.first<<" "<<l.second<<endl;
            //for  the moment we ask for number of links > 1 and we add the collection
            if(l.first > 0){
                fragmentMP f;
                f.edgeid=e->getEdge();
                f.pos= node2posinline[this->path[i]];

                if(node2posinline[this->path[i]] > node2posinline[this->path[i+1]]){
                    f.pos=node2posinline[this->path[i+1]];
                }
                f.d=abs(node2posinline[this->path[i]]-node2posinline[this->path[i+1]]);
		//we check that the distance is coherent with the circular one.
		if(f.d+f.pos > node2posinline[this->get_last_in_path()]){
			cout <<"CORP:"<<f.d<<" "<<f.pos<<" "<<node2posinline[this->get_last_in_path()]<<endl;
            		f.d=node2posinline[this->get_last_in_path()]-f.pos;
			cout <<"CORA:"<<f.d<<" "<<f.pos<<" "<<node2posinline[this->get_last_in_path()]<<endl;
        	}
                f.edge=true;
                frags.push_back(f);
               /* cout << "Pos edge "<<f.edgeid<<" in Line "<<f.pos<<" "<<f.d<<" "<<f.pos+f.d<<endl;*/
            }else{
                //we save the non validated edges, they should be validated by a reduced edges or are the candidate for split
                fragmentMP f;
                f.edgeid=e->getEdge();
                f.pos= node2posinline[this->path[i]];

                if(node2posinline[this->path[i]] > node2posinline[this->path[i+1]]){
                    f.pos=node2posinline[this->path[i+1]];
                }
                f.d=abs(node2posinline[this->path[i]]-node2posinline[this->path[i+1]]);

                f.edge=true;
                dfrags.push_back(f);
            }
        }
    }
    //we add the last ctg to the line
    int start=node2posinline[this->path[i]],end=node2posinline[this->path[i+1]];
    //it overlap to the left with an edge
    if(i > 0 && node2posinline[this->path[i]] <= node2posinline[this->path[i-1]])
        start=node2posinline[this->path[i-1]]+1;
    auto u=g->get_node(this->path[i]);
    //it overlap to the rigth with an edge

    //we save last contig in the array
    fragmentMP f;
    f.edgeid=u->getCtg_id();
    f.pos= start;
    if(start > end){
        f.pos=end;
    }
    f.d=abs(start-end);
	//we check that we dont exceed the size of the line this happend when we have some mate-edge overlaps
	if(f.d+f.pos > node2posinline[this->get_last_in_path()]){
		cout <<"CORP:"<<f.d<<" "<<f.pos<<" "<<node2posinline[this->get_last_in_path()]<<endl;
            	f.d=node2posinline[this->get_last_in_path()]-f.pos;
		cout <<"CORA:"<<f.d<<" "<<f.pos<<" "<<node2posinline[this->get_last_in_path()]<<endl;
        }
    f.contig=true;
    frags.push_back(f);



    /*cout << "Pos ctg "<<u->getId()<< " in Line " << start<<" "<<end<<" "<<abs(start-end)
         <<" "<<f.pos<<" "<<f.d<<" "<<f.pos+f.d<<" "<<f.edgeid<<endl;*/

}
//function that determine LQI segments within scaffolds
int SLines::LQI2break_lines(GraphS* g, int min_long_edge, vector<int>& cov ) {

    //we sort the fragments by position
    sort( frags.begin( ), frags.end( ), [ ]( const fragmentMP a, const fragmentMP b)
    {
        return a.pos < b.pos;
    });

    //we print the sorted fragments, just to check that is a good one
   /*for(auto f:frags){
       cout <<"FRAGMENTS: POS="<< f.pos<<" D="<<f.d<<" STOP="<<f.pos+f.d<<" ID="<<f.edgeid<<" C="<<f.contig<<" E="<<f.edge<<" P="<<f.path<<endl;
   }*/
    //we create the container that will hold all the queries
    clock_t begin = clock();
    //we fill the vector of counts wiht the fragments
    int maxpos=node2posinline[this->get_last_in_path()];
    cout << "Maxpos of fragment:"<<maxpos<<endl;
    assert(maxpos > 0 and maxpos <cov.size());
    //vector<int> cov(maxpos,0);//todo change for a single vector like in intervalMiss
    std::fill(cov.begin(),cov.end(),0);//takes 95 seconds
    int  avg_fragment_length=0;
    for(auto f:frags) {
        int stop=f.pos+f.d;
        avg_fragment_length+=f.d;
        for (auto it = cov.begin() + f.pos, end = cov.begin() + stop; it != end; ++it) { *it += 1; }//takes 71.5984 secs
    }

    //variables to determine LQI segments
    uint32_t start=0;
    uint32_t stop=0;
    uint32_t current=0;
    uint32_t minimum_cov=1;//this could be a variable, for the moment 0 is the target
    vector<LQI> results;
    uint32_t good_bases=0;
    //we detect the LQI segments
    for(uint32_t i=0;i<maxpos; ++i) {
        //cout << "FC: "<<i<<" "<<cov[i]<<endl;
        //we detect the LQI segments
        if(cov[i] < minimum_cov){
            if(current==0){
                current = i;
                start = i;
            }
            //allow a max dist of 10
            if(i < current+10){
                current=i;
                stop=i;
            }else{
                //is a new interval
                //cout << "LQI " << current_contig<<" "<<start <<" "<<stop<<" "<<stop-start<<endl;
                //we have to determine the edgeID to deleted
                LQI r;
                //r.ctgid=ctgid;
                r.start=start;
                r.stop=stop;
                //the mis assembly is at the interior of the contig length
                r.type=2;//means internal type
                //header type
                //we store the current LQI in the array
                results.push_back(r);
                current=0;
                start=i;
                stop=i;
            }
        }else{
            //are bases that are supportted by at lest two fragments
            good_bases++;
        }

    }

    //we clear the buffer
    std::fill(cov.begin(),cov.end(),0);
    //todo: move this stats to the ValidateBackbone class
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    avg_fragment_length=int(avg_fragment_length/frags.size());
    cout << "Total fragments " << frags.size() <<endl;
    cout << "Average fragment length " << avg_fragment_length<<endl;
    cout << "Total Examinated bases " << maxpos<<endl;
    cout << "Total Validated bases by at least "<<minimum_cov<<" fragments " << good_bases<<endl;
    cout << "Percentaje of validated bases  " << ((float)good_bases/(float)maxpos * 100)<<endl;
    cout << "Time spent in interval check " << elapsed_secs <<" secs "<<endl;
    cout << "Total LQI intervals detected "<< results.size()<<endl;
    int number_of_deleted_edges=0;
    //if there are dangereusse edges and LQI segments
    if(results.size() > 0) {
        //we use a simple array to detect the
        for(auto lqi:results){
            for (auto d:dfrags) {
                //we check if there is an overlap with the danger edges
                if((lqi.start>=d.pos  and lqi.start <=d.pos+d.d) or (lqi.stop>=d.pos  and lqi.stop <=d.pos+d.d)){
                    auto edge=g->get_edge(d.edgeid);
                   /* cout <<"POTENTIAL QUIMERIC EDGES ARRAY: EID="<<d.edgeid<<" LQI_S="<<lqi.start<<" LQI_E="<<lqi.stop<<" EID="
                         <<d.edgeid<<" E_S="<<d.pos<<" E_E="<<d.pos+d.d<<" "<<edge->getBw()<<" NLR="<<edge->getLong_reads().size()<<endl;*/

                    if(edge->getLong_reads().size() < min_long_edge) {
                        edge->setDeleted(true);
                        cout <<"DELETED EDGE ARRAY: EID="<<d.edgeid<<" LQI_S="<<lqi.start<<" LQI_E="<<lqi.stop<<" EID="
                             <<d.edgeid<<" E_S="<<d.pos<<" E_E="<<d.pos+d.d<<" "<<edge->getBw()<<" NLR="<<edge->getLong_reads().size()<<endl;
                        number_of_deleted_edges++;
                        continue;//we pass to the next resultt
                    }
                }
            }
        }

    }

    //we clean all the containers to use again this class
    frags.clear();
    dfrags.clear();
    // we erase the hash
    //node2posinline.erase(node2posinline.begin(),node2posinline.end());
    node2posinline.clear();


    return number_of_deleted_edges;
}

void SLines::save_edges_id(GraphS* g,vector<uint64_t> &edgesid) {
    //we iterate the edges only
    for(int i=1; i< this->get_path_size() - 1 ; i+=2){
        auto u=g->get_node(this->path[i]);
        auto v=g->get_node(this->path[i+1]);
        //is a contig
        if(this->_check_edge(u->getCtg_id(),v->getCtg_id())){
            continue;
        }else{
            //is a edge
            auto e=g->get_edge(u,v);
            edgesid.push_back(e->getEdge());
        }
    }
}
//function that shot the Sline fragments
void SLines::sort_fragments() {
    //we sort the fragments by position
    sort( frags.begin( ), frags.end( ), [ ]( const fragmentMP a, const fragmentMP b)
    {
        return a.pos < b.pos;
    });
    //we sort the danger edges
    sort( dfrags.begin( ), dfrags.end( ), [ ]( const fragmentMP a, const fragmentMP b)
    {
        return a.pos < b.pos;
    });
}
