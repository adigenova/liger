//
// Created by Alex Digenova on 2/20/19.
//

#include "GPolisher.h"

GPolisher::GPolisher(GraphS *reduced, GraphS *lines, int ncpu) {
        this->reduced=reduced;
        this->lines=lines;
        this->cpu=ncpu;

        //we create the reduced lemon graph
        _create_lemon_graph(this->reduced,this->bp,this->g2l,this->l2g);
        uint32_t totalseq2polish=0;
        //we mark the contigs already used in edge lines, as well as the edges
    for (int i = 0; i <lines->get_number_edges(); ++i) {
        auto e = this->lines->get_edge(i);
        auto u = this->lines->get_node(e->getSource());
        auto v = this->lines->get_node(e->getTarget());
        //we use the mate-edges to mark which contigs are already in scaffolds
        if(u->getCtg_id() !=v->getCtg_id()) {
            this->ctgsinlines[u->getCtg_id()] = 1;
            this->ctgsinlines[v->getCtg_id()] = 1;
            //we get the ctg to mark the header an queue of each contig
            auto ctga=this->lines->getContigs()->getcontig(u->getCtg_id());
            auto ctgb=this->lines->getContigs()->getcontig(v->getCtg_id());
            //we mark the nodes already in scaffold sequences
            this->nodesinlines[ctga.header]=1;
            this->nodesinlines[ctga.tail]=1;
            this->nodesinlines[ctgb.header]=1;
            this->nodesinlines[ctgb.tail]=1;
            //it should be longer than the minimal contig length
            if(e->getEdge_overlap()  > 300) {
                this->edges2polish.push_back(e->getEdge());
                totalseq2polish+=e->getEdge_overlap();
            }
        }
    }
    //we print some information
    printf("A total of %d edges can be polished, the total edge sequence length is:%d\n", (int)this->edges2polish.size(),totalseq2polish);
    uint64_t iterations=0;
    for(auto eid:this->edges2polish){
        auto e=this->lines->get_edge(eid);
        //check if this work
        auto iter=this->get_alternative_path(this->reduced,bp,e,this->g2l,this->l2g);
        iterations+=iter;
    }
    printf("A total of %lld iterations were done\n",iterations);
    //we print the path that we can polish
    printf("A total of %d paths were found for polishing the %d edges\n",(int)this->path2polish.size(),(int)this->edges2polish.size());
    for(auto p:this->path2polish)
        p->print_path();


}
//method that polish each edge using the paths
void GPolisher::polish_paths(EdgeS *e, PathS *p) {
    auto e_s=e->getSource();
    auto e_t=e->getTarget();
    auto nodesinpath=p->get_nodes_in_path();
    //the source and target are ok
    if(nodesinpath[0] != e_s and nodesinpath[p->get_path_size()-1] !=e_t)
        reverse(nodesinpath.begin(),nodesinpath.end());

    //we assert that the e_source is equal to the header of the path and e_target is equal to the tail of the path
    assert(nodesinpath[0] == e_s and nodesinpath[p->get_path_size()-1] == e_t);
    //todo: we have to check if the edge has gap sequence, in principle yes
    //the edge goes from source to target
    auto e_seq=e->getEdge_gapseq();
    CJoiner aln;//this class will became my aligment engine
    float min_iden=0.8;
    if(e->getEdge_avg_cns() >=5)
        min_iden=0.9;
    if(e->getEdge_avg_cns() >=10)
        min_iden=0.95;
    bool polished=false;
    //we iter the contigs of the path to aling each of them to the edge sequence
    for(int i=0; i< p->get_path_size() - 1 ; i++){
        auto u=this->reduced->get_node(nodesinpath[i]);
        auto v=this->reduced->get_node(nodesinpath[i+1]);
        //we are in a contig
        if(u->getCtg_id() == v->getCtg_id()){
            auto contig=this->reduced->getContigs()->getcontig(u->getCtg_id());
            string strand="-";
            //we have to revcomp the contig sequence
            //means that we found a match
            if(v->getId()==contig.header) {
                strand = "+";
            }
            //we mark that the contig was used in polishing
            contig.used_in_polishing=true;
            auto c_seq= strand.compare("+") == 0 ? contig.seq : revcomp(contig.seq); //we add the contig in the correct orientation
            //computing alignment of the contig sequence with the edge sequence
            //printf("CTG=%d, ctg.h=%d ctg.t=$d strand=%s seq_ctg=%s seq_edge=%")
            printf("ALN CTG %d and EDGE %d AVG_CNS_DEPTH=%d\n",(int)c_seq.length(),(int)e_seq.length(), e->getEdge_avg_cns());
            //we have to trim aligments like in the WConsensus class, due to errors in beg and end
            auto hit=aln.compute_optimal_aligment_trim(c_seq,e_seq,false,e->getEdge_avg_cns());

            //we decide if we can polish the edge with the contig sequence
            if(hit.qcov >=75 and hit.iden >= min_iden){
                //we can polish the edge sequence
                //we mask the target sequence to find another aligment
                //e_seq.replace(hit.tstart,hit.tstop,abs(hit.tstop-hit.tstart),'N');
                //we extract the portion of the query sequence involved in polishing
                auto p_seq=c_seq.substr(hit.qstart,abs(hit.qstop-hit.qstart)+1);//the illumina seq
                //we replace the sequence of the target seq by the illumina sequence
                e_seq.replace(hit.tstart,abs(hit.tstart-hit.tstop)+1,p_seq);
                polished=true;
                //we save that the contig was used in polishing
                this->ctgused2polish[contig.ctgid]++;
            }
        }
    }
    //we have to update the polished sequence of the edge
    if(polished){
        //we save the gapseq in the strand of the edge, so we then revert when we traverse the line
        e->set_edge_gapseq(e_seq.c_str());
        e->setEdge_overlap(e_seq.length());
        e->setPolished(true);
    }
    //we return because sequence was polished
}

void GPolisher::polish() {

    //there is no path to polish
    if(this->path2edges.size() == 0)
        return;

    for(auto i=0; i<this->path2polish.size(); i++)
        this->polish_paths(this->lines->get_edge(this->path2edges[i]),this->path2polish[i]);
}




//we find the best path for the edge based on long reads and overlaps
int  GPolisher::get_alternative_path(GraphS* g, ListGraph &bp,EdgeS* e, simplehash &gToLemon,simplehash &lemonTog){
    //we obtain the nodes from the edge
    auto u=g->get_node(e->getSource());
    auto v=g->get_node(e->getTarget());
    //init path from u-> source
    auto p=new PathS(u->getId());
    //paths cointainers
    vector<PathS*> paths;//paths to look for
    vector<PathS*> vpaths;//paths founds as valid in length with respect to the reduced edge

    int iter=0;
    paths.push_back(move(p));
    //max iteration in while todo:change to a variable related to the edge size
    int max_iter_p=5000000;

    //while there are elements in the stack
    while(!paths.empty()){
        //we got a copy of the last added path
        auto pc=paths.back();
        iter++;
	//we stop the se
	if(iter > max_iter_p)
		break;
        //we remove the path from the
        paths.pop_back();
        auto lp=pc->get_last_in_path();
        //we founded a valid path
        if(lp == v->getId()){
            if(pc->get_path_size() >=4) {
                //we validate that the lengh of the path is within the expected distance
                auto plen = pc->compute_path_length(g);
                //if the lengh of the path is similar to the length of the edge so we found a candidate
                if (abs(plen - e->getBd()) <= 4 * max(e->getBd_std(), pc->getPvar())) {
                    vpaths.push_back(pc);
                    //we stop after 100 paths founds
                    if(vpaths.size() >=100)
                        break;

                } else {
                    delete pc;
                }
            }else{
                //we found the same edge
                delete pc;
            }
        }else{
            //we start our search
            auto plen=pc->compute_path_length(g);
            //the length of the path is too large with respect to the reduced edge
            if(plen > e->getBd() and (abs(plen-e->getBd()) > 4*max(e->getBd_std(),pc->getPvar()))){
                delete pc;
                continue;
            }
            //we give up due that the path lenght is too large
            if(pc->get_path_size() > 80){
                delete pc;
                continue;
            }
            //auto ng=g->get_node(lp);
            auto nl=bp.nodeFromId(gToLemon[lp]);
            //we try to extend the current path
            for (ListGraph::IncEdgeIt el(bp, nl); el != INVALID; ++el) {
                //we got the node that is oposite to the current one
                auto n2=bp.oppositeNode(nl, el);
                //we got the id of the lemon node in graph node
                auto nn=lemonTog[bp.id(n2)];
                //the node is already in path
                if(pc->is_in_path(nn)){
                    continue;
                }

                auto gn=g->get_node(nn);
                //we have to skypt nodes already in the in the lines
                /*if(gn->is_short() or gn->is_repeat()){
                    continue;
                }*/
                if(this->ctgsinlines.count(gn->getCtg_id())>0 and gn->getId() != v->getId()){
                    continue;//we cannot use nodes that are already in the Lines
                }
                //we  create the new path
                auto npath=new PathS(pc);
                npath->add_node(nn);
                //we add the path if the new path is a valid one
                if(npath->is_valid(g)){
                    //the partial path should be a long read valid, this reduce the search a lot
                    if(npath->check_long_reads(g,e) == 3) {
                        paths.push_back(npath);
                    }else{
                        delete npath;
                    }
                }else{
                    delete npath;
                }
            }
            //we delete the current partial path
            delete pc;
        }

    }


    //we found some valids paths that we can reduce if they satisfied the long read criteria
    if(!vpaths.empty()){
        int validated=0;
        auto stmp=new PathS(u->getId());
        //I have to assing a rank for each path
        for(auto pv:vpaths){
            //we compute the final length
            pv->compute_path_length(g);
            //pv->print_path();
            int score=pv->check_long_reads(g,e);
            //we save the paths to sort them and choose the best one to reduce
            if(score == 3){
                //we choose the path with more support in the long reads
                if(pv->get_path_hits() >= stmp->get_path_hits()){
                    stmp=pv;
                }
                validated++;
            }
        }

        if(validated>0){
            //stmp->print_path();
            auto bpath=new PathS(stmp);
            bpath->compute_path_length(g);
            bpath->check_long_reads(g,e);
            bpath->print_path();
            //we store the path and the id of the edge that can be polished
            this->path2polish.push_back(bpath);
            this->path2edges.push_back(e->getEdge());
        }else{
          //  cout <<"NO paths founds for edge: "<<e->getBd()<<" "<<e->getBd_std()<<" NODES(e): "<<u->getId()<<" "<<v->getId()<<endl;
        }
    }else{
        //cout <<"NO paths founds for edge: "<<e->getBd()<<" "<<e->getBd_std()<<" NODES(e): "<<u->getId()<<" "<<v->getId()<<endl;
    }

    //we release the memory allocated by the paths
    for(auto tt:paths)
        delete tt;

    for(auto tt:vpaths)
        delete tt;

    //we have to clear the paths containers
    paths.erase(paths.begin(),paths.end());
    vpaths.erase(vpaths.begin(),vpaths.end());

    return iter;
}

//create the lemon graph from the GraphS
void GPolisher::_create_lemon_graph(GraphS *g, ListGraph &bp, simplehash &gToLemon,simplehash &lemonTog){

    //we add the nodes to the lemon graph create the correct graph
    for (int j = 0; j < g->get_number_nodes() ; j++) {
        auto nd=bp.addNode();
        auto no=g->get_node(j);
        gToLemon[no->getId()] = bp.id(nd);
        lemonTog[bp.id(nd)] = no->getId();
    }
    //we add all the edges to the graph
    //ListGraph::EdgeMap<int> eid(bp);
    for (int i = 0; i <g->get_number_edges() ; ++i) {
        auto e = g->get_edge(i);
        auto u = g->get_node(e->getSource());
        auto v = g->get_node(e->getTarget());
        auto leftLemonNode = bp.nodeFromId(gToLemon[u->getId()]);
        auto rightLemonNode = bp.nodeFromId(gToLemon[v->getId()]);
        auto edge = bp.addEdge(leftLemonNode, rightLemonNode);
        //eid[edge]=e->getEdge_id();
    }
}




