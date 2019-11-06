//
// Created by Alex Digenova on 7/5/18.
//

#include "TReduction.h"

GraphS* TReduction::transitive_reduction(int lme,GraphS *g) {
    //we create the Lemon graph to perform the sets of operations that we want to compute
    ListGraph bp;
    simplehash gToLemon;
    simplehash lemonTog;
    //we add the nodes to the lemon graph create the correct graph
    for (int j = 0; j < g->get_number_nodes() ; j++) {
        auto nd=bp.addNode();
        auto no=g->get_node(j);
        gToLemon[no->getId()] = bp.id(nd);
        lemonTog[bp.id(nd)] = no->getId();
    }
    //we add all the edges to the graph
    ListGraph::EdgeMap<int> eid(bp);
    for (int i = 0; i <g->get_number_edges() ; ++i) {
        auto e = g->get_edge(i);
        auto u = g->get_node(e->getSource());
        auto v = g->get_node(e->getTarget());
        auto leftLemonNode = bp.nodeFromId(gToLemon[u->getId()]);
        auto rightLemonNode = bp.nodeFromId(gToLemon[v->getId()]);
        auto edge = bp.addEdge(leftLemonNode, rightLemonNode);
        eid[edge]=e->getEdge_id();
    }
   //we compute the biconnected components of the graph
    ListGraph::EdgeMap<int> bcc(bp);
    int nbcc=biNodeConnectedComponents(bp,bcc);
    cout << "Number of BCC: " <<nbcc<<endl;
    //hash of bcc 2 nodes
    hashhash bcc2nodes;
    //hash of bcc 2 edges
    hash2vec bcc2reduce;
    //report some variables related to the reduced edges
    uint32_t candidate_reduce_e=0;
    uint32_t candidate_reduce_w=0;
    uint64_t candidate_reduce_l=0;
    //we check what's edges can be reduced
    for (ListGraph::EdgeIt e(bp); e != INVALID; ++e){
        auto u=g->get_node(lemonTog[bp.id(bp.u(e))]);
        auto v=g->get_node(lemonTog[bp.id(bp.v(e))]);
        //we mark the nodes in the bcc
        bcc2nodes[bcc[e]][u->getId()]=1;
        bcc2nodes[bcc[e]][v->getId()]=1;
        //we dont reduce repetitives edges
        if(u->is_repeat() or v->is_repeat()){
            continue;
        }
        //we dont reduce the edge if one of the two contigs is short
        if(u->is_short() or v->is_short()){
            continue;
        }
        //we dont reduce short edges
        //cout << " LG: "<< bp.id(bp.u(e)) << "-" << bp.id(bp.v(e)) <<" "<<u->getName()<<"--"<<v->getName()<<" "<< bcc[e] <<endl;
            //we got the edge from the graph
        auto es=g->get_edge_by_id(eid[e]);
        //we dont reduce dummy edges
        if(es->is_dummy()){
            continue;
        }
        //we dont reduce short edges
        if(es->getBd() < 2000){
            continue;
        }
        //we have an edge that can be reduced
        bcc2reduce[bcc[e]].push_back(es->getEdge_id());
        //cout <<"Candidate Reduce: "<<bcc[e]<<" "<<es->getEdge_id()<<" "<<es->getBw()<<" "<<es->getBd()<<" "<<es->getBd_std()<<endl;
        candidate_reduce_e++;
        candidate_reduce_w+=es->getBw();
        candidate_reduce_l+=es->getBd();
    }
    cout <<"STARTING TRANSITIVE REDUCTION: "<<endl;
    cout << "A total of "<<candidate_reduce_e<<" edges can be reduced, the total weight is "
         << candidate_reduce_w<< " ,whit a total lengh of "<<candidate_reduce_l<<" bases"<<endl;


    //we reduce the edges
    int total_iterations=0;
    //todo: build a threaded version when  the scaffodling graph is large and there are egdes larger than 120kb, as is the case of the low-memory pipeline.
    //step one sort the edges by distance within each BCC
     for(auto bi:bcc2reduce){
        //we sort the edges by distance
         sort(bi.second.begin(),bi.second.end(),[&g](const int a,const int b){return g->get_edge_by_id(a)->getBd() < g->get_edge_by_id(b)->getBd();});
         //we attempt to reduce the edges of the BCC in increasing order and distance
         for(auto e:bi.second){
             auto es=g->get_edge_by_id(e);
             total_iterations+=this->look_for_a_path(lme,g,bp,es,bi.first,bcc2nodes,gToLemon,lemonTog);
             //cout <<"SORTING edges: "<<bi.first<<" "<<es->getEdge()<<" "<<es->getBw()<<" "<<es->getBd()<<" "<<es->getBd_std()<<endl;
         }
     }
    cout << "Number of iterations in edge reduction: "<< total_iterations<<endl;
    bp.clear();
    //we need to change g to return a reduced graph
    auto reduced=new GraphS(g);
    cout <<"REDUCED GRAPH " <<reduced->get_number_edges()<<" ORIGINAL GRAPH "<<g->get_number_edges()<<endl;
    cout <<"ENDING TRANSITIVE REDUCTION"<<endl;
    return reduced;
}


//we reduce the edges
int  TReduction::look_for_a_path(int lme,GraphS* g, ListGraph &bp,EdgeS* e, int bi, hashhash &bcc2nodes,simplehash &gToLemon,simplehash &lemonTog){
    //we obtain the nodes from the edge
    auto u=g->get_node(e->getSource());
    auto v=g->get_node(e->getTarget());
    //init path from u-> source
    auto p=new PathS(u->getId());
    //std::unique_ptr<auto> p = make_unique<PathS>(u->getId());
    //std::unique_ptr<PathS> p( new PathS(u->getId()) );
    //unique_ptr<PathS> p(new PathS(u->getId()));
    //paths cointainers
    vector<PathS*> paths;//paths to look for
    vector<PathS*> vpaths;//paths founds as valid in length with respect to the reduced edge
    //vector<unique_ptr<PathS* > > paths;
    //vector<unique_ptr<PathS* > > vpaths;
    bool ulong=0;
	if(e->getBd() >= lme )
	ulong=1;

    int iter=0;
    paths.push_back(move(p));
    //while there are elements in the stack
    while(!paths.empty()){
        //we got a copy of the last added path
        auto pc=paths.back();
        iter++;
	//for ultra long edges we perform upto 20Million iterations, to keep the combinatory under control.
    //if(iter >=20000000)
	if(ulong ==true and iter >=20000000)
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
                    //for ultra long edges we just consider the fisrt 100 valid paths, to reduce the search
                    if(vpaths.size() >=100 and ulong == true)
                        break;
		    //we look upto 1000 paths, ussually is a repeat not masked, again to keep the combinatories under control
                    if(vpaths.size() >=1000)
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
                //not in BCC
                if(bcc2nodes[bi].count(nn) == 0){
                    continue;
                }
                //not a repeat node or a shorter one this was a critical desition increase NGA50 to 5Mb
                auto gn=g->get_node(nn);
                if(gn->is_short() or gn->is_repeat()){
                    continue;
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
    //todo: change this by getEdgefromID
    //funtion to obtain the edge to remove from the graph
        auto getEdge = [&bp](lemon::ListGraph::Node n1, lemon::ListGraph::Node n2)
        {
            for (lemon::ListGraph::IncEdgeIt edgeIt(bp, n1);
                 edgeIt != lemon::INVALID; ++edgeIt)
            {
                if (bp.oppositeNode(n1, edgeIt) == n2) return edgeIt;
            }
            return lemon::ListGraph::IncEdgeIt(lemon::INVALID);
        };

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
                //todo:choose the one wiht more hits?
                //pv->print_path();
                /*if(pv->get_path_size() >= stmp->get_path_size()){
                    if(pv->getClen() > stmp->getClen()){
                        stmp=pv;
                    }
                }*/
                //we choose the path with more support in the long reads
                if(pv->get_path_hits() >= stmp->get_path_hits()){
                    stmp=pv;
                }
                validated++;
            }
        }

        if(validated>0){
            cout << "STMP PATH: "<<validated <<endl;
            //todo save this to a pathdb for later use
            stmp->print_path();
            //cout <<endl;
            //we select the edge to remove from the lemon graph
            auto del_edge=getEdge(bp.nodeFromId(gToLemon[u->getId()]),bp.nodeFromId(gToLemon[v->getId()]));

            if(!bp.valid(del_edge)){
                cout << "Invalid edge"<<u->getId()<<" "<<v->getId()<<" "<<e->getEdge()<<endl;
                exit(1);
            }
            //increase the w of the paths edges
            stmp->reduce_edge(g,e);
            //we remove the edge from the lemon graph
            bp.erase(del_edge);
            //we save the reduced edge and the selected path
            //reduced_edges[e->getEdge()]=stmp->getPlen();
            //we save the distance and the var of the path
            reduced_edges[e->getEdge()].first=stmp->getPlen();
            reduced_edges[e->getEdge()].second=stmp->getPvar();
        }else{
            cout <<"NO paths founds for edge: "<<bi<<" "<<e->getEdge()<<" "<<e->getBw()<<" "<<e->getBd()<<" "<<e->getBd_std()<<" NODES(e): "<<u->getId()<<" "<<v->getId()<<endl;
        }
    }else{
        cout <<"NO paths founds for edge: "<<bi<<" "<<e->getEdge()<<" "<<e->getBw()<<" "<<e->getBd()<<" "<<e->getBd_std()<<" NODES(e): "<<u->getId()<<" "<<v->getId()<<endl;
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

void TReduction::print_reduced_edges() {
    cout << "Reduced edge(ID) and path_length"<<endl;
    for(auto e:reduced_edges){
        cout << e.first<<" "<<e.second.first<<" "<<e.second.second<<endl;
    }
}
