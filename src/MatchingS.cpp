//
// Created by Alex Digenova on 7/10/18.
//

#include "MatchingS.h"

GraphS* MatchingS::matching_cover(GraphS* g) {

    //Step 1: We create the lemon graph for running the matching without the dummy edges
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
    ListGraph::EdgeMap<int64_t> eid(bp);
    ListGraph::EdgeMap<int64_t> wt(bp);
    int64_t totalw=0;
    int64_t totale=0;
    int64_t totalnd=0;
    for(int k=0; k<g->get_number_edges(); k++){
        //we have to add the dummy edges and the selected by the greedy method
        auto es=g->get_edge(k);
        //we enchure that all edges will be considered in the matching
        es->setDeleted(false);
        //we skypt the dummy edges or low supported edges
        if(es->is_dummy() or es->getBw() < 5){
            continue;
        }
        //cout << "IN GRAPH:"<<es->getEdge()<<" "<<es->getBw()<<" "<<es->is_dummy()<<endl;
        auto u=g->get_node(es->getSource());
        auto v=g->get_node(es->getTarget());
        totalnd++;
        if(u->is_repeat() or v->is_repeat()){
            continue;
        }
        //we dont scaff if one of the two contigs is short due to low coverage
        if(u->is_short() or v->is_short()){
            continue;
        }
        //is an input edge of the matching should not be removed
        auto leftLemonNode = bp.nodeFromId(gToLemon[u->getId()]);
        auto rightLemonNode = bp.nodeFromId(gToLemon[v->getId()]);
        auto edge = bp.addEdge(leftLemonNode, rightLemonNode);
        eid[edge]=es->getEdge_id();
        wt[edge]=es->getBw();//for a dummy edge the wt is 0
        totalw+=wt[edge];
        totale++;
    }

    //Step 2: we run the matching cover
    //we create the matching object
    MaxWeightedMatching<ListGraph, ListGraph::EdgeMap<int64_t> > mwpm(bp, wt);
    mwpm.run();
    cout << "Matching Done" << endl;
    //Step 3: we remove the edges not present in the matching
    int64_t n_edge_match=0;
    simplehash matching_edges;
    for (ListGraph::EdgeIt e(bp); e != INVALID; ++e) {
        auto eg=g->get_edge_by_id(eid[e]);
        if(mwpm.matching(e)) {
            n_edge_match++;
            matching_edges[eid[e]]++;
        }else{
            bp.erase(e);
        }
    }
    cout << "Total Edges :"<<totalnd<<" input M "<<totale<<" Edges in matching: "<<n_edge_match<<" MCW: " << mwpm.matchingWeight() <<" Total W:"<<totalw<<" % "<<(float)mwpm.matchingWeight()/totalw << endl;

    //Step 4 : we check for cicles in the covered graph
    //we add the dummy edges to build the graph and check for cicles
    for(int k=0; k<g->get_number_edges(); k++){
        //we have to add the dummy edges and the selected by the greedy method
        auto es=g->get_edge(k);
        //is not a dummy edge and was not selected
        if(!es->is_dummy()){
            //the edge is not part of the matching
            if(matching_edges.count(es->getEdge_id())==0){
                es->setDeleted(true);
            }
            continue;
        }

        auto u=g->get_node(es->getSource());
        auto v=g->get_node(es->getTarget());

        auto leftLemonNode = bp.nodeFromId(gToLemon[u->getId()]);
        auto rightLemonNode = bp.nodeFromId(gToLemon[v->getId()]);
        auto edge = bp.addEdge(leftLemonNode, rightLemonNode);
        eid[edge]=es->getEdge_id();
        wt[edge]=es->getBw();//for a dummy edge the wt is 0
    }

    //we check for simple cycles at this point
    if(!acyclic(bp)){
        cout << "the graph has cicles and we should remove some edges"<<endl;
        //we compute the biconnected components of the graph to remove the remaining cicles
        ListGraph::EdgeMap<int64_t> bcc(bp);
        int nbcc=biNodeConnectedComponents(bp,bcc);
        //cout << "Number of nbcc "<<nbcc<<endl;
        unordered_map<int64_t, ListGraph::Edge> bcc2remove;
        simplehash bccsize;
        //I have to save the lowest w edge for each BCC
        for (ListGraph::EdgeIt e(bp); e != INVALID; ++e){
           // cout <<"EDGES in cicles "<<eid[e]<<" "<<wt[e]<<" "<<bcc[e]<<endl;
            //is not a dummy edge
            if(wt[e] > 0){
                //we save a value in the bcc2remove
                if(bcc2remove.count(bcc[e]) == 0){
                    bcc2remove[bcc[e]]=e;
                }else{
                    auto stored=bcc2remove[bcc[e]];
                    //the w of stored is larger than the new
                    if(wt[stored] > wt[e]){
                        bcc2remove[bcc[e]]=e;
                    }
                }
            }
            //store the size of the BCC
            bccsize[bcc[e]]++;
        }

        //we iterate the edge to remove from the graph
        for(auto bi:bcc2remove){
            //is not a cicle, the min size of a cicle is 2 edges
            if(bccsize[bi.first] < 2){
                continue;
            }
            auto el=bi.second;
            cout << "Deleting edge: "<<bi.first<<" "<<bccsize[bi.first]<<" "<<wt[el]<<" "<<eid[el]<<endl;
            auto eg=g->get_edge_by_id(eid[el]);

            //we mark the edge has deleted in our graph
            eg->setDeleted(true);
            //we mark the nodes to save that this scf is circular
            circular[eg->getTarget()]=eg->getBd();
            circular[eg->getSource()]=eg->getBd();
            //we remove the edges
            bp.erase(el);
        }
        //we check that the cicles have been removed
        if(!acyclic(bp)){
            cout << "ERROR in Matching COVER: at this point we removed all the cycles, however we found cicles still"<<endl;
            exit(1);
        }
    }else{
        cout << "there is no cicles in the graph"<<endl;
    }
    //we erase the lemon graph
    bp.clear();
    //we perform the copy of the graph
    //step 5: we return the matching covered new graph
    auto mcover=new GraphS(g);
    cout <<"Matching COVER GRAPH " <<mcover->get_number_edges()<<" ORIGINAL GRAPH "<<g->get_number_edges()<<endl;
    return mcover;
}

GraphS *MatchingS::matching_cover(GraphS *g, float max_repeat_factor, int short_ctg_length) {

    //Step 1: We create the lemon graph for running the matching without the dummy edges
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
    ListGraph::EdgeMap<int64_t> eid(bp);
    ListGraph::EdgeMap<int64_t> wt(bp);
    int64_t totalw=0;
    int64_t totale=0;
    int64_t totalnd=0;
    for(int k=0; k<g->get_number_edges(); k++){
        //we have to add the dummy edges and the selected by the greedy method
        auto es=g->get_edge(k);
        //we enchure that all edges will be considered in the matching
        es->setDeleted(false);
        //we skypt the dummy edges or low supported edges
        if(es->is_dummy() or es->getBw() < 5){
            continue;
        }
        //cout << "IN GRAPH:"<<es->getEdge()<<" "<<es->getBw()<<" "<<es->is_dummy()<<endl;
        auto u=g->get_node(es->getSource());
        auto v=g->get_node(es->getTarget());
        totalnd++;

        //we skip short contigs
        if(u->getnode_len() < short_ctg_length or v->getnode_len() < short_ctg_length){
            continue;
        }

        //we fetch the contigs and skip the repetitive ones.
        auto cu=g->getContigs()->getcontig(u->getCtg_id());
        auto cv=g->getContigs()->getcontig(v->getCtg_id());
        //we skip repeat contigs
        if(cu.coverage > max_repeat_factor or cv.coverage > max_repeat_factor){
            continue;
        }

        //is an input edge of the matching should not be removed
        auto leftLemonNode = bp.nodeFromId(gToLemon[u->getId()]);
        auto rightLemonNode = bp.nodeFromId(gToLemon[v->getId()]);
        auto edge = bp.addEdge(leftLemonNode, rightLemonNode);
        eid[edge]=es->getEdge_id();
        wt[edge]=es->getBw();//for a dummy edge the wt is 0
        totalw+=wt[edge];
        totale++;
    }

    //Step 2: we run the matching cover
    //we create the matching object
    MaxWeightedMatching<ListGraph, ListGraph::EdgeMap<int64_t> > mwpm(bp, wt);
    mwpm.run();
    cout << "Matching Done" << endl;
    //Step 3: we remove the edges not present in the matching
    int64_t n_edge_match=0;
    simplehash matching_edges;
    for (ListGraph::EdgeIt e(bp); e != INVALID; ++e) {
        auto eg=g->get_edge_by_id(eid[e]);
        if(mwpm.matching(e)) {
            n_edge_match++;
            matching_edges[eid[e]]++;

        }else{

            bp.erase(e);
        }
    }
    cout << "Total Edges :"<<totalnd<<" input M "<<totale<<" Edges in matching: "<<n_edge_match<<" MCW: " << mwpm.matchingWeight() <<" Total W:"<<totalw<<" % "<<(float)mwpm.matchingWeight()/totalw << endl;

    //Step 4 : we check for cicles in the covered graph
    //we add the dummy edges to build the graph and check for cicles
    for(int k=0; k<g->get_number_edges(); k++){
        //we have to add the dummy edges and the selected by the greedy method
        auto es=g->get_edge(k);
        //is not a dummy edge and was not selected
        if(!es->is_dummy()){
            //the edge is not part of the matching
            if(matching_edges.count(es->getEdge_id())==0){
                es->setDeleted(true);
            }
            continue;
        }

        auto u=g->get_node(es->getSource());
        auto v=g->get_node(es->getTarget());

        auto leftLemonNode = bp.nodeFromId(gToLemon[u->getId()]);
        auto rightLemonNode = bp.nodeFromId(gToLemon[v->getId()]);
        auto edge = bp.addEdge(leftLemonNode, rightLemonNode);
        eid[edge]=es->getEdge_id();
        wt[edge]=es->getBw();//for a dummy edge the wt is 0
    }

    //we check for simple cycles at this point
    if(!acyclic(bp)){
        cout << "the graph has cicles and we should remove some edges"<<endl;
        //we compute the biconnected components of the graph to remove the remaining cicles
        ListGraph::EdgeMap<int64_t> bcc(bp);
        int nbcc=biNodeConnectedComponents(bp,bcc);
        //cout << "Number of nbcc "<<nbcc<<endl;
        unordered_map<int64_t, ListGraph::Edge> bcc2remove;
        simplehash bccsize;
        //I have to save the lowest w edge for each BCC
        for (ListGraph::EdgeIt e(bp); e != INVALID; ++e){
            // cout <<"EDGES in cicles "<<eid[e]<<" "<<wt[e]<<" "<<bcc[e]<<endl;
            //is not a dummy edge
            if(wt[e] > 0){
                //we save a value in the bcc2remove
                if(bcc2remove.count(bcc[e]) == 0){
                    bcc2remove[bcc[e]]=e;
                }else{
                    auto stored=bcc2remove[bcc[e]];
                    //the w of stored is larger than the new
                    if(wt[stored] > wt[e]){
                        bcc2remove[bcc[e]]=e;
                    }
                }
            }
            //store the size of the BCC
            bccsize[bcc[e]]++;
        }

        //we iterate the edge to remove from the graph
        for(auto bi:bcc2remove){
            //is not a cicle, the min size of a cicle is 2 edges
            if(bccsize[bi.first] < 2){
                continue;
            }
            auto el=bi.second;
            cout << "Deleting edge: "<<bi.first<<" "<<bccsize[bi.first]<<" "<<wt[el]<<" "<<eid[el]<<endl;
            auto eg=g->get_edge_by_id(eid[el]);

            //we mark the edge has deleted in our graph
            eg->setDeleted(true);
            //we mark the nodes to save that this scf is circular
            circular[eg->getTarget()]=eg->getBd();
            circular[eg->getSource()]=eg->getBd();
            //we remove the edges
            bp.erase(el);
        }
        //we check that the cicles have been removed
        if(!acyclic(bp)){
            cout << "ERROR in Matching COVER: at this point we removed all the cycles, however we found cicles still"<<endl;
            exit(1);
        }
    }else{
        cout << "there is no cicles in the graph"<<endl;
    }
    //we erase the lemon graph
    bp.clear();
    //we perform the copy of the graph
    //step 5: we return the matching covered new graph
    auto mcover=new GraphS(g);
    cout <<"Matching COVER GRAPH " <<mcover->get_number_edges()<<" ORIGINAL GRAPH "<<g->get_number_edges()<<endl;
    return mcover;

}
