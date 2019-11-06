//
// Created by Alex Digenova on 7/10/18.
//

#include "GreedyS.h"

GraphS* GreedyS::greedy_cover(GraphS* g) {
    //store the edges
    vector<pair<int,int>> edges2w;
    //store the contigs used
    unordered_map<uint32_t ,bool> cused;

    //Step 1: we select the edges to reduce
    for (int i = 0; i <g->get_number_edges() ; ++i) {
        auto e=g->get_edge(i);
        //we skypt the dummy edges or low supported edges
        if(e->is_dummy() or e->getBw() < 5){
            continue;
        }
        //we continue because is not a dummy edge
        auto u=g->get_node(e->getSource());
        auto v=g->get_node(e->getTarget());
        //we dont scaff repetitives contigs
        if(u->is_repeat() or v->is_repeat()){
            continue;
        }
        //we dont scaff if one of the two contigs is short due to low coverage
        if(u->is_short() or v->is_short()){
            continue;
        }
        pair<int,int> ec(e->getEdge_id(),e->getBw());
        edges2w.push_back(ec);
    }

    //step 2: we sort the edges by w
    sort(edges2w.begin(),edges2w.end(),[](const pair<int,int> a,const pair<int,int> b){return a.second > b.second;});
    //step 3, we add the sorted edges by w to the new graph
    int totalw=0;
    int greedyw=0;
    unordered_map<int,bool> selected_edges;
    for(auto te:edges2w){
        //cout << te.first <<" "<<te.second<<endl;
        auto e=g->get_edge_by_id(te.first);
        //we can add the edge to the graph
        if(cused.count(e->getSource()) == 0 and cused.count(e->getTarget()) == 0){
            cused[e->getSource()]=true;
            cused[e->getTarget()]=true;
            greedyw+=te.second;
            selected_edges[te.first]=true;
        }/*else{
            //we delete the
            //e->setDeleted(true);
        }*/
        totalw+=te.second;
    }

    //we performed the matching
    cout << "GREEDY MATHING: " << greedyw <<" "<< totalw<<" "<<(float)greedyw/totalw<<" #EDGES "<<selected_edges.size()<<endl;

    //Step 3: We check for simple cicles and destroy them
    //we create the new graph we need to remove cicles for chains
    //we build the lemon graph to check for cicles
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
    ListGraph::EdgeMap<int> wt(bp);
    int i=0;
    for(int k=0; k<g->get_number_edges(); k++){
        //we have to add the dummy edges and the selected by the greedy method
        auto es=g->get_edge(k);
         es->setDeleted(false);
        //is not a dummy edge and was not selected
        if(!es->is_dummy() and !selected_edges.count(es->getEdge_id())){
            es->setDeleted(true);//this edge was not selected by the greedy matching
            continue;
        }
        //cout << "IN GRAPH:"<<es->getEdge()<<" "<<es->getBw()<<" "<<es->is_dummy()<<endl;
        auto u=g->get_node(es->getSource());
        auto v=g->get_node(es->getTarget());
        auto leftLemonNode = bp.nodeFromId(gToLemon[u->getId()]);
        auto rightLemonNode = bp.nodeFromId(gToLemon[v->getId()]);
        auto edge = bp.addEdge(leftLemonNode, rightLemonNode);
        eid[edge]=es->getEdge_id();
        wt[edge]=es->getBw();//for a dummy edge the wt is 0
        i++;
    }

    //we check for simple cycles at this point
    if(!acyclic(bp)){
        cout << "the graph has cicles and we should remove some edges"<<endl;
        //we compute the biconnected components of the graph to remove the remaining cicles
        ListGraph::EdgeMap<int> bcc(bp);
        int nbcc=biNodeConnectedComponents(bp,bcc);
        cout << "Number of nbcc "<<nbcc<<endl;
        unordered_map<int, ListGraph::Edge> bcc2remove;
        simplehash bccsize;
        //I have to save the lowest w edge for each BCC
        for (ListGraph::EdgeIt e(bp); e != INVALID; ++e){
            //cout <<"EDGES in cicles "<<eid[e]<<" "<<wt[e]<<" "<<bcc[e]<<endl;
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
            //is not a cicle, the min size of a cicle is 3 edges
            if(bccsize[bi.first] < 3){
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
            cout << "ERROR in GREEDY COVER: at this point we removed all the cicles, however we found cicles still"<<endl;
            exit(1);
        }

    }else{
        cout << "there is no cicles in the graph"<<endl;
    }
    //we clear the tmp graph
    bp.clear();
    //step 4: we return the greedy covered new graph
    auto gcover=new GraphS(g);
    cout <<"GREEDY COVER GRAPH " <<gcover->get_number_edges()<<" ORIGINAL GRAPH "<<g->get_number_edges()<<endl;
    return gcover;
}
