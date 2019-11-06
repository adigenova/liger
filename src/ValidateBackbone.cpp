//
// Created by Alex Digenova on 9/19/18.
//

#include "ValidateBackbone.h"


GraphS* ValidateBackbone::Check_Backbone(int min_frag_len,int min_long_edge, GraphS *cover, hash2pair64 &reduced_edges,simplehash &circular_nodes) {
    //STEP 1: we create the scaffold lines and save them as paths of the cover graph
    create_backbone_lines(cover,circular_nodes);

    //STEP 2: we  add the contig seqs and the non reduced edges as genomic fragments, edges are added only if they are supported by
    //synthetic libraries larger than min_frag_len
    uint64_t scaff_length=0;
    int max_scaff=0;
    for (auto p:lines) {
       // p->print_path();
        p->compute_pos_in_line(cover);
        //add the non-reduced edges and contigs as fragments to compute the coverage
        p->add_edges_contigs_to_line(cover,min_frag_len);
        if(p->getPlen()>=min_frag_len){
            scaff_length+=p->getPlen();
            if(p->getPlen() >= max_scaff)
                max_scaff=p->getPlen();
        }
    }

    //STEP 3: we add the physical coverage from reduced edges to the lines
    uint32_t source,target;
    int sp,tp;
    int total_edges_added=0;
    int total_edges_att=0;
    uint64_t total_path_coverage=0;
    for(auto e:reduced_edges){
        if(e.second.first < min_frag_len){
            continue;
        }
        total_edges_att++;
        //this are the reduced edges
        depairf(e.first,source,target);
        sp=-1;
        tp=-1;
        if(nodes2lines.count(source))
            sp=nodes2lines[source];

        if(nodes2lines.count(target))
            tp=nodes2lines[target];
        //both nodes are in the same SLine
        if(sp>=0 and sp==tp) {
         //we attempt to add the reduced edge as a fragment to the line
            auto added=lines[sp]->add_fragment_to_line(e.first,e.second.first,e.second.second);
            if(added){
                total_edges_added+=added;
                total_path_coverage+=e.second.first;
            }


        }
    }

    //we print some stats of the process
    cout << "Total Reduced edges >= "<<min_frag_len<<" (bp) = "<<total_edges_att<<endl;
    cout <<"A total of "<<total_edges_added<<" ("<<100 * ((float)total_edges_added/(float)total_edges_att)<<"%) reduced edges were converted to fragments"<<endl;
    cout << "The Physical Coverage represented in the Paths is: "<<(float)total_path_coverage/(float)scaff_length<<endl;
    cout << "Scaffolds length larger than >= "<< min_frag_len<<" = "<<scaff_length<<endl;
    cout << "Max scaffold length ="<< max_scaff<<endl;

    //STEP 4: we compute the LQI and determine posibble miss-assemblies at the level of edges, contigs are already good quality?
   // int total_deleted_edges=0;
    //we reserve the vector
   // vector<int> cov(max_scaff+1,0);
    //cov.reserve(max_scaff+1);
    //std::fill(cov.begin(),cov.end(),0);
   /* for (auto p:lines) {
        //we compute the LQI segments and breaks the scaffols lines if necessary, we break when physical coverage drops to 0 and we mark the nearest
        //edge as deleted
        total_deleted_edges+=p->LQI2break_lines(cover,min_long_edge,cov);
    }
    cov.clear();
    //we delete the containers of the class
   *//* nodes2lines.erase(nodes2lines.begin(),nodes2lines.end());
    lines.erase(lines.begin(),lines.end());*//*
    lines.clear();
    nodes2lines.clear();*/
    //we validate the backbone here
    // 1. we sort the lines

    auto total_deleted_edges=this->validate_backbone(max_scaff+1,min_long_edge,cover);
    //we clear the vector before it is destroyed
    if(total_deleted_edges > 0){
       // cout << "A total of " << total_deleted_edges <<" edges were deleted in the validation step"<<endl;
        //we copy the edges and return the new graph
        auto validated=new GraphS(cover);
        cout <<"VALIDATED GRAPH " <<validated->get_number_edges()<<" ORIGINAL GRAPH "<<cover->get_number_edges()<<endl;
        return validated;
    }else{
        //we don't deleted the edges
        cout <<"VALIDATED GRAPH  " <<cover->get_number_edges()<<" ORIGINAL GRAPH "<<cover->get_number_edges()<<" are equals"<<endl;
        return cover;
    }

}


int ValidateBackbone::validate_backbone(int max_scaff_size,int min_long_edge,GraphS *g) {

    clock_t begin = clock();
    //Step 1: we sort the fragments and dedges inside each Sline;
    for (auto p:lines)
        p->sort_fragments();

    //we create the container to validate each line
    int total_deleted_edges=0;
    uint64_t total_exam_bases=0;
    uint64_t total_good_bases=0;
    uint64_t  avg_fragment_length=0;
    uint64_t  total_fragments=0;
    //we reserve the vector
    vector<uint8_t> cov(max_scaff_size+1,0); //max value of coverage is 2^8=256
    //results vector
    vector<LQI> results;

    for(auto p:lines) {
        //clock_t begin = clock();
        //we fill the vector of counts wiht the fragments
        int maxpos = p->get_line_size();
        //cout << "Maxpos of fragment:" << maxpos << endl;
        assert(maxpos > 0 and maxpos < cov.size());
        total_exam_bases+=maxpos;
        std::fill(cov.begin(), cov.end(), 0);//takes 95 seconds

        for (auto f:p->get_fragments()) {
            int stop = f.pos + f.d;
            avg_fragment_length += f.d;
            //we check the the fragment don't exceed the size of the container
            //cout <<"BUG: "<<maxpos<<" "<<stop<<" "<<f.pos<<" "<<f.d<<" "<<f.contig<<" "<<f.edge<<" "<<f.path<<endl;
            assert(stop <= maxpos and stop >=0);
            total_fragments++;//this count every fragment
            //we store upto a maximum coverage not necessary more
            for (auto it = cov.begin() + f.pos, end = cov.begin() + stop;
                 it != end; ++it) { if(*it < 250){*it += 1;} }//takes 71.5984 secs
        }

        //variables to determine LQI segments
        uint32_t start = 0;
        uint32_t stop = 0;
        uint32_t current = 0;
        uint32_t minimum_cov = 1;//this could be a variable, for the moment 0 is the target

        //uint32_t good_bases = 0;
        //we detect the LQI segments
        for (uint32_t i = 0; i < maxpos; ++i) {
            //cout << "FC: "<<i<<" "<<cov[i]<<endl;
            //we detect the LQI segments
            if (cov[i] < minimum_cov) {
                if (current == 0) {
                    current = i;
                    start = i;
                }
                //allow a max dist of 10
                if (i < current + 10) {
                    current = i;
                    stop = i;
                } else {
                    //is a new interval
                    //cout << "LQI " << current_contig<<" "<<start <<" "<<stop<<" "<<stop-start<<endl;
                    //we have to determine the edgeID to deleted
                    LQI r;
                    //r.ctgid=ctgid;
                    r.start = start;
                    r.stop = stop;
                    //the mis assembly is at the interior of the contig length
                    r.type = 2;//means internal type
                    //header type
                    //we store the current LQI in the array
                    results.push_back(r);
                    current = 0;
                    start = i;
                    stop = i;
                }
            } else {
                //are bases that are supportted by at lest two fragments
                //good_bases++;
                total_good_bases++;
            }

        }

        //we clear the buffer
        //fill(cov.begin(),cov.begin()+maxpos,0);//takes 95 seconds
        //cov.erase(cov.begin(),cov.end());
        std::fill(cov.begin(), cov.end(), 0);

        /* for(auto lqi:results){
             cout <<"LQI="<< lqi.edgeid<<" "<<lqi.start<<" "<<lqi.stop<<" "<<lqi.type<<endl;
         }*/
        //if there are dangereusse edges and LQI segments
        if (results.size() > 0) {
            //we use a simple array to detect the
            for (auto lqi:results) {
                for (auto d:p->get_fragments_dangereus()) {
                    //we check if there is an overlap with the danger edges
                    if ((lqi.start >= d.pos and lqi.start <= d.pos + d.d) or
                        (lqi.stop >= d.pos and lqi.stop <= d.pos + d.d)) {
                        auto edge = g->get_edge(d.edgeid);
                        /* cout <<"POTENTIAL QUIMERIC EDGES ARRAY: EID="<<d.edgeid<<" LQI_S="<<lqi.start<<" LQI_E="<<lqi.stop<<" EID="
                              <<d.edgeid<<" E_S="<<d.pos<<" E_E="<<d.pos+d.d<<" "<<edge->getBw()<<" NLR="<<edge->getLong_reads().size()<<endl;*/

                        if (edge->getLong_reads().size() < min_long_edge) {
                            edge->setDeleted(true);
                            cout << "DELETED EDGE ARRAY: EID=" << d.edgeid << " LQI_S=" << lqi.start << " LQI_E="
                                 << lqi.stop << " EID="
                                 << d.edgeid << " E_S=" << d.pos << " E_E=" << d.pos + d.d << " " << edge->getBw()
                                 << " NLR=" << edge->getLong_reads().size() << endl;
                            total_deleted_edges++;
                            continue;//we pass to the next resultt
                        }
                    }
                }
            }

        }
        //we clear the LQI detected
        results.clear();
    }
    //here we print the stats
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    avg_fragment_length = int(avg_fragment_length / total_fragments);
    cout << "Total lines examinated "<< lines.size()<<endl;
    cout << "Total fragments " << total_fragments << endl;
    cout << "Average fragment length " << avg_fragment_length << endl;
    cout << "Total Examinated bases " << total_exam_bases << endl;
    cout << "Total Validated bases by at least " << 1 << " fragments " << total_good_bases << endl;
    cout << "Percentaje of validated bases  " << ((float) total_good_bases / (float) total_exam_bases * 100) << endl;
    cout << "A total of " << total_deleted_edges <<" edges were deleted in the validation step"<<endl;
    cout << "Time spent in lines validation " << elapsed_secs << " secs " << endl;

    return total_deleted_edges;
}



void ValidateBackbone::create_backbone_lines(GraphS *g,simplehash &circular_nodes){
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
    ListGraph::EdgeMap<int> eid(bp);
    simplehash lemonE2g;
    for(int k=0; k<g->get_number_edges(); k++){
        //we have to add the dummy edges and the selected by the greedy method
        auto es=g->get_edge(k);
        //we mark the edge as live
        es->setDeleted(false);
        //cout << "IN GRAPH:"<<es->getEdge()<<" "<<es->getBw()<<" "<<es->is_dummy()<<endl;
        auto u=g->get_node(es->getSource());
        auto v=g->get_node(es->getTarget());
        //is an input edge of the matching should not be removed
        auto leftLemonNode = bp.nodeFromId(gToLemon[u->getId()]);
        auto rightLemonNode = bp.nodeFromId(gToLemon[v->getId()]);
        auto edge = bp.addEdge(leftLemonNode, rightLemonNode);
        eid[edge]=es->getEdge_id();
        lemonE2g[bp.id(edge)]=es->getEdge_id();
    }

    //we compute the connected components of the graph
    ListGraph::NodeMap<int> cc(bp);
    int ncc=connectedComponents(bp,cc);
    cout << "Number of CC "<< ncc << endl;
    //hash2vec nodes2cc;
    simplehash scc;//store the size of each cc
    hash2vec headers;

    //we iterate the nodes of the graph to add them to the cc and search for a header/tail of each scf
    for (ListGraph::NodeIt n(bp); n != INVALID; ++n) {
        int cnt = 0;
        for (ListGraph::IncEdgeIt e(bp, n); e != INVALID; ++e) {
            cnt++;
        }
        //cout << "deg(" << bp.id(n) << ") = " << cnt <<" " << cc[n] <<endl;
        //we found a header
        if(cnt == 1){
            headers[cc[n]].push_back(bp.id(n));
        }
        scc[cc[n]]++;
    }

    //object to perform the path search
    Bfs<ListGraph> bfs(bp);
    int current_path=0;
    //we iterate the cc
    for(auto c:headers){
        assert(c.second.size() == 2);
        //we validate long contigs
        if(scc[c.first] > 2){
            //we obtain the shortest path from the bfs
            auto h=bp.nodeFromId(c.second[0]);
            auto t=bp.nodeFromId(c.second[1]);
            //we found the path
            bfs.run(h,t);
            vector<int> path;
            path.push_back(bp.id(t));
            auto prev = bfs.predNode(t);
            while (prev != INVALID) {
                path.push_back(bp.id(prev));
                prev = bfs.predNode(prev);
            }
            reverse(path.begin(),path.end());
            //we have to translate the path to the graph nodes
            //auto p=new PathS(lemonTog[path[0]]);
            auto p=new SLines(lemonTog[path[0]]);
            nodes2lines[lemonTog[path[0]]]=current_path;
            //we save is the line is circular or not
            if(circular_nodes.count(lemonTog[path[0]]))
                p->set_as_circular(lemonTog[path[0]]);

            for (int i = 1; i <path.size() ; ++i) {
                p->add_node(lemonTog[path[i]]);
                nodes2lines[lemonTog[path[i]]]=current_path;
                if(circular_nodes.count(lemonTog[path[i]]))
                    p->set_as_circular(circular_nodes[lemonTog[path[i]]]);
            }
            //we print the path
            p->compute_path_length(g);
            //we save the path in the graph
            //cout << "Printing graph lines"<<endl;
            //p->print_path();
            lines.push_back(p);
            current_path++;
        }
    }
    //we distroy the tmp map
    bp.clear();

}


