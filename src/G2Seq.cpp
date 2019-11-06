//
// Created by Alex Digenova on 7/12/18.
//

#include "G2Seq.h"

//we create an AGP file the Illumina backbone
void G2Seq::graph2agp(GraphS *g, string file) {


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


    Bfs<ListGraph> bfs(bp);
    //we create the ouput file
    ofstream agp;
    agp.open(file);
    //rcon << "#Ctg_id Name length"<<endl;


    //we iterate the cc
    for(auto c:headers){
        assert(c.second.size() == 2);
        //we can create the AGP file
        //cout <<" Headers cc("<<c.first <<") size "<<scc[c.first]<<" H("<<c.second[0]<<") T("<<c.second[1]<<")"<<endl;
        if(scc[c.first] > 2){
            //we obtain the shortest path from the bfs
            auto h=bp.nodeFromId(c.second[0]);
            auto t=bp.nodeFromId(c.second[1]);
            //we found the path
            bfs.run(h,t);
            vector<int> path;
            //bfs style of shortest path
            //cout << bp.id(t);
            path.push_back(bp.id(t));
            auto prev = bfs.predNode(t);
            while (prev != INVALID) {
              //  cout << "<-" << bp.id(prev);
                path.push_back(bp.id(prev));
                prev = bfs.predNode(prev);
            }
            reverse(path.begin(),path.end());
            //we print the path in g
            int start=1;
            int counter=1;
            for (int i=0; i<path.size()-1; i++) {
                //cout << g->get_node(lemonTog[p])->getName() <<"->";
                auto u =g->get_node(lemonTog[path[i]]);
                auto v =g->get_node(lemonTog[path[i+1]]);
                //is dummy edge
                if(u->getCtg_id() == v->getCtg_id()){
                    //we got the contig from the node
                    auto contig=g->getContigs()->getcontig(u->getCtg_id());
                    string strand="-";
                    //means that we found a match
                    //if(v->getName().find("+")!=std::string::npos)
                    if(v->getId()==contig.header)
                        strand="+";
                    agp<<join("\t",vector<string>({"SCF"+to_string(c.first+1),to_string(start),to_string(start+contig.length),to_string(counter),"W",contig.namectg,"1",to_string(contig.length),strand}))<<endl;
                    counter++;
                    start+=contig.length;
                }else{
                    //is a real edge
                    auto e=g->get_edge(u,v);
                    int gap=e->getBd();
                    //we dont merge contigs in the backbone
                    if(gap < 20){
                        gap=20;
                    }
                    agp<<join("\t",vector<string>({"SCF"+to_string(c.first+1),to_string(start),to_string(start+gap),to_string(counter),"N",to_string(gap),"scaffold","yes","paired-end"}))<<endl;
                    counter++;
                    start+=gap;
                }
            }

        }else{
            //is a singleton, we should output the name of the scaffold
            //auto u=bp.nodeFromId(c.second[0]);
            //print  AGP join("\t","SCF".$c,1,1+$contigs->{$tctg}->{len},"1","W",$tctg, 1,$contigs->{$tctg}->{len},"+")."\n";
            auto u =g->get_node(lemonTog[c.second[0]]);
            auto contig=g->getContigs()->getcontig(u->getCtg_id());
            agp<<join("\t",vector<string>({"SCF"+to_string(c.first+1),"1",to_string(contig.length+1),"1","W",contig.namectg,"1",to_string(contig.length),"+"}))<<endl;
        }
    }

    //we distroy the tmp map
    bp.clear();
    //we close the file
    agp.close();
}

void G2Seq::graph2seq(GraphS *g, string prefix, bool singletons) {

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


    Bfs<ListGraph> bfs(bp);
    //we create the ouput file
    ofstream fasta;
    fasta.open(prefix+".asm.wengan.fasta");
    //rcon << "#Ctg_id Name length"<<endl;
    ofstream bedfile;
    bedfile.open(prefix+".asm.wengan.bed");
    //SC = super contig
    //UC = unique contig
    //we iterate the cc
    for(auto c:headers){
        assert(c.second.size() == 2);
        //we can create the AGP file
        //cout <<" Headers cc("<<c.first <<") size "<<scc[c.first]<<" H("<<c.second[0]<<") T("<<c.second[1]<<")"<<endl;
        if(scc[c.first] > 2){
            //we obtain the shortest path from the bfs
            auto h=bp.nodeFromId(c.second[0]);
            auto t=bp.nodeFromId(c.second[1]);
            //we found the path
            bfs.run(h,t);
            vector<int> path;
            //bfs style of shortest path
            //cout << bp.id(t);
            path.push_back(bp.id(t));
            auto prev = bfs.predNode(t);
            while (prev != INVALID) {
                //  cout << "<-" << bp.id(prev);
                path.push_back(bp.id(prev));
                prev = bfs.predNode(prev);
            }
            reverse(path.begin(),path.end());
            //we print the path in g
            int start=1;
            int counter=1;
            string scfseq="";
            for (int i=0; i<path.size()-1; i++) {
                //cout << g->get_node(lemonTog[p])->getName() <<"->";
                auto u =g->get_node(lemonTog[path[i]]);
                auto v =g->get_node(lemonTog[path[i+1]]);
                //is contig edge
                if(u->getCtg_id() == v->getCtg_id()){
                    //we got the contig from the node
                    auto contig=g->getContigs()->getcontig(u->getCtg_id());
                    string strand="-";
                    //means that we found a match
                    //if(v->getName().find("+")!=std::string::npos)
                    if(v->getId()==contig.header) {
                        strand = "+";
                    }

                    counter++;
                    start+=contig.length;
                    auto beg=scfseq.length();
                    scfseq+= strand.compare("+") == 0 ? contig.seq : revcomp(contig.seq) ; //we add the contig in the correct orientation
                    auto end=scfseq.length();

                    bedfile <<join("\t",vector<string>({"WSC"+to_string(c.first+1),"CTG",to_string(beg),to_string(end),strand,to_string(contig.length),contig.namectg}))<<endl;

                }else{
                    //is a real edge
                    auto e=g->get_edge(u,v);
                    int gap=e->getBd();
                    //we ask if the edge has sequence
                    if(e->isEdge_has_seq()){
                            string gseq=e->getEdge_gapseq();
                            //transform(gseq.begin(), gseq.end(), gseq.begin(), ::tolower);
                            auto beg=scfseq.length();
                            scfseq += e->getSource() == u->getId() ? gseq : revcomp(gseq);
                            auto end=scfseq.length();
                        bedfile << join("\t",vector<string>({"WSC"+to_string(c.first+1),"GAPS",to_string(beg),to_string(end),e->getSource() == u->getId() ? "+":"-",to_string(gseq.length()),to_string(e->getEdge()),e->getSource() == u->getId() ? gseq : revcomp(gseq)}))<<endl;
                            //cout<<join("\t",vector<string>({"SCF"+to_string(c.first+1),to_string(start),to_string(start+gap),to_string(counter),to_string(gap),e->getEdge_gapseq()}))<<endl;
                    }else{
                        //is an overlap and we have to reduce the current seq to the length of the overlap
                        //cout << "overlap:"<<e->getEdge_overlap()<<endl;
                        auto beg=scfseq.length();
                        scfseq=scfseq.substr(0,scfseq.length()+e->getEdge_overlap());
                        auto end=scfseq.length();
                        bedfile << join("\t",vector<string>({"WSC"+to_string(c.first+1),"GAPO",to_string(beg),to_string(end),e->getSource() == u->getId() ? "+":"-",to_string(e->getEdge_overlap()),to_string(e->getEdge())}))<<endl;
                    }
                    counter++;
                    start+=gap;
                }
            }
            //we print the scaffold sequence
            fasta << ">WSC"<<c.first+1<<endl;
            fasta << seqformat(60,scfseq);

        }else{
            //is a singleton, we should output the name of the scaffold
            if(singletons) {
                auto u = g->get_node(lemonTog[c.second[0]]);
                auto contig = g->getContigs()->getcontig(u->getCtg_id());
                //we report the contig only if it is longer than 5kb and not used for polishing edges
                if (contig.length >= 5000 and !contig.used_in_polishing) {
                    fasta << ">WUC" << c.first + 1 << endl;
                    fasta << seqformat(60, contig.seq);
                    bedfile << join("\t", vector<string>(
                            {"WUC" + to_string(c.first + 1), "CTG", to_string(1), to_string(contig.length), "+",
                             to_string(contig.length), contig.namectg})) << endl;

                }
            }

        }
    }
    cout << prefix+".asm.wengan.fasta file created\n"<<endl;
    //we distroy the tmp map
    bp.clear();
    //we close the file
    fasta.close();
    bedfile.close();
}


