//
// Created by Alex Digenova on 10/2/18.
//

#include "BFilling.h"


BFilling::BFilling(GraphS *lines, string fastalongreads, string prefix) {

    //this->reducedG=reduced;
    this->linesG=lines;
    this->filename=fastalongreads;
    this->prefix=prefix;
    //we determine the order of the long reads by determining the order how we will proceess the eddes
    //we build the paths to determine the edge order -> long read order for building the long read bank
    //STEP 1: We build  and order the lines by length
    build_lines();
    //we print the lines order, we can order the paths by length before going to the shorter paths
    //for(auto p:linesP){
      //  p->print_path();
    //}
    //we sort the paths from largest to shortest
    sort( linesP.begin( ), linesP.end( ), [ ]( const SLines* a, const SLines* b)
    {
        return a->getPlen() > b->getPlen();
    });

    //STEP2: having the lines ordered, we can determine the edge ordering and subsequence the long read ordering
    /*for(auto p:linesP){
      p->print_path();
    }*/
    //we determine the edges order

    for(auto p:linesP){
        p->save_edges_id(lines,eids);
    }

    //STEP3: We determine the order of the LONG READS and which of them should be stored in the cache memory to avoid slow reads from disk.
    //save the number of times that a long reads is used in the edges
    unordered_map<uint32_t ,uint32_t > countle;
    //save the order in which the long read should be stored in the LonReadBank
    unordered_map<uint32_t,uint32_t> lorder;
    //store the long reads that should be cached to avoid huge jumps while reading the longread databank
    unordered_map<uint32_t,uint32_t> lcache;
    //current store the pos in the longread bank
    uint32_t  current=0;
    vector<int> longreads;
    for(auto i=0; i<eids.size(); ++i){
        //cout <<i<<" "<<eids[i]<<endl;
        //we assing the long-reads order and we save somes reads in a cached fashion wich are the reads more accessed
        auto e=lines->get_edge(eids[i]);
        //auto longreads=e->getLong_reads();
        e->get_bestn_long_reads(20,longreads);
        for(auto l:longreads){
            if(lorder.count(l) == 0){
                lorder[l]=current;
                current++;
            }else{
                //should be something like 100X or 200X
                // we should cache this read because is involved in distant edges i.e circular scaffodls or edge errors?
                if(current-lorder[l] > 100){
                    //cout <<"CACHE "<<l.first<<" C "<<lorder[l.first]<<" SAW "<<current<<" READS after "<<current-lorder[l.first]<<endl;
                    lcache[l]=current-lorder[l];
                }
            }
            //cout <<"LRO "<<i<<" "<<eids[i]<<" "<<l.first<<" Order "<<lorder[l.first]<<endl;
            countle[l]++;
        }
    }
    //we print some information like the reads that should be cached in order to avoid random access to the disk
    cout <<"Total number of long reads to store "<<lorder.size()<<endl;
    cout << "Total number of long reads that should be cached "<< lcache.size()<<" frac: "<< (float)lcache.size()/(float)lorder.size()<<endl;

    //STEP 4: We create the bank object using the order of the long reads and we cache some long reads in memory to avoid slow disk usage
    clock_t begin = clock();
    cout << "Building the long reads Database"<<endl;
    longreadDB=new BSeqDB(filename,lorder,lcache);
    cout << "Total number of long reads stored in the database : "<<longreadDB->getTotalseqs()<<endl;
    cout << "Total number of long reads cached in memory : "<<longreadDB->getNumberLRCahed()<<endl;
    //we assert that the number of reads stored is equal to the number of reads that we want to store
    assert(lorder.size() == longreadDB->getTotalseqs());
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Time spent building the longread database " << elapsed_secs <<" seconds"<<endl;
    //longreadDB->testRead(lorder);
    //longreadDB->DBdump();
    //longreadDB->DumpCache();
    //End of the constructor of BFilling class, we determined the edge ordering-> long -read ordering and we create the scaffolds lines.
    //STEP 5: we are ready to create the consensus for each edge and fill the gaps using the long-reads an the scaffodling graph
    //we create a log file for output
    //for file writing

    //todo: change the contig name file
    bflog.open(this->prefix+".BFilling.log");
    bflog << "EID Proper_Gap Low_Iden Exact_Overlap Estimated_distance Estimated_distance_std Gap_Consensus_Coverage Coverage_source_end Identity_source_end Coverage_target_end Identity_target_end"<<endl;

}



//multicore way of building edges

//threads variables
struct concurrentedges {
    pthread_mutex_t mut;
    pthread_mutex_t log;//to print some information from threads
    vector<uint64_t> *task;//vector of edges to fill gaps
    BFilling *bf;//object Bfilling
    //Contig *a;
};

void * threaded_edges(void* args) {
    //casting of variables
    auto *tmp = (concurrentedges *) args;
    pthread_mutex_t *mutex = &tmp->mut;
    pthread_mutex_t *log = &tmp->log;
    auto bf = (BFilling *) tmp->bf;
    //auto  contigs=(Contig *) tmp->a;
    auto tasks = (vector<uint64_t >*) tmp->task;
    bool dojob = 1;//
    uint64_t eid = 1;//edge id
    vector<string> lseqs;
    vector<linking>  lin;
    int lid=0;
    string lseq="";
    //fillseqs2(uint64_t eid, int n_reads,vector<string> & lseqs,   , int &lid, string &lseq)
    while (dojob) {
        //we got the long reads passing by the edge
        pthread_mutex_lock(mutex);
        if (tasks->size() > 0) {
            eid = tasks->back();
            tasks->pop_back();
            //we load the long_reads needed here, because after is not possible
            //bf->fillseqs(eid,15,lseqs);
            bf->fillseqs2(eid,20,lseqs,lin,lid,lseq);
        } else {
            dojob = 0;
        }
        pthread_mutex_unlock(mutex);
        //we can build the edges in parallel
        if (dojob) {
            //bf->edge2consensus(eid,lseqs);
            bf->edge2consensus2(eid,lseqs,lin,lid,lseq);
            //we print some information from
            pthread_mutex_lock(log);
            //we print some information from the edge to a log file
            bf->print_cns_info(eid);
            pthread_mutex_unlock(log);
            //we clean the local container of sequences
            //lseqs.erase(lseqs.begin(),lseqs.end());
            lseqs.clear();
            lseqs.shrink_to_fit();
            //lin.erase(lin.begin(),lin.end());
            lin.clear();
            lin.shrink_to_fit();
            lseq.clear();
        }

    }
    return NULL;
}

//function for save information related to the gap filled among edges
void BFilling::print_cns_info(uint64_t eid){
    auto e=linesG->get_edge(eid);
    //I have to store more information as the avg identity and coverage of the flanking nodes
    bflog <<eid<<" "<<e->isEdge_proper_gap()<<" "<<e->isEdge_low_identity()<<" "<<e->getEdge_overlap()<<" "<<e->getBd()
          <<" "<<e->getBd_std()<<" "<<e->getEdge_avg_cns()<<" "<<e->getCovs()<<" "<<e->getIdens() <<" "<<e->getCovt()<<" "<<e->getIdent()<<endl;

}

void BFilling::edge2consensus2(uint64_t eid,vector<string> & lseqs, vector<linking> & lin , int &lid, string &lseq) {

    //we are at each edge, we have to create the edge sequence or made the corresponding overlap
    auto e = linesG->get_edge(eid);
    string seq_end_source = "", seq_end_target = "";//
    //STEP 1: we create the consensus sequence of the edge using the lead seq and orient it according to the edge direction
    auto wcns = new WConsensus(lin, lseqs, lseq, lid, 500);
    auto consensus = wcns->consensusfromlead(lin, lseqs, lseq, lid);
    //auto coverage = wcns->get_coverage_cns();
    auto coverage_s=wcns->get_coverage_cns_string();

    assert(coverage_s.length() == consensus.length());
    //we orient the edge sequence accordingly
    if (edge_dir.count(pairf(e->getSource(), e->getTarget())) == 0) {
        consensus = revcomp(consensus);
        //we have to reverse the value of the coverage vector
        //reverse(coverage.begin(), coverage.end());
        reverse(coverage_s.begin(), coverage_s.end());
    }

    //STEP 2: we map the contigs end to the consensus sequence
    get_edgeseq_ends(e, 500, seq_end_source, seq_end_target);

    //STEP 3: revcomp leadseq
    bool leadrev=false;
    for(auto l:lin) {
        if (l.lonread_id == lid) {
            leadrev = l.switched;
            break;
        }
    }
    //STEP:4 we merge the contigs ends using the consensus sequence
    CJoiner joiner;//object that compute the best fit of the contigs-ends in the edge
    joiner.contigs_best_fit(seq_end_source, seq_end_target, consensus, coverage_s, lseq, e,
                         edge_dir.count(pairf(e->getSource(), e->getTarget())) != 0, leadrev);
    //we end the building of the contig sequence
    delete wcns;

}

//first test of the concurrent gapfilling of edges
void BFilling::create_edge_consensus(int number_threads) {
    pthread_mutex_t mut;
    pthread_mutex_t log;
    pthread_mutex_init(&mut, NULL);
    pthread_mutex_init(&log, NULL);
    //we copy an reverse the set of edges
    vector<uint64_t > cnstask;
    cnstask=eids;
    //we reverse the vector because we use pop->back() to get an edge
    reverse(cnstask.begin(),cnstask.end());
    concurrentedges tmp;
    tmp.mut=mut;//MUTEX to control access to files
    tmp.task=&cnstask;
    tmp.log=log;
    tmp.bf=this;
    //tmp.a=a;
    //we create the threads
    pthread_t *tab_threads= new pthread_t [number_threads];
    //create the threads
    for(int ii=0;ii<number_threads;ii++) {
        pthread_create(&tab_threads[ii], NULL, threaded_edges, &tmp);
    }
    //wait for thread to finish
    for(int ii=0;ii<number_threads;ii++)
    {
        pthread_join(tab_threads[ii], NULL);
    }
    //at this point all the edges were processed
}


void BFilling::create_edge_consensus() {
    vector<string> lseqs;
    vector<linking>  lin;
    int lid=0;
    string lseq="";
    //
    for(auto e:eids) {
        this->fillseqs2(e,20,lseqs,lin,lid,lseq);
        this->edge2consensus2(e,lseqs,lin,lid,lseq);
        this->print_cns_info(e);
        lseqs.clear();
        lseqs.shrink_to_fit();
        lin.clear();
        lin.shrink_to_fit();
        lseq.clear();
    }
}

//function that obtain the ends of contigs accoording to the given orientation in the line
void BFilling::get_edgeseq_ends(EdgeS* e,int d,string &seq_end_source, string &seq_end_target){

    uint32_t  source=e->getSource(),target=e->getTarget();

    if(edge_dir.count(pairf(source,target)) == 0){
        source=e->getTarget();
        target=e->getSource();
    }

    auto source_seq=linesG->get_seq_ctg_node(source);
    auto target_seq=linesG->get_seq_ctg_node(target);

    if(d > source_seq.length())
        d=source_seq.length();
    if(d > target_seq.length())
        d=target_seq.length();

    //string seq_o_source="";
    if(nodes2orientations[source] == 1){
        //+
        seq_end_source=source_seq.substr(source_seq.length()-d,d);
    }else{
        //-
        seq_end_source=revcomp(source_seq.substr(0,d));
    }
    //cout <<target<< " "<<nodes2orientations[target]<<" "<<target_seq.substr(0,abs(e->getBd()))<<endl;
    //string seq_o_target="";
    if(nodes2orientations[target] == 1){
        //+
        seq_end_target=target_seq.substr(0,d);
    }else{
        //-
        seq_end_target=revcomp(target_seq.substr(target_seq.length()-d,d));
    }

}

//return the best N long reads that pass by the
void BFilling::fillseqs2(uint64_t eid, int n_reads,vector<string> & lseqs,  vector<linking> & lin , int &lid, string &lseq){
   //links for each edge
    auto e=linesG->get_edge(eid);
    lid=e->get_leadseq_id(0);//we get the id of the lead seq

    if(lid == -1){
        printf("Error in getting leadseq for edge %d", e->getEdge_id());
        exit(1);
    }
    e->get_longread_links(n_reads,lin);//set the maximum of different long reads to make consensus return the 15 long reads whit more support on edges
    sort(lin.begin(),lin.end(),[](linking a,linking b){return a.p1<b.p1;});
    int sup=0;
    //int lid=0;
    //string lseq="";

    for(auto l:lin) {
        //should move this call to the calling method because we need to
        auto seq = longreadDB->ChunkGetSeq(l.lonread_id);
        //cout << l.c1<<" "<<l.c2<<" "<<l.lonread_id<<" "<<l.p1<<" "
        //    <<l.p2<<" "<<l.switched<<" "<<l.gaps<<" "<<l.gape<<" "<<l.ror1<<" "<<l.ror2<<" "<<seq.seq.substr(l.gaps,l.gape-l.gaps)<<endl;
        auto subseq = seq.seq.substr(l.gaps, l.gape - l.gaps);
        //I need to check the aligment strand for the long reads
        if (l.switched) {
            subseq = revcomp(subseq);
        }
        //we store the sequence in the local adaptor
        lseqs.push_back(subseq);
        if(sup < e->long_reads_pass(l.lonread_id)) {
            sup = e->long_reads_pass(l.lonread_id);
            lid=l.lonread_id;
            lseq=seq.seq;
        }
    }


}


void BFilling::build_lines() {
//Step 1: We create the lemon graph for running the matching without the dummy edges
    GraphS* g=linesG;
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
    //cout << "Number of CC "<< ncc << endl;
    //hash2vec nodes2cc;
    simplehash scc;//store the size of each cc
    hash2vec headers;

    //we iterate the nodes of the graph to add them to the cc and search for a header/tail of each scf
    for (ListGraph::NodeIt n(bp); n != INVALID; ++n) {
        int cnt = 0;
        for (ListGraph::IncEdgeIt e(bp, n); e != INVALID; ++e) {
            cnt++;
        }
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
            for (int i = 1; i <path.size() ; ++i) {
                p->add_node(lemonTog[path[i]]);
                nodes2lines[lemonTog[path[i]]]=current_path;
                //nodes2orientations[]
                /*if(circular_nodes.count(lemonTog[path[i]]))
                    p->set_as_circular(circular_nodes[lemonTog[path[i]]]);*/
            }
            //we set the orientation of the nodes in the line
            for (int i=0; i<path.size()-1; i++) {
                //cout << g->get_node(lemonTog[p])->getName() <<"->";
                auto u = g->get_node(lemonTog[path[i]]);
                auto v = g->get_node(lemonTog[path[i + 1]]);
                //is dummy edge
                if (u->getCtg_id() == v->getCtg_id()) {
                    //we got the contig from the node
                    auto contig = g->getContigs()->getcontig(u->getCtg_id());
                    int  strand = -1;
                    //means that we found a match
                    //if(v->getName().find("+")!=std::string::npos)
                    if (v->getId() == contig.header)
                        strand = 1;
                    //we set the contig orientations orientations
                    nodes2orientations[lemonTog[path[i]]]=strand;
                    nodes2orientations[lemonTog[path[i+1]]]=strand;
                    //we save the node two paths order
                    //nodes2pathpos[lemonTog[path[i]]]
                }else{//is an edge and we save the direction in were we are passing by it
                    edge_dir[pairf(lemonTog[path[i]],lemonTog[path[i+1]])]=1;
                }


            }
            //we print the path
            p->compute_path_length(g);
            //we save the path in the graph
            //cout << "Printing graph lines"<<endl;
            //p->print_path();
            linesP.push_back(p);
            current_path++;
        }
    }
    //we destroy the tmp lemon graph
    bp.clear();
}
