//
// Created by Alex Digenova on 3/5/19.
//

#include "SPolisher.h"


uint32_t SPolisher::_compute_edge_length(){
    uint32_t totalseq2polish=0;
    // STEP 1: we collect minimizers frequencies in a simple hash table
    for (int i = 0; i <g->get_number_edges() ; ++i) {
        auto e = g->get_edge(i);
        auto u = g->get_node(e->getSource());
        auto v = g->get_node(e->getTarget());
        //it is a mate-edge
        if(u->getCtg_id() !=v->getCtg_id()) {
            if(e->getEdge_overlap() >= 300 and ! e->isPolished()) {
                totalseq2polish += e->getEdge_overlap();
            }
        }
    }
    //
    return totalseq2polish;
}

void SPolisher::build_edge_index() {

    uint32_t totalseq2polish=0;
    //we set the mer-size and  the w size by default we work with (5,15)-minimizers
    sizeKmer=kmersize;
    kmer_type kmer, seed, seed_revcomp, minimizer;
    //we pick the minimizer every 5 bases
    int w=wsize;

    if (sizeKmer == (int)(sizeof(kmer_type)*4))
        kmerMask = -1;
    else
        kmerMask=(((kmer_type)1)<<(sizeKmer*2))-1;
    //to gather some stats
    uint32_t  number_edge_to_polish=0;
    uint32_t  total_length_edge_to_polish=0;
    //know indexing time
    clock_t begin = clock();
    // STEP 1: we collect minimizers frequencies in a simple hash table
    for (int i = 0; i <g->get_number_edges() ; ++i) {
        auto e = g->get_edge(i);
        auto u = g->get_node(e->getSource());
        auto v = g->get_node(e->getTarget());
        //it is a mate-edge
        if(u->getCtg_id() !=v->getCtg_id()) {
            //we save that this contigs are involved in edges and shold not be used for polishing
            ctginedges[u->getCtg_id()]++;
            ctginedges[v->getCtg_id()]++;
            //it should be longer than the minimal contig length and not marked as polished before
            if(e->getEdge_overlap() >= 300 and ! e->isPolished()) {
                //we have to index this sequence
                //this->edges2polish.push_back(e->getEdge());
                totalseq2polish+=e->getEdge_overlap();
                auto seq=e->getEdge_gapseq();
                number_edge_to_polish++;
                total_length_edge_to_polish+=e->getEdge_overlap();
                //we create a char* pointer to seq
                auto seqc=const_cast<char *>(seq.c_str());
                //we print the k-mers of this seq
                //we compute the canonical kmer of the sequence
                //uint64_t a=0,b=0;
                minimizer = extractKmerFromRead(seqc, 0, &seed, &seed_revcomp);
                int mp=0;
                bool ms= minimizer != seed;//false if fwd

                for (auto j=1; j<seq.length()-sizeKmer+1; j++) {
                    kmer = extractKmerFromRead(seqc, j, &seed, &seed_revcomp);
                    if(j%w == 0){
                        //printf("Selected in window=%d Minimizer=%d Mseq=%s\n",j/w-1,minimizer,seq.substr(mp,sizeKmer).c_str());
                        minimizers2counts[minimizer]++;
                        //we store the minimizer related information in an array
                        mmindex h;
                        h.pos= static_cast<uint32_t>(j);
                        h.mh=minimizer;
                        h.seid= static_cast<uint32_t>(i);
                        h.strand=ms;
                        //we store upto 100 locations for each minimizer, high copy number minimizers are removed
                        if(minimizers2counts[minimizer] <= maxfreq)
                            mpositions.push_back(h);

                        minimizer=kmer;
                        mp=j;
                        ms= minimizer != seed; //true if fwd
                    }else{
                        if(kmer < minimizer){
                            minimizer=kmer;
                            mp=j;
                            ms= minimizer != seed; //true if fwd
                        }
                    }
                }
                //log to know about the index construction
                if(number_edge_to_polish%10000==0){
                    printf("A total of %d edges have been indexed, total indexed seq=%d\n",number_edge_to_polish,totalseq2polish);
                }

            }
        }
    }
    //We compute an some minimizers stats
    simplehash2 hist;
    //simplehash2 kmer2index;
   // cout << "Total number of k-mers "<<kmers2counts.size()<<endl;
    //cout << " "<<minimizers2counts.size()<<endl;
    //cout << "Total number of (15,5)-minimizers positions "<<mpositions.size()<<endl;
    printf("A total of %d edges can be polish, the total length is %d\n",number_edge_to_polish,total_length_edge_to_polish);
    printf("Index information:\n");
    printf(" Total number of (%d,%d)-minimizers = %d\n",w,sizeKmer,(int)minimizers2counts.size());
    printf(" Total number of (%d,%d)-minimizers positions = %d\n",w,sizeKmer,(int)mpositions.size());

    uint32_t index=0;
    uint32_t maxcountmm=0;
    for ( auto it = minimizers2counts.begin(); it != minimizers2counts.end(); ++it ) {
          //cout << " " << it->first << ":" << it->second<<endl;
           hist[it->second]++;
        //we save the maximum freq observed
        if(maxcountmm < it->second)
            maxcountmm=it->second;
    }
    //sort(hist.begin(),hist.end(),[ ](const uint32_t))
    //we print the histogram
   /* for ( auto it = hist.begin(); it != hist.end(); ++it )
            cout << "HMM: " << it->first << " " << it->second<<endl;
    */
    //to not print more than 100bp
    if(maxcountmm > 100)
        maxcountmm=100;

    for(uint32_t k=0; k<maxcountmm; ++k) {
        if (hist.find(k) != hist.end())
            //cout << "HMM: " << k << " " << hist[k] << endl;
            printf("HMM: %d %d\n",k,hist[k]);
    }
    //Step 2, we sort the array of positions by minimizer,seq and pos
    sort( mpositions.begin( ), mpositions.end( ), [ ]( const mmindex a, const mmindex b)
    {
        return tie(a.mh,a.seid,a.pos) < tie(b.mh,b.seid,b.pos);
    });

    /*for(auto p:mpositions)
        printf("m=%d e=%d p=%d\n",p.mh,p.seid,p.pos);*/

    //we iterate the array of positions to store the minimizer start location in the array
    uint32_t rindex=0;
    kmer_type lmm=0;//last minimizer observed
    uint32_t  mm_excluded=0;
    for (auto &mmi : mpositions) {
        auto f=minimizers2counts[mmi.mh];
        if(mmi.mh != lmm) {
            if (f <= maxfreq) {
                minimizers2counts[mmi.mh] = rindex;
            }
            else{
                minimizers2counts[mmi.mh] = UINT32_MAX;
                mm_excluded++;
            }
            //we update the minimizer
            lmm= mmi.mh;
        }
      /*  if(minimizers2counts[mmi.mh] == UINT64_MAX)
            printf("D %d %d %d %d\n", mmi.mh,rindex,minimizers2counts[mmi.mh],f);
        else
            printf("K %d %d %d %d\n", mmi.mh,rindex,minimizers2counts[mmi.mh],f);*/
        rindex++;
    }
    //cout << "Total number of (5,15)-minimizers   "<<mm_excluded<<" F(k) > "<<maxfreq<<endl;
    printf(" Total number of (%d,%d)-minimizers excluded = %d\n",w,sizeKmer,mm_excluded);

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Time spent in minimizer index construction " << elapsed_secs <<" secs "<<endl;

    //end of minimizer indexing
}

SPolisher::SPolisher(GraphS* vmcg) {
    //max kmer size is 16 for the moment could be 32 but not necessary since we are just indexing edges
    kmersize=15;
    wsize=5;
    maxfreq=100;
    g=vmcg;
    //this code is for avoid expensive rehasing assume that each 15-mer is unique, only true for the vector of positions.
    auto s2p=_compute_edge_length();
    cout << "Total sequence length to polish = "<<s2p<<endl;
    minimizers2counts.reserve(s2p/wsize);
    cout << "Reserving  "<<s2p/wsize<<" buckets"<<endl;
    cout << "Bucket size = "<<minimizers2counts.bucket_count()<<endl;
    //we reserve the vector of positions
    mpositions.reserve(s2p/wsize);


}

SPolisher::SPolisher(int k, int w, int mf,GraphS* vmcg) {
    //max kmer size is 16 for the moment could be 32 but not necessary since we are just indexing edges
    kmersize=k;
    wsize=w;
    maxfreq=mf;
    g=vmcg;

    auto s2p=_compute_edge_length();
    cout << "Total sequence length to polish = "<<s2p<<endl;
    minimizers2counts.reserve(s2p/wsize);
    cout << "Reserving  "<<s2p/wsize<<" buckets"<<endl;
    cout << "Bucket size = "<<minimizers2counts.bucket_count()<<endl;
    //we reserve the vector of positions
    mpositions.reserve(s2p/wsize);
}




//we polish only the edges that were not polished in the Graph polishing step
void SPolisher::polish(unordered_map<int, int> & map) {

     //printing contigs already in edges
    /*for(auto i:ctginedges){
        printf("ctg=%d t=%d\n",i.first,i.second);
    }*/

    //determine which contigs are not in edges
    //STEP 1: we determine which contigs can we use for polishing the edges
    printf("Contigs not in edges:\n");
    vector<int> contigs2polish;
    uint32_t  clen=0;
    int cmin=1000000000;
    int cmax=0;
    for(auto c:g->getContigs()->get_all_contigs()){
        if(ctginedges.find(c.ctgid) == ctginedges.end()) {
            //printf("ctg=%d t=%d cov=%d\n", c.ctgid, c.length, c.coverage);
            contigs2polish.push_back(c.ctgid);
            clen+=c.length;
            if(c.length > cmax)
                cmax=c.length;
            if(c.length < cmin)
                cmin=c.length;
        }
    }

    printf("A total of %d contigs can be used for polishing edges\n",(int)contigs2polish.size());
    printf("Total contig length=%d  min=%d max=%d\n",clen,cmin,cmax);
    printf("A total of %d contigs were used in the Graph Polishing Step:\n",(int)map.size());
    // q: should we include contigs used in the GPolishing step?
    // A: yes, we have to use such contigs, we excluded the edges that were already polished in the GPolished step, to not over polish
    /*for(auto i:map){
        printf("ctg=%d t=%d\n",i.first,i.second);
    }*/

    //STEP 2: we map the contigs minimizer to the edge to determine potential contig edge pairs for polishing
    //we collect the hits of each contig
    clock_t begin = clock();
    vector<mmhit> hitsdb;//database of hits between edges and contigs
    uint32_t numberofctgs=0;
    for (auto cid:contigs2polish) {
        //auto c=g->getContigs()->getcontig(cid);
        //printf("Mapping contig %d:\n",cid);
        //we store upto a maximum of the copy number

        _map_contig2edges(cid, hitsdb);
        if(numberofctgs%10000==0){
            cout << "A total of "<<numberofctgs<<" contigs Have been mappad,  "<<hitsdb.size()<<" hits have been collected\n";
        }
        numberofctgs++;
    }
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Time spent in mapping contigs 2 edges :" << elapsed_secs <<" secs "<<endl;
    printf("A total of %d hits were collected\n",(int)hitsdb.size());
    //STEP 3: we sort the hits collected by edge, position, number of minimizers and approximate matching length
    sort( hitsdb.begin( ), hitsdb.end( ), []( const mmhit a, const mmhit b)
    {
        return tie(a.tid,a.rs,a.cnt,a.mlen) > tie(b.tid,b.rs,b.cnt,b.mlen);
    });

    /*//we print the hits
    for(auto h:hitsdb)
        printf("HIT: cid=%d eid=%d strand=%d qs=%d qe=%d rs=%d re=%d cnt=%d mlen=%d blen=%d\n",
               h.qid, h.tid, h.strand, h.qs, h.qe, h.rs, h.re, h.cnt, h.mlen, h.blen);
*/

    begin = clock();
    //we have to select the best cover for each edge (greedy algorithm)
    if(!hitsdb.empty())
        _hits2polish_edge(hitsdb);

    //we end the polishing by contig alignments
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Time spent in polishing edges :" << elapsed_secs <<" secs "<<endl;


}

void SPolisher::_map_contig2edges(int cid, vector<mmhit> &hitsdb) {
    //1-we get the contig sequence
    auto c=g->getContigs()->getcontig(cid);
    auto seqc=const_cast<char *>(c.seq.c_str());
    kmer_type kmer, seed, seed_revcomp, minimizer;
    int w=wsize;
    //we pick the minimizers for fill the vector of hits
    minimizer = extractKmerFromRead(seqc, 0, &seed, &seed_revcomp);
    auto cstrand=minimizer != seed;

    uint32_t mp=0;
    vector<mmhsp>hpos_fwd;//store + kmers
    vector<mmhsp>hpos_rev;//store - kmers

    for (auto j=1; j<c.seq.length()-sizeKmer+1; j++) {
        kmer = extractKmerFromRead(seqc, j, &seed, &seed_revcomp);
        if(j%w == 0){
            //we match the minimizer to the database
            if(minimizers2counts.find(minimizer) != minimizers2counts.end()) {
                //we skypt masked minimizers
                if(minimizers2counts[minimizer]==UINT32_MAX)
                    continue;
                    //we print the hits of the minimizer
                for (int i = minimizers2counts[minimizer]; i <mpositions.size() ; ++i) {
                        if(mpositions[i].mh == minimizer) {
                          //  printf("Matching minimizer in window=%d Minimizer=%d cid=%d eid=%d pose=%d posc=%d\n",
                            //        j / w, minimizer, cid, mpositions[i].seid, mpositions[i].pos, mp);
                            //we create an mmhsp
                            mmhsp h;
                            h.pos=mpositions[i].pos;
                            h.seid=mpositions[i].seid;
                            h.cpos=mp;
                            h.cid= static_cast<uint32_t>(cid);
                            //it foward if the strand is stored in the
                            h.strand=cstrand != mpositions[i].strand; //fwd if both were obtained from the same strad
                            //we store the hsp according to the strand
                            if(h.strand)
                                hpos_rev.push_back(h);
                            else
                                hpos_fwd.push_back(h);
                        }else{
                            break;
                        }
                }
            }
            //printf("Selected in window=%d Minimizer=%d Mseq=%s\n",j/w-1,minimizer,seq.substr(mp,sizeKmer).c_str());
            //minimizers2counts[minimizer]++;
            //we update the minimizer related information pos, kmer, cstrand
            minimizer=kmer;
            mp=j;
            cstrand=minimizer != seed;
        }else{
            //we ask if we have to change the current minimizer in the window
            if(kmer < minimizer){
                minimizer=kmer;
                mp=j;
                cstrand=minimizer != seed;
            }
        }
    }

    //we collect the hits of each contig
    vector<mmhit> localhits;
    //we collect the hits
    hsps2hits(hpos_fwd,localhits);
    hsps2hits(hpos_rev,localhits);

    /*//we sort the hits according to a minimizers and the coverage of the contig
    sort( localhits.begin( ), localhits.end( ), []( const mmhit a, const mmhit b)
    {
        return tie(a.cnt,a.mlen) > tie(b.cnt,b.mlen);
    });*/
    //we have the candidates edges for each contig
    for(auto h:localhits)
        //50% of contig aligned to the edge sequence
        if( (float) h.mlen/c.length > 0.5) {
            /*printf("HIT: cid=%d eid=%d strand=%d qs=%d qe=%d rs=%d re=%d cnt=%d mlen=%d blen=%d clen=%d ccov=%d rep=%d\n",
                   h.qid, h.tid, h.strand, h.qs, h.qe, h.rs, h.re, h.cnt, h.mlen, h.blen, c.length, c.coverage,
                   c.repeat);*/
            //we save the hit to the hitdb database
            hitsdb.push_back(h);
        }
    //we clean for the moment the local container
    localhits.erase(localhits.begin(),localhits.end());
}





void SPolisher::best_collinear2hits(vector<mmhsp> & hsps, vector<mmhit> & localhits, vector<int> &lis){
    //start of hit
    auto qstart=0;
    int max_g=1000;//maximum gap in a chain
    uint32_t cnt=1;
    uint32_t blen=0;
    uint32_t mlen=0;
    auto f=hsps[lis[0]];//we pick the first hsp
    //todo: compute the number of different minimizers used in the collinear blocks
    //we create a hit only for appropiate collinear blocks
    int i=1;
    for(i=1; i< lis.size(); ++i){
        //distance in edge
        int tl=abs((int)(hsps[lis[i-1]].pos-hsps[lis[i]].pos));
        //distance in query
        int ql=abs((int)(hsps[lis[i-1]].cpos-hsps[lis[i]].cpos));
        //we have to broke the line at i
        if(tl > max_g or ql > max_g){
            //we have to break the hit
            if(cnt >= 4){
                //we have to create the hit
                mmhit hit;
                hit.cnt=cnt;
                hit.mlen=mlen;
                hit.blen=blen;
                //todo: I have to check this there is some cases wehre qs is larger than qe
                hit.rs= f.strand ? hsps[lis[i-1]].pos:hsps[lis[qstart]].pos;
                hit.qs= f.strand ? hsps[lis[i-1]].cpos:hsps[lis[qstart]].cpos;
                hit.re= f.strand ? hsps[lis[qstart]].pos:hsps[lis[i-1]].pos;
                hit.qe= f.strand ? hsps[lis[qstart]].cpos:hsps[lis[i-1]].cpos;
                hit.qid = f.cid;
                hit.tid = f.seid;
                hit.strand = f.strand;//save the strand of the contig
                //we accomodate the reference sequence
                if(f.strand){
                    //we have to switch the ref strand
                    auto tmp=hit.rs;
                    hit.rs=hit.re;
                    hit.re=tmp;
                }
                //we kept only hits larger than 200
                if(hit.qs<hit.qe and hit.rs<hit.re and hit.qe-hit.qs >=200 and hit.re-hit.rs >=200)
                    localhits.push_back(hit);
            }

            cnt=1;
            qstart=i;//we update the start to i
            blen=0;
            mlen=0;

        }else{
            cnt++;
            blen += tl > ql ? tl : ql;
            mlen += ql < tl ? ql : tl;
        }
    }
    //we have to rescue the hit if the chain was not broken
    if(cnt >= 4){
        mmhit hit;
        hit.cnt=cnt;
        hit.mlen=mlen;
        hit.blen=blen;
        hit.rs= f.strand ? hsps[lis[i-1]].pos:hsps[lis[qstart]].pos;
        hit.qs= f.strand ? hsps[lis[i-1]].cpos:hsps[lis[qstart]].cpos;
        hit.re= f.strand ? hsps[lis[qstart]].pos:hsps[lis[i-1]].pos;
        hit.qe= f.strand ? hsps[lis[qstart]].cpos:hsps[lis[i-1]].cpos;
        hit.qid = f.cid;
        hit.tid = f.seid;
        hit.strand = f.strand;//save the strand of the contig
        if(f.strand){
            //we have to switch the ref strand
            auto tmp=hit.rs;
            hit.rs=hit.re;
            hit.re=tmp;
        }
        //we kept only hits larger than 200
        if(hit.qs<hit.qe and hit.rs<hit.re and hit.qe-hit.qs >=200 and hit.re-hit.rs >=200)
        localhits.push_back(hit);
    }
}


void SPolisher::hsps2hits(vector<mmhsp> & hsps, vector<mmhit> & localhits){

    if(hsps.size() < 4) // 4 * 15 at least 60bp alignemnt
        return;
    //we sort the hsps by edge and pos
    sort( hsps.begin( ), hsps.end( ), []( const mmhsp a, const mmhsp b)
    {
        return tie(a.seid,a.pos) < tie(b.seid,b.pos);
    });

    //vector<mmhit> localhits;

    //we search the longest chain by
    auto f=hsps[0];//first hsp
    vector<int> a;//array
    vector<int> lis;//list of index
    //we iter each of the potential edges
    int offset=0;
    for (int k=0; k<hsps.size(); ++k) {
        auto h=hsps[k];
        if(h.seid == f.seid) {
            a.push_back(h.cpos);
        }
        else{
            //we have to compute the largest collinear block
            if(f.strand)
                reverse(a.begin(),a.end());
            //we determine the largest collinear block
            find_largest_collinear_block(a,lis,offset);
            if(lis.size() >= 4 ){
                best_collinear2hits(hsps,localhits,lis);
            }
            //we have to erase the a array and lis array
            a.erase(a.begin(),a.end());
            lis.erase(lis.begin(),lis.end());
            f=h;
            a.push_back(f.cpos);
            offset=k;
        }
        //printf("HSP: c=%d e=%d cp=%d ep=%d strand=%d offset=%d\n", h.cid, h.seid, h.cpos, h.pos, h.strand, offset);
    }

    //we have more than 1 element
    if(a.size() >= 4){
        if(f.strand)
            reverse(a.begin(),a.end());

        find_largest_collinear_block(a,lis,offset);
        //we resue the last hits
        if(lis.size() >= 4 ){
            best_collinear2hits(hsps,localhits,lis);
        }
        //we have to erase the a
        a.erase(a.begin(),a.end());
        lis.erase(lis.begin(),lis.end());
    }

    /*//we print the hits if there are hits
    for(auto h:localhits)
        printf("HIT: cid=%d eid=%d strand=%d qs=%d qe=%d rs=%d re=%d cnt=%d mlen=%d blen=%d\n"
                ,h.qid,h.tid,h.strand,h.qs,h.qe,h.rs,h.re,h.cnt,h.mlen,h.blen);
*/
    //we clean for the moment the local container
    //localhits.erase(localhits.begin(),localhits.end());

}

/* Finds longest strictly increasing subsequence. O(n log k) algorithm. */
void SPolisher::find_largest_collinear_block(vector<int> &a, vector<int> &b, int offset){
    vector<int> p(a.size());
    int u, v;

    if (a.empty()) return;

    b.push_back(0);

    for (size_t i = 1; i < a.size(); i++)
    {
        // If next element a[i] is greater than last element of
        // current longest subsequence a[b.back()], just push it at back of "b" and continue
        if (a[b.back()] < a[i])
        {
            p[i] = b.back();
            b.push_back(i);
            continue;
        }

        // Binary search to find the smallest element referenced by b which is just bigger than a[i]
        // Note : Binary search is performed on b (and not a).
        // Size of b is always <=k and hence contributes O(log k) to complexity.
        for (u = 0, v = b.size()-1; u < v;)
        {
            int c = (u + v) / 2;
            if (a[b[c]] < a[i]) u=c+1; else v=c;
        }

        // Update b if new value is smaller then previously referenced value
        if (a[i] < a[b[u]])
        {
            if (u > 0) p[i] = b[u-1];
            b[u] = i;
        }
    }

    for (u = b.size(), v = b.back(); u--; v = p[v]) b[u] = v;

    //we have to update the offset of the array
    if(offset > 0){
        for (int &j : b) {
            j +=offset;
        }
    }
}


//we polish each each edge with the selected hits
void SPolisher::_hits2polish_edge(vector<mmhit> & hitsdb) {

    if(hitsdb.empty())
        return;

    //we collect the hits for each each
    auto f=hitsdb[0];
    vector<mmhit> edgehits;
    for (auto h:hitsdb) {
        if(h.tid==f.tid)
            edgehits.push_back(h);
        else{
            //we have to polish the current edge
            _polish_edge(edgehits);
            //we update the current first hit
            f=h;
            edgehits.clear();
            //we store the current
            edgehits.push_back(h);
        }
    }

    //we polish the last edge if necessary
    if(!edgehits.empty())
        _polish_edge(edgehits);

}

void SPolisher::_polish_edge(vector<mmhit> & edgehits) {


    //we print the hits
    printf("Hits at edge level:\n");

    //this method expect an int has parameter and refers to the edge position in the graph
    auto e=g->get_edge((int)edgehits[0].tid);
    auto e_seq=e->getEdge_gapseq();

    CJoiner aln;//this class will became my aligment engine
    //identity is computing without taking into account the SNP, should we s
    float min_iden=0.8;//base identity
    if(e->getEdge_avg_cns() >=2)
        min_iden=0.85;
    if(e->getEdge_avg_cns() >=3)
        min_iden=0.91;
    if(e->getEdge_avg_cns() >=4)
        min_iden=0.93;
    //we have to be more strict because is only based on aligments
    if(e->getEdge_avg_cns() >=5)
        min_iden=0.94;
    if(e->getEdge_avg_cns() >=10)
        min_iden=0.96;
    if(e->getEdge_avg_cns() >=15)
        min_iden=0.97;
    if(e->getEdge_avg_cns() >=20)
        min_iden=0.985;

    // I have to select the best hits based on minimizer matches
    bool polished=false;
    //we select the best hits based on a simple greedy algorithm
    vector<mmhit> bestctgs;
    _get_best_ctg_for_edge(edgehits,bestctgs);

    printf("A total of %d ctgs were selected for polishing from %d canditates\n",(int)bestctgs.size(),(int)edgehits.size());
    //we do the aligments to polish the edge sequence
    for(auto h:bestctgs) {
        printf("HIT: cid=%d eid=%d strand=%d qs=%d qe=%d rs=%d re=%d cnt=%d mlen=%d blen=%d min_iden=%f\n",
               h.qid, h.tid, h.strand, h.qs, h.qe, h.rs, h.re, h.cnt, h.mlen, h.blen, min_iden);

        auto c=g->getContigs()->getcontig(h.qid);
        //g->getContigs()->set_used_in_polishing()
        c.used_in_polishing=true;//we set that the conntig was used in polishing
        //we trim the contig sequence to the estimated by the minimizers
        //auto c_seq=c.seq.substr(h.qs,(h.qe-h.qs)+1);
        auto c_seq= h.strand == 0 ? c.seq.substr(h.qs,(h.qe-h.qs)+1) : revcomp(c.seq.substr(h.qs,(h.qe-h.qs)+1)); //we add the contig in the correct orientation
        //we revcomp the contig seq if the match was in the reverse order
        auto match=aln.compute_optimal_aligment_trim(c_seq,e_seq,false,e->getEdge_avg_cns());
        //we polish only if min identity and only if 300 bases and (h.qe-h.qs >= 300)
        if(match.qcov >=75 and match.iden >= min_iden){
            //we can polish the edge sequence
            //we mask the target sequence to find another aligment
            //e_seq.replace(hit.tstart,hit.tstop,abs(hit.tstop-hit.tstart),'N');
            //we extract the portion of the query sequence involved in polishing
            auto p_seq=c_seq.substr(match.qstart,abs(match.qstop-match.qstart)+1);//the illumina seq
            //we replace the sequence of the target seq by the illumina sequence
            e_seq.replace(match.tstart,abs(match.tstart-match.tstop)+1,p_seq);
            polished=true;
            //we save that the contig was used in polishing
            //this->ctgused2polish[contig.ctgid]++;
        }
    }

    //we check if the edge was polished and we update it;
    if(polished){
        //we save the gapseq in the strand of the edge, so we then revert when we traverse the line
        e->set_edge_gapseq(e_seq.c_str());
        e->setEdge_overlap(e_seq.length());
        e->setPolished(true);
    }


}

void SPolisher::_get_best_ctg_for_edge(vector<mmhit> &edgehits, vector<mmhit> &bestctgs) {

    //we sort the contig by the number of exact match
    sort( edgehits.begin( ), edgehits.end( ), []( const mmhit a, const mmhit b)
    {
        return tie(a.cnt,a.mlen) > tie(b.cnt,b.mlen);
    });

    for (auto h:edgehits) {
         if(bestctgs.empty())
             bestctgs.push_back(h);
         else
             //check if the contig that we want to add dont  have overlap with a previous added contig
             if(!_compute_overlap(h,bestctgs))
                 bestctgs.push_back(h);
                    //we can add only if dont a significant overlap with a previous added contig
    }

}

bool SPolisher::_compute_overlap(mmhit &h, vector<mmhit> &bestctgs) {
    int overlap=0;
    auto add=true;
    for(auto b:bestctgs){
        //we hit is contained in a previous best hit
        if(h.rs >= b.rs and h.re <= b.re )
            return 1;
        //the hit has an overlap at 3' with a previous added contig
        if(h.rs <= b.rs and h.re <= b.re and h.re >= b.rs)
                overlap+=(h.re-b.rs);
        //the hit has a overlap at 5' with a previous added contig
        if(h.rs >= b.rs and h.rs <= b.re and h.re >= b.re)
            overlap+=(b.re-h.rs);
        //the hit covers entierly a previous one
        if(h.rs < b.rs and h.re > b.re and h.rs < b.re)
            return 1;
    }
    //we compute the alignment length
    int not_covered=abs((int)h.re-(int)h.rs)-overlap;
    //the contig add 100bp to the polising
    return not_covered < 200;
}

//constructor for testing the mapping of ctg vs edges
SPolisher::SPolisher(int k, int w, int mf, string fctg, string fref){

    //auto s2p=_compute_edge_length();
    kmersize=k;
    wsize=w;
    maxfreq=mf;
    //we read and load the fasta file of edges
    if(k > 32){
        cout << "Max minimizer size is 32"<<endl;
        exit(1);
    }


    this->read_fasta(fctg,queries);
    //we store the contigs index
    uint32_t i=0;
    for(auto c:queries){
        q2pos[c.ctgid]=i;
        i++;
    }
    this->read_fasta(fref,targets);
    uint32_t j=0;
    for(auto c:targets){
        t2pos[c.ctgid]=j;
        j++;
    }
    cout << "Number of queries:"<<queries.size()<<endl;
    cout << "Number of targets:"<<targets.size()<<endl;
    uint64_t s2p=0;
    for(auto c:queries)
            s2p+=c.length;
    cout << "Total sequence length to polish = "<<s2p<<endl;
    minimizers2counts.reserve(s2p/wsize);
    cout << "Reserving  "<<s2p/wsize<<" buckets"<<endl;
    cout << "Bucket size = "<<minimizers2counts.bucket_count()<<endl;
    //we reserve the vector of positions
    mpositions.reserve(s2p/wsize);

}

//Build edge index
//index is working
void SPolisher::build_index(){
    uint32_t totalseq2polish=0;
    //we set the mer-size and  the w size by default we work with (5,15)-minimizers
    sizeKmer=kmersize;
    kmer_type kmer, seed, seed_revcomp, minimizer;
    //we pick the minimizer every 5 bases
    int w=wsize;

    if (sizeKmer == (int)(sizeof(kmer_type)*4))
        kmerMask = -1;
    else
        kmerMask=(((kmer_type)1)<<(sizeKmer*2))-1;
    //to gather some stats
    uint32_t  number_edge_to_polish=0;
    uint32_t  total_length_edge_to_polish=0;
    //know indexing time
    clock_t begin = clock();
    //fisrt pass inndexinng the k-mers
    for(auto c:targets){
        //we skypt short edges
        if(c.length < 300)
             continue;
        //we build the index
        totalseq2polish+=c.length;
        auto seq=c.seq;
        number_edge_to_polish++;
        total_length_edge_to_polish+=c.length;
        //we create a char* pointer to seq
        auto seqc=const_cast<char *>(seq.c_str());
        //we print the k-mers of this seq
        //we compute the canonical kmer of the sequence
        //uint64_t a=0,b=0;
        minimizer = extractKmerFromRead(seqc, 0, &seed, &seed_revcomp);
        int mp=0;
        bool ms= minimizer != seed;//false if fwd

        for (auto j=1; j<seq.length()-sizeKmer+1; j++) {
            kmer = extractKmerFromRead(seqc, j, &seed, &seed_revcomp);
            if(j%w == 0){
                //printf("Selected in window=%d Minimizer=%d Mseq=%s\n",j/w-1,minimizer,seq.substr(mp,sizeKmer).c_str());
                minimizers2counts[minimizer]++;
                //we store the minimizer related information in an array
                mmindex h;
                h.pos= static_cast<uint32_t>(j);
                h.mh=minimizer;
                h.seid= static_cast<uint32_t>(c.ctgid);
                h.strand=ms;
                //we store upto 100 locations for each minimizer, high copy number minimizers are removed
                if(minimizers2counts[minimizer] <= maxfreq)
                    mpositions.push_back(h);

                minimizer=kmer;
                mp=j;
                ms= minimizer != seed; //true if fwd
            }else{
                if(kmer < minimizer){
                    minimizer=kmer;
                    mp=j;
                    ms= minimizer != seed; //true if fwd
                }
            }
        }

        if(number_edge_to_polish%10000==0){
            printf("A total of %d edges have been indexed, total indexed seq=%d maxfreq=%d\n",number_edge_to_polish,totalseq2polish,maxfreq);
        }

    }

    //We compute an some minimizers stats
    simplehash2 hist;
    printf("A total of %d edges can be polish, the total length is %d\n",number_edge_to_polish,total_length_edge_to_polish);
    printf("Index information:\n");
    printf(" Total number of (%d,%d)-minimizers = %d\n",w,sizeKmer,(int)minimizers2counts.size());
    printf(" Total number of (%d,%d)-minimizers positions = %d\n",w,sizeKmer,(int)mpositions.size());

    uint32_t index=0;
    uint32_t maxcountmm=0;
    for ( auto it = minimizers2counts.begin(); it != minimizers2counts.end(); ++it ) {
        //cout << " " << it->first << ":" << it->second<<endl;
        hist[it->second]++;
        //we save the maximum freq observed
        if(maxcountmm < it->second)
            maxcountmm=it->second;
    }
    //to not print more than 100bp
    if(maxcountmm > 100)
        maxcountmm=100;

    for(uint32_t k=0; k<maxcountmm; ++k) {
        if (hist.find(k) != hist.end())
            //cout << "HMM: " << k << " " << hist[k] << endl;
            printf("HMM: %d %d\n",k,hist[k]);
    }
    //Step 2, we sort the array of positions by minimizer,seq and pos
    sort( mpositions.begin( ), mpositions.end( ), [ ]( const mmindex a, const mmindex b)
    {
        return tie(a.mh,a.seid,a.pos) < tie(b.mh,b.seid,b.pos);
    });

    /*for(auto p:mpositions)
        printf("m=%d e=%d p=%d\n",p.mh,p.seid,p.pos);*/

    //we iterate the array of positions to store the minimizer start location in the array
    uint32_t rindex=0;
    kmer_type lmm=0;//last minimizer observed
    uint32_t  mm_excluded=0;
    for (auto &mmi : mpositions) {
        auto f=minimizers2counts[mmi.mh];
        if(mmi.mh != lmm) {
            if (f <= maxfreq) {
                minimizers2counts[mmi.mh] = rindex;
            }
            else{
                minimizers2counts[mmi.mh] = UINT32_MAX;
                mm_excluded++;
            }
            //we update the minimizer
            lmm= mmi.mh;
        }
        rindex++;
    }
    //cout << "Total number of (5,15)-minimizers   "<<mm_excluded<<" F(k) > "<<maxfreq<<endl;
    printf(" Total number of (%d,%d)-minimizers excluded = %d\n",w,sizeKmer,mm_excluded);

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Time spent in minimizer index construction " << elapsed_secs <<" secs "<<endl;

}

//now the mapping function that we want to evaluate
void SPolisher::map_contigs(){
    //we map the contigs to the edges
    clock_t begin = clock();
    vector<mmhit> hitsdb;//database of hits between edges and contigs
    uint32_t numberofctgs=0;
    for (auto c:queries) {
        //auto c=g->getContigs()->getcontig(cid);
        //printf("Mapping contig %d:\n",cid);
        //we store upto a maximum of the copy number
        //_map_contig2edges(cid, hitsdb);
        //cout << "Mapping contig "<<c.namectg<<" length:"<<c.length<<endl;
        _map_contig2edges(c,hitsdb);
        if(numberofctgs%10000==0){
            cout << "A total of "<<numberofctgs<<" contigs Have been mappad,  "<<hitsdb.size()<<" hits have been collected\n";
        }
        numberofctgs++;
    }
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Time spent in mapping contigs 2 edges :" << elapsed_secs <<" secs "<<endl;
    printf("A total of %d hits were collected\n",(int)hitsdb.size());

    //we polish the edges
    //STEP 3: we sort the hits collected by edge, position, number of minimizers and approximate matching length
    sort( hitsdb.begin( ), hitsdb.end( ), []( const mmhit a, const mmhit b)
    {
        return tie(a.tid,a.rs,a.cnt,a.mlen) > tie(b.tid,b.rs,b.cnt,b.mlen);
    });

    begin = clock();
    //we have to select the best cover for each edge (greedy algorithm)
    if(!hitsdb.empty())
        _hits2polish_edge2(hitsdb);

    //we end the polishing by contig alignments
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Time spent in polishing edges :" << elapsed_secs <<" secs "<<endl;





}

//we polish each each edge with the selected hits
void SPolisher::_hits2polish_edge2(vector<mmhit> & hitsdb) {

    if(hitsdb.empty())
        return;

    //we collect the hits for each each
    auto f=hitsdb[0];
    vector<mmhit> edgehits;
    for (auto h:hitsdb) {
        if(h.tid==f.tid)
            edgehits.push_back(h);
        else{
            //we have to polish the current edge
            _polish_edge2(edgehits);
            //we update the current first hit
            f=h;
            edgehits.clear();
            //we store the current
            edgehits.push_back(h);
        }
    }

    //we polish the last edge if necessary
    if(!edgehits.empty())
        _polish_edge2(edgehits);

}

void SPolisher::_polish_edge2(vector<mmhit> & edgehits) {


    //we print the hits
    printf("Hits at edge level:\n");

    //this method expect an int has parameter and refers to the edge position in the graph
    //auto e=g->get_edge((int)edgehits[0].tid);
    //auto e_seq=e->getEdge_gapseq();
    assert(edgehits[0].tid == targets[t2pos[edgehits[0].tid]].ctgid);
    auto e_seq=targets[t2pos[edgehits[0].tid]].seq;

    CJoiner aln;//this class will became my aligment engine
    //identity is computing without taking into account the SNP, should we s
    float min_iden=0.8;//base identity
    //we have to be more strict because is only based on aligments
    /*if(e->getEdge_avg_cns() >=5)
        min_iden=0.9;
    if(e->getEdge_avg_cns() >=10)
        min_iden=0.95;
    if(e->getEdge_avg_cns() >=15)
        min_iden=0.97;
    if(e->getEdge_avg_cns() >=20)
        min_iden=0.98;
*/
    // I have to select the best hits based on minimizer matches
    bool polished=false;
    //we select the best hits based on a simple greedy algorithm
    vector<mmhit> bestctgs;
    _get_best_ctg_for_edge(edgehits,bestctgs);

    printf("A total of %d ctgs were selected for polishing from %d canditates\n",(int)bestctgs.size(),(int)edgehits.size());
    //we do the aligments to polish the edge sequence
    for(auto h:bestctgs) {
        printf("HIT: cid=%d eid=%d strand=%d qs=%d qe=%d rs=%d re=%d cnt=%d mlen=%d blen=%d min_iden=%f\n",
               h.qid, h.tid, h.strand, h.qs, h.qe, h.rs, h.re, h.cnt, h.mlen, h.blen, min_iden);

        //auto c=g->getContigs()->getcontig(h.qid);
        auto c=queries[q2pos[h.qid]];
        //we trim the contig sequence to the estimated by the minimizers
        //auto c_seq=c.seq.substr(h.qs,(h.qe-h.qs)+1);
        auto c_seq= h.strand == 0 ? c.seq.substr(h.qs,(h.qe-h.qs)+1) : revcomp(c.seq.substr(h.qs,(h.qe-h.qs)+1)); //we add the contig in the correct orientation
        //we skypt short
      /*  if(h.qs > h.qe){
            cout << "Warning:: h.qs > h.qe"<<endl;
        }
        //todo: this is a provisory solution, I have to check why there is hits having the previous warning
        if(c_seq.length() < 100 || e_seq.length() < 100)
            continue;*/
        //we revcomp the contig seq if the match was in the reverse order
        //auto match=aln.compute_optimal_aligment_trim(c_seq,e_seq,false,e->getEdge_avg_cns());
        auto match=aln.compute_optimal_aligment_trim(c_seq,e_seq,false,3);
        //we polish only if min identity and only if 300 bases and (h.qe-h.qs >= 300)
        if(match.qcov >=75 and match.iden >= min_iden){
            //we can polish the edge sequence
            //we mask the target sequence to find another aligment
            //e_seq.replace(hit.tstart,hit.tstop,abs(hit.tstop-hit.tstart),'N');
            //we extract the portion of the query sequence involved in polishing
            auto p_seq=c_seq.substr(match.qstart,abs(match.qstop-match.qstart)+1);//the illumina seq
            //we replace the sequence of the target seq by the illumina sequence
            e_seq.replace(match.tstart,abs(match.tstart-match.tstop)+1,p_seq);
            polished=true;
            //we save that the contig was used in polishing
            //this->ctgused2polish[contig.ctgid]++;
        }
    }

    //we check if the edge was polished and we update it;
    if(polished){
        //we save the gapseq in the strand of the edge, so we then revert when we traverse the line
        /*e->set_edge_gapseq(e_seq.c_str());
        e->setEdge_overlap(e_seq.length());
        e->setPolished(true);*/
    }


}



void SPolisher::_map_contig2edges(ctg &c, vector<mmhit> &hitsdb) {
    //1-we get the contig sequence
    //auto c=cid;
    auto seqc=const_cast<char *>(c.seq.c_str());
    if(c.seq.length() <100)
        return;
    //we set the mer-size and  the w size by default we work with (5,15)-minimizers
    sizeKmer=kmersize;
    kmer_type kmer, seed, seed_revcomp, minimizer;
    //we pick the minimizer every 5 bases
    int w=wsize;

    if (sizeKmer == (int)(sizeof(kmer_type)*4))
        kmerMask = -1;
    else
        kmerMask=(((kmer_type)1)<<(sizeKmer*2))-1;


    /*kmer_type kmer, seed, seed_revcomp, minimizer;
    int w=wsize;*/
    //we pick the minimizers for fill the vector of hits
    minimizer = extractKmerFromRead(seqc, 0, &seed, &seed_revcomp);
    auto cstrand=minimizer != seed;

    uint32_t mp=0;
    vector<mmhsp>hpos_fwd;//store + kmers
    vector<mmhsp>hpos_rev;//store - kmers

    for (auto j=1; j<c.seq.length()-sizeKmer+1; j++) {
        kmer = extractKmerFromRead(seqc, j, &seed, &seed_revcomp);
        if(j%w == 0){
            //we match the minimizer to the database
            if(minimizers2counts.find(minimizer) != minimizers2counts.end()) {
                //we skypt masked minimizers
                if(minimizers2counts[minimizer]==UINT32_MAX)
                    continue;
                //we print the hits of the minimizer
                for (int i = minimizers2counts[minimizer]; i <mpositions.size() ; i++) {
                    if(mpositions[i].mh == minimizer) {
                        /* printf("Matching minimizer in window=%d Minimizer=%d cid=%d eid=%d pose=%d posc=%d cnanme=%s\n",
                               j / w, minimizer, c.ctgid, mpositions[i].seid, mpositions[i].pos, mp,c.namectg.c_str());*/
                        //we create an mmhsp
                        mmhsp h;
                        h.pos=mpositions[i].pos;
                        h.seid=mpositions[i].seid;
                        h.cpos=mp;
                        h.cid= static_cast<uint32_t>(c.ctgid);
                       // h.m=minimizer;//we save the minimizer
                        //it foward if the strand is stored in the
                        h.strand=cstrand != mpositions[i].strand; //fwd if both were obtained from the same strad
                        //we store the hsp according to the strand
                        if(h.strand)
                            hpos_rev.push_back(h);
                        else
                            hpos_fwd.push_back(h);
                    }else{
                        break;
                    }
                }
            }
            //printf("Selected in window=%d Minimizer=%d Mseq=%s\n",j/w-1,minimizer,c.seq.substr(mp,sizeKmer).c_str());
            //minimizers2counts[minimizer]++;
            //we update the minimizer related information pos, kmer, cstrand
            minimizer=kmer;
            mp=j;
            cstrand=minimizer != seed;
        }else{
            //we ask if we have to change the current minimizer in the window
            if(kmer < minimizer){
                minimizer=kmer;
                mp=j;
                cstrand=minimizer != seed;
            }
        }
    }
    //we check if there is minimizers problems
   /* cout << "Hits F:"<<hpos_fwd.size()<<" R:"<<hpos_rev.size()<<endl;
    hpos_fwd.clear();
    hpos_rev.clear();
    return;*/
    //we collect the hits of each contig
    vector<mmhit> localhits;
    //we collect the hits
    hsps2hits(hpos_fwd,localhits);
    hsps2hits(hpos_rev,localhits);

    /*//we sort the hits according to a minimizers and the coverage of the contig
    sort( localhits.begin( ), localhits.end( ), []( const mmhit a, const mmhit b)
    {
        return tie(a.cnt,a.mlen) > tie(b.cnt,b.mlen);
    });*/
    //we have the candidates edges for each contig
    for(auto h:localhits)
        //50% of contig aligned to the edge sequence
        if( (float) h.mlen/c.length > 0.5) {
            /*printf("HIT: cid=%d eid=%d strand=%d qs=%d qe=%d rs=%d re=%d cnt=%d mlen=%d blen=%d clen=%d ccov=%d rep=%d\n",
                   h.qid, h.tid, h.strand, h.qs, h.qe, h.rs, h.re, h.cnt, h.mlen, h.blen, c.length, c.coverage,
                   c.repeat);*/
            //we save the hit to the hitdb database
            hitsdb.push_back(h);
        }
    //we clean for the moment the local container
    localhits.erase(localhits.begin(),localhits.end());

}

//Construct contigs from FASTA file
void SPolisher::read_fasta(string filename, vector<ctg> & seqs){
    std::ifstream input(filename);
    if(!input.good()){
        std::cerr << "Error opening '"<<filename<<"'. Bailing out." << std::endl;
        exit(1);
    }
    int id=0;
    uint32_t ori=0;
    std::string line, name, content;
    while(getline( input, line ).good() ){
        if( line.empty() || line[0] == '>' ){ // Identifier marker
            if( !name.empty() ){ // Print out what we read from the last entry
                //std::cout << name << ":" << content << std::endl;
                //cout << name << ":" <<id<<":"<< content << std::endl;
                //todo: improve the parsing of contig name for the moment we expect and space
                ctg tmp;
                tmp.ctgid=id;
                tmp.namectg=name.substr(0, name.find(' '));
                tmp.length= static_cast<int>(content.length());
                tmp.seq=content;
                seqs.push_back(tmp);
                id++;
                name.clear();

            }
            if( !line.empty() ){
                name = line.substr(1);//pick name before the end;
            }
            content.clear();
        } else if( !name.empty() ){
            if( line.find(' ') != std::string::npos ){ // Invalid sequence--no spaces allowed
                name.clear();
                content.clear();
            } else {
                content += line;
            }
        }
    }

    if( !name.empty() ){ // Print out what we read from the last entry
        //cout << name << ":" <<id<<":"<< content << std::endl;
        ctg tmp;
        tmp.ctgid=id;
        tmp.namectg=name.substr(0, name.find(' '));
        tmp.length= static_cast<int>(content.length());
        tmp.seq=content;
        seqs.push_back(tmp);
    }
}



//minimap steps:
//For each cluster, sort (qi,ti) by qi and solve a longest increasing sequence problem for ti. This finds the longest co-linear matching chain. Break the chain whenever there is a gap longer than -g=10000.

//Output the start and end of the chain if it contains -c=4 or more minimizer matches and the matching length is no less than -L=40.

/*static inline void mm_cal_fuzzy_len(mm_reg1_t *r, const mm128_t *a)
{
    int i;
    r->mlen = r->blen = 0;
    if (r->cnt <= 0) return;
    r->mlen = r->blen = a[r->as].y>>32&0xff;
    for (i = r->as + 1; i < r->as + r->cnt; ++i) {
        int span = a[i].y>>32&0xff;
        int tl = (int32_t)a[i].x - (int32_t)a[i-1].x;
        int ql = (int32_t)a[i].y - (int32_t)a[i-1].y;
        r->blen += tl > ql? tl : ql;
        r->mlen += tl > span && ql > span? span : tl < ql? tl : ql;
    }
}*/

//float identity = (float)r->mlen / r->blen;

/* API DOCS
typedef struct {
int32_t id;             // ID for internal uses (see also parent below)
int32_t cnt;            // number of minimizers; if on the reverse strand
int32_t rid;            // reference index; if this is an alignment from inversion rescue
int32_t score;          // DP alignment score
int32_t qs, qe, rs, re; // query start and end; reference start and end
int32_t parent, subsc;  // parent==id if primary; best alternate mapping score
int32_t as;             // offset in the a[] array (for internal uses only)
int32_t mlen, blen;     // seeded exact match length; seeded alignment block length
int32_t n_sub;          // number of suboptimal mappings
int32_t score0;         // initial chaining score (before chain merging/spliting)
 */

//FAST-SG channing
/*vector<hit> compute_score(vector<hit> &tfwd, int seqlen){
    //we sort the array by ctg star
    sort(tfwd.begin(),tfwd.end(),compare_by_ctg_pos);
    hit mfwd=tfwd[0];
    int mc=0;
    vector<hit> index;
    int add_last=0;
    for (int i = 0; i < tfwd.size() ; i++) {
        //means that both kmers comes from the same region
        if((mfwd.ctg == tfwd[i].ctg) && (tfwd[i].pos-mfwd.pos <= seqlen)){
            mc++;
        }else{  //we add the last element to the array
            mfwd.score=mc;
            index.push_back(mfwd);
            //we add a new mfwd
            mc=1;
            mfwd=tfwd[i];
            mfwd.score=1;
            if(i == tfwd.size() - 1){
                add_last=1;
            }
        }
    }
    //we add the last kmer to the set because is has an score of 1
    if(add_last){
        index.push_back(mfwd);
    }
    //we add the kmer to the set because it has score max
    if(mc == tfwd.size()){
        mfwd.score=mc;
        index.push_back(mfwd);
    }
    //we sort by score in order to print the output
    sort(index.begin(),index.end(),compare_score);
    return index;
}*/


/*minimap

static void proc_intv(mm_tbuf_t *b, int which, int k, int min_cnt, int max_gap)
{
	int i, j, l_lis, rid = -1, rev = 0, start = b->intv.a[which].y, end = start + b->intv.a[which].x;

	// make room for arrays needed by LIS (longest increasing sequence)
	if (end - start > b->m) {
		b->m = end - start;
		kv_roundup32(b->m);
		b->a = (uint64_t*)realloc(b->a, b->m * 8);
		b->b = (size_t*)realloc(b->b, b->m * sizeof(size_t));
		b->p = (size_t*)realloc(b->p, b->m * sizeof(size_t));
	}

	// prepare the input array _a_ for LIS
	b->n = 0;
	for (i = start; i < end; ++i)
		if (b->coef.a[i].x != UINT64_MAX)
			b->a[b->n++] = b->coef.a[i].y, rid = b->coef.a[i].x << 1 >> 33, rev = b->coef.a[i].x >> 63;
	if (b->n < min_cnt) return;
	radix_sort_64(b->a, b->a + b->n);

	// find the longest increasing sequence
	l_lis = rev? ks_lis_low32gt(b->n, b->a, b->b, b->p) : ks_lis_low32lt(b->n, b->a, b->b, b->p); // LIS
	if (l_lis < min_cnt) return;
	for (i = 1, j = 1; i < l_lis; ++i) // squeeze out minimizaers reused in the LIS sequence
		if (b->a[b->b[i]]>>32 != b->a[b->b[i-1]]>>32)
			b->a[b->b[j++]] = b->a[b->b[i]];
	l_lis = j;
	if (l_lis < min_cnt) return;

	// convert LISes to regions; possibly break an LIS at a long gaps
	for (i = 1, start = 0; i <= l_lis; ++i) {
		int32_t qgap = i == l_lis? 0 : ((uint32_t)b->mini.a[b->a[b->b[i]]>>32].y>>1) - ((uint32_t)b->mini.a[b->a[b->b[i-1]]>>32].y>>1);
		if (i == l_lis || (qgap > max_gap && abs((int32_t)b->a[b->b[i]] - (int32_t)b->a[b->b[i-1]]) > max_gap)) {
			if (i - start >= min_cnt) {
				uint32_t lq = 0, lr = 0, eq = 0, er = 0, sq = 0, sr = 0;
				mm_reg1_t *r;
				kv_pushp(mm_reg1_t, b->reg, &r);
				r->rid = rid, r->rev = rev, r->cnt = i - start, r->rep = 0;
				r->qs = ((uint32_t)b->mini.a[b->a[b->b[start]]>>32].y>>1) - (k - 1);
				r->qe = ((uint32_t)b->mini.a[b->a[b->b[i-1]]>>32].y>>1) + 1;
				r->rs = rev? (uint32_t)b->a[b->b[i-1]] : (uint32_t)b->a[b->b[start]];
				r->re = rev? (uint32_t)b->a[b->b[start]] : (uint32_t)b->a[b->b[i-1]];
				r->rs -= k - 1;
				r->re += 1;
				for (j = start; j < i; ++j) { // count the number of times each minimizer is used
					int jj = b->a[b->b[j]]>>32;
					b->mini.a[jj].y += 1ULL<<32;
					kv_push(uint32_t, b->reg2mini, jj); // keep minimizer<=>reg mapping for derep
				}
				for (j = start; j < i; ++j) { // compute ->len
					uint32_t q = ((uint32_t)b->mini.a[b->a[b->b[j]]>>32].y>>1) - (k - 1);
					uint32_t r = (uint32_t)b->a[b->b[j]];
					r = !rev? r - (k - 1) : (0x80000000U - r);
					if (r > er) lr += er - sr, sr = r, er = sr + k;
					else er = r + k;
					if (q > eq) lq += eq - sq, sq = q, eq = sq + k;
					else eq = q + k;
				}
				lr += er - sr, lq += eq - sq;
				r->len = lr < lq? lr : lq;
			}
			start = i;
		}
	}
}*/