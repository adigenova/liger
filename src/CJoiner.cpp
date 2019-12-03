//
// Created by Alex Digenova on 2/15/19.
//

#include "CJoiner.h"



cehits CJoiner::TrimALN(EdlibAlignResult  &result, int targetL, int queryL, int wsize, float tiden) {
    //Step 1: compute the offset of the window
    int w=wsize; //window for computing identity of the two sequences
    int good_w=0;
    int total_w=0;
    string wcigar="";
    //int bad_w=0;
    //aligment variables
    int m=0,i=0,d=0,s=0;//match, insertion, deletion, snps : bases in target= m+s+d bases in query = m+s+i  //0 based
    char moveCodeToChar[] = {'M', 'I', 'D', 'S'};
    int alnstart=0;
    int qb=0;//qbases
    int tb=result.startLocations[0];//tbases

    m=0,i=0,d=0,s=0;
    //temporal vector to store windows
    vector<cehits> localw;
    cehits rtrim;
    //int j=0;
    //step 2: check each aligment carefully by window in the target
    for(auto j=alnstart; j<result.alignmentLength; j++){
        switch(result.alignment[j]) {
            case 0 :  //Match
                qb++; tb++;
                m++;
                break;
            case 1 : i++;//Insertion
                qb++;
                break;
            case 2 : //Deletion
                d++;
                tb++;
                break;
            case 3 : //s++;subtitution
                s++;
                qb++; tb++;
                break;
        }
        //means that we will change of window in the nex base, so  we compute the relevant variables
        if((tb)%w == 0 && (m+d+s)>0){
            int dq=m+i+s;
            int dt=m+d+s;
            int daln=m+s+i+d;
            bool type= (float)(m) / (daln) >= tiden;
            bool partial =  dt == w;
            total_w++;
            good_w+=type;
            wcigar+= type == true ? "M":"S";

            cehits tmp;
            tmp.w=(tb-1)/w;
            tmp.qstart=qb-dq;
            tmp.qlen=dq;
            tmp.tstart=tb-dt;
            tmp.tlen=dt;
            tmp.alive=type;
            tmp.partial=partial;
            tmp.iden=(float)m/daln;
            //we store in the current vector container
            localw.push_back(tmp);
            //we comback this variables to 0
            m=0,i=0,d=0,s=0;
        }
    }

    //we have to check the last window of the sequence
    if((m+d+s)>0){
        int dq=m+i+s;
        int dt=m+d+s;
        int daln=m+s+i+d;

        bool type= (float)(m) / (daln) >= tiden;
        bool partial =  dt == w;

        cehits tmp;
        tmp.w=(tb-1)/w;
        tmp.qstart=qb-dq;
        tmp.qlen=dq;
        tmp.tstart=tb-dt;
        tmp.tlen=dt;
        tmp.alive=type;
        tmp.partial=partial;
        tmp.iden=(float)m/daln;
        //we store in the current vector container
        localw.push_back(tmp);
        //we comback this variables to 0
        m=0,i=0,d=0,s=0;
        total_w++;
        good_w+=type;
        wcigar+= type == true ? "M":"S";
    }
    //step 3: add good windows to a vector
    cout << "TOTALW="<<total_w<<" GOODW="<<good_w<<" "<<" BADW="<<total_w-good_w<<" CIGARW="<<wcigar<<" SW="<<result.startLocations[0]/w<<endl;
    if(good_w > 0){
        //we remove the missmatch from the begining and end of the vector
        int32_t begin = 0, end = static_cast<int32_t>(localw.size() - 1);
        for (; begin < localw.size(); ++begin) {
            if (localw[begin].alive)
                break;
        }
        for (; end >= 0; --end) {
            if (localw[end].alive)
                break;
        }
        //we store the good query windows in the global array of windows
        cout << "TOTALW="<<total_w<<" GOODW="<<good_w<<" "<<" BADW="<<total_w-good_w<<" CIGARC="<<wcigar.substr(begin,(end-begin)+1)<<" SW="<<result.startLocations[0]/w
             <<" B="<<begin<<" END="<<end<<endl;

        rtrim.tstart=localw[begin].tstart;
        rtrim.tstop=localw[end].tstart+localw[end].tlen;
        rtrim.qstart=localw[begin].qstart;
        rtrim.qstop=localw[end].qstart+localw[end].qlen;
        for(auto k=begin; k<=end; k++)
            rtrim.iden+=localw[k].iden;

        rtrim.iden=rtrim.iden/(abs(begin-end)+1);
        rtrim.qcov=(100 *abs(rtrim.qstart-rtrim.qstop)/queryL);
        //we clean the local container
        localw.erase(localw.begin(),localw.end());
        return rtrim;
        /*for(auto k=begin; k<=end; k++)
            this->windows[localw[k].w].push_back(localw[k]);*/
        //we have to create a new cehits with is the union of the
    }
    //we clean the container
    localw.erase(localw.begin(),localw.end());
    //return good_w;
    return rtrim;
}



cehits CJoiner::compute_optimal_aligment_trim(string &q, string &t, bool log, int cns_cov) {
    //return cehits();
    //cehits hit;
    cehits rtrim;

    EdlibAlignResult result;
    result = edlibAlign(q.c_str(), int(q.length()), t.c_str(),
                        int(t.length()),
                        edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL,
                                            0));
    //we need to compute the cigar line to adjust the variables in the query sequence
    if (result.status == EDLIB_STATUS_OK) {
        //for visualization only
        if(log) {
            this->printAlignment(q.c_str(), t.c_str(), result.alignment, result.alignmentLength, *(result.endLocations),
                                 EDLIB_MODE_HW);
        }

        //todo: trim alingment when the query sequence start/end with gaps, likely a partial alignment or homopolymer seq useful with nanopore
        float min_iden=0.8;
        if(cns_cov >=5)
            min_iden=0.9;
        if(cns_cov >=10)
            min_iden=0.95;

        //we trim the alignment in order to specify the minimum identity
        rtrim=this->TrimALN(result,t.length(),q.length(),30,min_iden);
        printf("qs=%d qe=%d qc=%d ts=%d te=%d iden=%f\n",rtrim.qstart,rtrim.qstop,rtrim.qcov,rtrim.tstart,rtrim.tstop,rtrim.iden);



    } else {
        printf("Error computing edit distance on CJoiner class\n");
        exit(1);
    }
    //we clean the edlib object
    edlibFreeAlignResult(result);
    //we the boundaries and identity for the current location
    return rtrim;
}


cehits CJoiner::compute_optimal_aligment(string &q, string &t, bool log) {
    //return cehits();
    cehits hit;
    EdlibAlignResult result;
        result = edlibAlign(q.c_str(), int(q.length()), t.c_str(),
                            int(t.length()),
                            edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL,
                                                0));
    //we need to compute the cigar line to adjust the variables in the query sequence
    if (result.status == EDLIB_STATUS_OK) {
        //for visualization only
        if(log) {
            this->printAlignment(q.c_str(), t.c_str(), result.alignment, result.alignmentLength, *(result.endLocations),
                                 EDLIB_MODE_HW);
        }

        //todo: trim alingment when the query sequence start/end with gaps, likely a partial alignment or homopolymer seq useful with nanopore
        hit.tstart=result.startLocations[0];
        hit.tstop=result.endLocations[0];
        //identity
        hit.identity=100-(100*result.editDistance/result.alignmentLength);
        hit.qcov=(100 * abs(hit.tstart-hit.tstop)/q.length());
        hit.edit_distance=result.editDistance;
        //due to gaps the lenght of the alignment could be 100
        if(hit.qcov>100){
            hit.qcov=100;
        }

    } else {
        printf("Error computing edit distance on CJoiner class\n");
        exit(1);
    }
    //we clean the edlib object
    edlibFreeAlignResult(result);
    //we the boundaries and identity for the current location
    return hit;
}

//function that compute N-Optimal alignment of the query in the target sequence, it stop when the identity drop to less than 0.65
//it mask any optimal aligment
vector<cehits> CJoiner::compute_n_optimal_aligment(string &query, string &target, int n) {
    //return vector<cehits>();
    vector<cehits> all_aligments;
    string targetseq=target;
    //find the best optimal alignment between the query and target seqs
    auto best=this->compute_optimal_aligment(query, targetseq, false);
    all_aligments.push_back(best);
    //we stop sercing because the best is already bad
    if(best.identity < 75)
        return all_aligments;
    //we mask the target seq to find another hit
    targetseq.replace(best.tstart,best.tstop,abs(best.tstop-best.tstart),'N');
    int i=0;
    while(i<n) {
        cehits next = this->compute_optimal_aligment(query, targetseq, false);
        //we stop searching because identity is really bad
        if(next.identity < 75)
            return all_aligments;
        all_aligments.push_back(next);
        //we mask the target sequence to find another aligment
        targetseq.replace(next.tstart,next.tstop,abs(next.tstop-next.tstart),'N');
        i++;
    }
    return all_aligments;
}



//we try to fix the edge here
void CJoiner::contigs_best_fit(string &source_end, string &target_end, string &consensus, string &coverage_s, string &longread, EdgeS *e, bool revcns, bool leadrev) {
    //we maps the end of the contigs to the consensus sequence

    auto s2c = this->compute_optimal_aligment(source_end, consensus, false);
    auto t2c = this->compute_optimal_aligment(target_end, consensus, false);

    //we compute the onserved distance
    int d_obs=this->get_obs_distance_hits(&s2c,&t2c);

    //we compare the d_obs to the expected one considering the std or a maximum desviation of 10%
    // u+- 3*std= 99.7% of cases
    if(abs(d_obs-e->getBd()) < 4 * max(e->getBd_std(),(int)(e->getBd()*0.1+0.5))){
                //the edge is ok, we have to fill the variables
                this->fill_edge_variables(d_obs,revcns,e,consensus,coverage_s,s2c,t2c);
    }else{
        //we have to rescue the edge, first we use multiple mapping of the ends
        auto hsource=this->compute_n_optimal_aligment(source_end,consensus,5);
        auto htarget=this->compute_n_optimal_aligment(target_end,consensus,5);
        //we check all the possible combinations and choose the most close to the expected distance
        int selected_a=-1;
        int selected_b=-1;
        int best_dist=d_obs;
        //we find the best pair
        for (auto a=0; a<hsource.size(); ++a){
            //we compute the distance among the pairs
            for(auto b=0; b<htarget.size(); ++b){
                int d2=this->get_obs_distance_hits(&hsource[a],&htarget[b]);
                 //we pick the new distance if it better fit the previous one
                if(abs(d2-e->getBd()) < abs(best_dist-e->getBd())){
                    best_dist=d2;
                    selected_a=a;
                    selected_b=b;
                }
            }
        }
        //we check if some of the distances were selected and they are in the expected range
        if((selected_a >= 0 and selected_b >=0) and (abs(best_dist-e->getBd()) < 4 * max(e->getBd_std(),(int)(e->getBd()*0.1+0.5)))) {
            //cout << "S A=" << selected_a << " S B=" << selected_b << "S D=" << best_dist << endl;
            this->fill_edge_variables(best_dist,revcns,e,consensus,coverage_s,hsource[selected_a],htarget[selected_b]);
        }else{
            //our last attempt is to try the lead sequence
            auto local_lonread=longread;

            int selected_al=-1;
            int selected_bl=-1;
            int best_distl=d_obs;
            vector<cehits> hsourcel, htargetl;
            if(e->getBd() > 0){
                 if(leadrev == 1)
                    revcomp(local_lonread);

                hsourcel=this->compute_n_optimal_aligment(source_end,local_lonread,5);
                htargetl=this->compute_n_optimal_aligment(target_end,local_lonread,5);
                //we check all the possible combinations and choose the most close to the expected distance
                //we find the best path
                for (auto a=0; a<hsourcel.size(); ++a){
                    //we compute the distance among the pairs
                    for(auto b=0; b<htargetl.size(); ++b){
                        //we skypt low quality matchs
                        if(hsourcel[a].identity < 75 || htargetl[b].identity < 75)
                            continue;
                        int d2=this->get_obs_distance_hits(&hsourcel[a],&htargetl[b]);
                        //we pick the new distance if it better fit the previous one
                        if(abs(d2-e->getBd()) < abs(best_distl-e->getBd())){
                            best_distl=d2;
                            selected_al=a;
                            selected_bl=b;
                        }
                    }
                }
            }
            //we check if we found a good result
             if((selected_al>=0 and selected_bl >=0) and (abs(best_distl-e->getBd()) < 4 * max(e->getBd_std(),(int)(e->getBd()*0.1+0.5)))){


                 char c=(char)1;//there is only one long read
                 string long_reads_coverage_s (local_lonread.length(), c);
                 this->fill_edge_variables(best_distl,revcns,e,local_lonread,long_reads_coverage_s,hsourcel[selected_al],htargetl[selected_bl]);
             }
             else{
                    //we pick the one from multiples alignments to the consensus or the optimal alignments
                    if(selected_a >=0 and selected_b >=0)
                        this->fill_edge_variables(best_dist,revcns,e,consensus,coverage_s,hsource[selected_a],htarget[selected_b]);
                    else
                        this->fill_edge_variables(d_obs,revcns,e,consensus,coverage_s,s2c,t2c);
             }
        }

    }

}


int CJoiner::get_obs_distance_hits(cehits *s2c, cehits *t2c){
    int d_obs=0;
    // no overlap
    if(s2c->tstop < t2c->tstart)
        d_obs = abs(t2c->tstart - s2c->tstop);
    else
        //overlap shorter than the contigs ends < 0.5kb
    if (s2c->tstop >= t2c->tstart && s2c->tstop <= t2c->tstop)
        d_obs = -abs(s2c->tstop - t2c->tstart);
    else    //weird case when the source is aligned after the target could happen when the overlap is larger than the 500bp
    if (t2c->tstart < s2c->tstart)
        d_obs = -abs(t2c->tstart - s2c->tstop);
    else
        if (s2c->tstop >= t2c->tstart && s2c->tstop >= t2c->tstop)
            d_obs = -abs(s2c->tstop - t2c->tstart);
        else
            printf("Warning case not take into account s.b=%d s.e=%d t.b=%d t.e=%d\n",s2c->tstart,s2c->tstop,t2c->tstart,t2c->tstop);

    return d_obs;
}


void CJoiner::fill_edge_variables(int distance_obs,bool revcns, EdgeS* e , string &consensus, string &coverage_s,cehits sh, cehits th){


    e->setCovs(sh.qcov);
    e->setCovt(th.qcov);
    e->setIdens(sh.identity);
    e->setIdent(th.identity);

    if(sh.identity < 70 || th.identity < 70)
        e->setEdge_low_identity(true);
    else
        e->setEdge_proper_gap(true);


    if (distance_obs <= 0) {
        e->setEdge_overlap(distance_obs);
    }else{
        //we have to assert that the source is previous to the target
      

        auto cstart=sh.tstop+1;
        auto cstop=th.tstart-1;

        //we mark that the edge is lowIdentity

        //we save the consensus sequence
        auto gapseq=consensus.substr(cstart,abs(cstart-cstop)+1);
        auto gapdepth=coverage_s.substr(cstart,abs(cstart-cstop)+1);

        //we save the gapseq in the strand of the edge, so we then revert when we traverse the line
        if(revcns == 0) {
            e->set_edge_gapseq(revcomp(gapseq).c_str());
            //todo: set the edge coverage sequence for report quality values later
            reverse(gapdepth.begin(),gapdepth.end());
        }
        else {
            e->set_edge_gapseq(gapseq.c_str());
        }

        int avg_cov=0;
        //we compute the average coverage and we create another vector holding the average depth
        for (int i = 0; i <gapdepth.length() ; ++i) {
            //avg_cov+=coverage[i];
            avg_cov+=(int)gapdepth[i];
        }
        avg_cov= static_cast<int>((float)avg_cov / (float)(gapdepth.length()));

        e->setEdge_avg_cns(avg_cov);
        e->setEdge_overlap(gapseq.length());
    }



}



void CJoiner::printAlignment(const char* query, const char* target, const unsigned char* alignment, const int alignmentLength, const int position, const EdlibAlignMode modeCode) const {
    int tIdx = -1;
    int qIdx = -1;

    if (modeCode == EDLIB_MODE_HW) {
        tIdx = position;
        for (int i = 0; i < alignmentLength; i++) {
            if (alignment[i] != EDLIB_EDOP_INSERT)
                tIdx--;
        }
    }
    printf("\n");
    for (int start = 0; start < alignmentLength; start += 50) {
        // target
        printf("T: ");
        int startTIdx;
        for (int j = start; j < start + 50 && j < alignmentLength; j++) {
            if (alignment[j] == EDLIB_EDOP_INSERT)
                printf("-");
            else
                printf("%c", target[++tIdx]);
            if (j == start)
                startTIdx = tIdx;
        }
        printf(" (%d - %d)\n", max(startTIdx, 0), tIdx);
        // match / mismatch
        printf("   ");
        for (int j = start; j < start + 50 && j < alignmentLength; j++) {
            printf(alignment[j] == EDLIB_EDOP_MATCH ? "|" : " ");
        }
        printf("\n");
        // query
        printf("Q: ");
        int startQIdx = qIdx;
        for (int j = start; j < start + 50 && j < alignmentLength; j++) {
            if (alignment[j] == EDLIB_EDOP_DELETE)
                printf("-");
            else
                printf("%c", query[++qIdx]);
            if (j == start)
                startQIdx = qIdx;
        }
        printf(" (%d - %d)\n\n", max(startQIdx, 0), qIdx);
    }
    printf("\n");
}
