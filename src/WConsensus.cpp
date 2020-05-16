//
// Created by Alex Digenova on 2/1/19.
//

#include "WConsensus.h"


//constructor of the class
WConsensus::WConsensus(vector<linking> &lin, vector<string> &seqs, string &leadseq, uint32_t lid, int wsize) {
    //we set the window size
    this->windowsize=wsize;
    //we prepare the lead sequence
    int index=-1;
    string ls;
    for(int i=0; i<lin.size(); i++){
        auto l=lin[i];
        //we found the lead sequence
        if(l.lonread_id == lid){
            int start= l.gaps-1000 > 0 ? l.gaps-1000:l.gaps;
            int stop= l.gape+1000 < leadseq.length() ? l.gape+1000:l.gape;
            ls=leadseq.substr(start,stop-start);
            if(l.switched)
                ls=revcomp(ls);
            index=i;
            break;
        }
    }
     	
    //we populate the vector with the consensus vector
    for(int i=0; i<ls.length(); i+=wsize){
        cwindow tmp;
        tmp.w=i/wsize;
        tmp.qstart=-1;
        tmp.qlen=-1;
        tmp.tstart=i;
        //we check that we dont pass the leadseq length
        tmp.tlen= static_cast<int>(tmp.tstart + wsize < ls.length() ? wsize : ls.length() - tmp.tstart);
        tmp.alive=1;
        tmp.partial=0;
        tmp.iden=1.0;
        tmp.qindex=-1;
        vector<cwindow> v;
        v.push_back(tmp);
        //we store in the current vector container
        this->windows.push_back(v);
    }

    /*for( auto w:this->windows){
        for( auto item:w){
            cout << item.w <<" "<<item.tstart<<" "<<item.tlen<<" "<<item.alive<<" "<<item.iden<<" "<<ls.length()<<" "<<leadseq.length()<<endl;
        }
    }*/

}


//it create a consensus from a lead sequence, it compute an aligment of the non-lead sequences agains the lead to determine the boundaries
string WConsensus::consensusfromlead(vector<linking> &lin, vector<string> &seqs, string &leadseq, uint32_t lid) {

     //we prepare the lead sequence
    int lmin=lin[0].p1;
    int lmax=0;
    int index=0;
    string ls="";
     for(int i=0; i<lin.size(); i++){
         auto l=lin[i];
         if(l.p1 < lmin )
             lmin=l.p1;
          if(l.p1+l.gape-l.gaps > lmax)
              lmax=l.p1+l.gape-l.gaps;
         //we found the lead sequence
         if(l.lonread_id == lid){
             int start= l.gaps-1000 > 0 ? l.gaps-1000:l.gaps;
             int stop= l.gape+1000 < leadseq.length() ? l.gape+1000:l.gape;
             ls=leadseq.substr(start,stop-start);
             if(l.switched)
                 ls=revcomp(ls);
             index=i;
             break;
            /* cout << "*****"<<l.c1 <<" "<<l.c2<<" P1="<<l.p1<<" P2="<<l.p2<<" START="<<l.p1<<" END="<<l.p1+(l.gape-l.gaps)<<" LG="
                  <<(l.gape-l.gaps)<<" LR="<<l.lonread_id<<" SLR="<<l.switched<<" LRGS="<<l.gaps<<" LRGE="<<l.gape<<endl;*/
         }/*else{
             cout << l.c1 <<" "<<l.c2<<" P1="<<l.p1<<" P2="<<l.p2<<" START="<<l.p1<<" END="<<l.p1+(l.gape-l.gaps)<<" LG="
                  <<(l.gape-l.gaps)<<" LR="<<l.lonread_id<<" SLR="<<l.switched<<" LRGS="<<l.gaps<<" LRGE="<<l.gape<<endl;
         }*/
     }

    //step two we align the other seqs to the lead sequence for creating windows of length k
    for(int i=0; i<lin.size(); i++){
        //we skypt the lead sequence aligned to itself
        if(lin[i].lonread_id == lid)
            continue;

        EdlibAlignResult result;
        result = edlibAlign(seqs[i].c_str(), int(seqs[i].length()), ls.c_str(),
                                int(ls.length()),
                                edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL,
                                                    0));

        if (result.status == EDLIB_STATUS_OK) {
            /*printf("ALN leadseq with: %d TL=%d QL=%d ED=%d ALN_L=%d TS=%d TE=%d \n", i, leadseq.length(), seqs[i].length(),
                   result.editDistance,result.alignmentLength, result.startLocations[0],result.endLocations[0]);
*/
            auto alnlen=this->CheckAlignment3(result,leadseq.length(),seqs[i].length(), ls, seqs[i], i);
            //this->printAlignment(seqs[i].c_str(), leadseq.c_str(), result.alignment, result.alignmentLength,result.alignmentLength,
            //                     *(result.endLocations), EDLIB_MODE_HW);
        }
        edlibFreeAlignResult(result);
    }

    string consensus2="";
    //we build the consensus sequence
    for( auto i=0; i<this->windows.size(); ++i){
        //cout << "W: "<<i<<" Size="<<this->windows[i].size()<<endl;
        consensus2+=this->generate_consensus2(this->windows[i],seqs,ls);
    }


    //todo: we have to clean the local containers
    //windows.clear();
     return consensus2;
}

string WConsensus::generate_consensus2(vector<cwindow> & window, vector<string> &seqs, string & leadseq){


    //main diference is the type of alignment on the window, NW is global and SW is local
    //we create the alignment type, probably shoutl be in the previous method
/*
    auto alignment_local = spoa::createAlignmentEngine(spoa::AlignmentType::kSW,
                                                        5, -4, -8);
*/

    //racon use the following aln engines
    //spoa::createAlignmentEngine(spoa::AlignmentType::kNW, match, mismatch, gap);
    //where match, mismatch ang gap are     int8_t match = 5; int8_t mismatch = -4; int8_t gap = -8;
    //in linux I cannot have the two modes of alignment
   auto alignment_global = spoa::createAlignmentEngine(spoa::AlignmentType::kSW,
                                                        5, -4, -8);
    //we prealloc the memory of the alignment engine same as in racon
    alignment_global->prealloc(this->windowsize, 5);

    string consensus_="";
    if (window.size() < 3) {
        consensus_ = leadseq.substr(window[0].tstart, window[0].tlen);
        for(auto j=0; j< window[0].tlen; j++ )
            this->coverage.push_back(1);
        return consensus_;//we return the best read sequence
    }

    //we build a consensus
    auto graph = spoa::createGraph();

    /*auto alignment = alignment_global->align_sequence_with_graph(leadseq.substr(window[0].tstart, window[0].tlen), graph);
    cout <<leadseq.length()<<" "<<window[0].tstart<<" "<<window[0].tlen<<endl;
    graph->add_alignment(alignment, leadseq.substr(window[0].tstart, window[0].tlen));
    */
    assert(window[0].tstart+window[0].tlen <= leadseq.length());

    auto lsub=leadseq.substr(window[0].tstart, window[0].tlen);
    //graph->add_alignment(spoa::Alignment(), lsub, 1);
    auto alignment = alignment_global->align_sequence_with_graph(lsub, graph);
    graph->add_alignment(alignment,lsub);


    int global=0;
    int local=0;
    int sbegin=0;
    int send=0;
    for(auto item:window){
        //we have already aligned the reference sequence at the begining
        if(item.qindex < 0)
            continue;
        //we don't shorter sequences
        if(item.qlen < 50)
            continue;
        //the index should be lower than the number of seqs
        assert(item.qindex < seqs.size());
        //the length of the substring should be shorter than the total length
        assert(item.qstart+item.qlen <= seqs[item.qindex].length());

        auto s2g=seqs[item.qindex].substr(item.qstart, item.qlen);
     /*   cout << item.w <<" "<<item.tstart<<" "<<item.tlen<<" "<<item.alive<<" "<<item.iden
             << item.qstart <<" "<<item.qlen<<" "<<item.qindex<<" "<<item.partial<<endl;
*/
        //true means a complete aligned sequence
        spoa::Alignment alignment;
        if(item.partial) {
            alignment = alignment_global->align_sequence_with_graph(s2g, graph);
            global++;
            send++;
            sbegin++;
        }else{
            //we have perform a global alignment restricted to the local window
            alignment = alignment_global->align_sequence_with_graph(s2g, graph);
            local++;
            if(item.w*this->windowsize == item.tstart)
                sbegin++;
            else
                send++;
        }
        graph->add_alignment(alignment, s2g);
    }

   // cout << "sbegin="<<sbegin<<" send="<<send<<endl;

    vector<uint32_t> coverages;
    consensus_ = graph->generate_consensus(coverages);

    uint32_t average_begin = (sbegin) / 2;

    int32_t begin = 0, end = consensus_.size() - 1;
    for (; begin < static_cast<int32_t>(consensus_.size()); ++begin) {
        if (coverages[begin] >= average_begin) {
            break;
        }
    }

    uint32_t average_end = (send) / 2;
    for (; end >= 0; --end) {
        if (coverages[end] >= average_end) {
            break;
        }
    }

    if (begin >= end) {
        fprintf(stderr, "[Wengan::WConsensus::generate_consensus2] warning: "
                "lead layout might be chimeric in window !\n");
        //we return the lead sequence window unpolished
        consensus_ = leadseq.substr(window[0].tstart, window[0].tlen);
        for(auto j=0; j< window[0].tlen; j++ )
            this->coverage.push_back(1);
        //return consensus_;//we return the best read sequence
    } else {
        consensus_ = consensus_.substr(begin, end - begin + 1);
        for(auto j=begin; j<=end ; j++ )
            this->coverage.push_back(coverages[j]);
       // cout << "BEGIN= "<<begin<<" END="<<end<<" CNS="<<consensus_<<endl;
    }

    //we release the memory allocated by the spoa graph those are smart pointers, in theory reset mark
    graph.reset();
    alignment_global.reset();
    return consensus_;
}


int WConsensus::CheckAlignment3(EdlibAlignResult  &result, int targetL, int queryL, string  &leadseq, string &queryseq, int qindex) {
    //Step 1: compute the offset of the window
    int w=this->windowsize; //window for computing identity of the two sequences
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
    vector<cwindow> localw;
    //int j=0;
    //step 2: check each aligment carefully by window in the target
    for(auto j=alnstart; j<result.alignmentLength; j++){
        //cout <<moveCodeToChar[result.alignment[j]]<<" "<<j<<" "<<tb/w<<" "<<tb%w<<" ";
        switch(result.alignment[j]) {
            case 0 :  //Match
                //cout <<queryseq[qb]<<"-"<<leadseq[tb]<<" "<<qb<<" "<<tb<<endl;
                qb++; tb++;
                m++;
                break;
            case 1 : i++;//Insertion
                //cout <<queryseq[qb]<<" "<<"-"<<" "<<qb<<" "<<tb<<endl;
                qb++;
                break;
            case 2 : //Deletion
                //cout <<"-"<<" "<<leadseq[tb]<<" "<<qb<<" "<<tb<<endl;
                d++;
                tb++;
                break;
            case 3 : //s++;subtitution
                //cout <<queryseq[qb]<<" "<<leadseq[tb]<<" "<<qb<<" "<<tb<<endl;
                s++;
                qb++; tb++;
                break;
        }
        //means that we will change of window in the nex base, so  we compute the relevant variables
        if((tb)%w == 0 && (m+d+s)>0){
            int dq=m+i+s;
            int dt=m+d+s;
            int daln=m+s+i+d;
            bool type= (float)(m) / (daln) >= 0.65;
            bool partial =  dt == w;
            total_w++;
            good_w+=type;
            wcigar+= type == true ? "M":"S";
            //cout <<endl;
            /*cout << "T="<<tb-dt<<"-"<<tb<<" Q="<<qb-dq<<"-"<<qb<<" "
                 <<"W="<<(tb-1)/w<<" BT="<<dt<<" BQ="<<dq<<" ALN="<<daln
                 <<" I1="<< (float)m/daln<<" TYPE="<<type<<" PARTIAL="<<partial<<endl;*/
            //cout <<"T="<<leadseq.substr(tb-dt,dt)<<endl;
            //cout <<"Q="<<queryseq.substr(qb-dq,dq)<<endl;
            //cout <<endl;
            //we store the window in the
            cwindow tmp;
            tmp.w=(tb-1)/w;
            tmp.qstart=qb-dq;
            tmp.qlen=dq;
            tmp.tstart=tb-dt;
            tmp.tlen=dt;
            tmp.alive=type;
            tmp.partial=partial;
            tmp.iden=(float)m/daln;
            tmp.qindex=qindex;
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

     /* cout <<"LastW:"<<endl;
        cout << "T="<<tb-dt<<"-"<<tb<<" Q="<<qb-dq<<"-"<<qb<<endl;
        cout <<"W="<<(tb-1)/w<<" BT="<<dt<<" BQ="<<dq<<" ALN="<<daln<<" I1="<< (float)m/daln<<endl;
        cout <<"T="<<leadseq.substr(tb-dt,dt)<<endl;
        cout <<"Q="<<queryseq.substr(qb-dq,dq)<<endl;
        cout <<endl;
     */
        bool type= (float)(m) / (daln) >= 0.65;
        bool partial =  dt == w;
        //cout <<endl;
        /*cout << "T="<<tb-dt<<"-"<<tb<<" Q="<<qb-dq<<"-"<<qb<<" "
             <<"W="<<(tb-1)/w<<" BT="<<dt<<" BQ="<<dq<<" ALN="<<daln
             <<" I1="<< (float)m/daln<<" TYPE="<<type<<" PARTIAL="<<partial<<endl;*/
        cwindow tmp;
        tmp.w=(tb-1)/w;
        tmp.qstart=qb-dq;
        tmp.qlen=dq;
        tmp.tstart=tb-dt;
        tmp.tlen=dt;
        tmp.alive=type;
        tmp.partial=partial;
        tmp.iden=(float)m/daln;
        tmp.qindex=qindex;
        //we store in the current vector container
        localw.push_back(tmp);
        //we comback this variables to 0
        m=0,i=0,d=0,s=0;
        total_w++;
        good_w+=type;
        wcigar+= type == true ? "M":"S";
    }
    //step 3: add good windows to a vector
    //cout << "TOTALW="<<total_w<<" GOODW="<<good_w<<" "<<" BADW="<<total_w-good_w<<" CIGARW="<<wcigar<<" SW="<<result.startLocations[0]/w<<endl;
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
        /*cout << "TOTALW="<<total_w<<" GOODW="<<good_w<<" "<<" BADW="<<total_w-good_w<<" CIGARC="<<wcigar.substr(begin,(end-begin)+1)<<" SW="<<result.startLocations[0]/w
             <<" B="<<begin<<" END="<<end<<endl;*/
        for(auto k=begin; k<=end; k++)
                this->windows[localw[k].w].push_back(localw[k]);

    }

    return good_w;

}


void WConsensus::printAlignment(const char* query, const char* target, const unsigned char* alignment,  int alignmentLength, const int correctedaln,  const int position, const EdlibAlignMode modeCode) const {
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
    alignmentLength=correctedaln;
    for (int start = 0; start < alignmentLength; start += 100) {
        // target
        printf("T: ");
        int startTIdx;
        for (int j = start; j < start + 100 && j < alignmentLength; j++) {
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
        for (int j = start; j < start + 100 && j < alignmentLength; j++) {
            printf(alignment[j] == EDLIB_EDOP_MATCH ? "|" : " ");
        }
        printf("\n");
        // query
        printf("Q: ");
        int startQIdx = qIdx;
        for (int j = start; j < start + 100 && j < alignmentLength; j++) {
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

WConsensus::~WConsensus() {
    //we erase the containers
    /*for(auto w:this->windows)
         w.erase(w.begin(),w.end());*/
    this->windows.clear();
    this->windows.shrink_to_fit();
        //we erase the container
   // this->windows.erase(this->windows.begin(),this->windows.end());
    //this->coverage.erase(this->coverage.begin(),this->coverage.end());
    this->coverage.clear();
    this->coverage.shrink_to_fit();
}

string WConsensus::get_coverage_cns_string() {

    char c=(char)1;//there is only one long read
    string cons_d (this->coverage.size(), c);
    for(auto i=0; i<this->coverage.size(); ++i)
        cons_d[i]=(char)(this->coverage[i]);

    return cons_d;
}






