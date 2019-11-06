//
// Created by Alex Digenova on 2/20/18.
//

#include "MPLibrary.h"


//default construct
MPLibrary::MPLibrary(string file) {
    //we init the defaults variables
    mobs=10000;
    avg_insert_size=0;
    std_insert_size=0;
    rank=0;
    pout=0;
    this->file=file;
    //we infer the insert size of the library
    this->_infer_insert_size();
}


MPLibrary::MPLibrary(string file, int avg_ins, int std_ins ) {
    mobs=10000;
    avg_insert_size=avg_ins;
    std_insert_size=std_ins;
    rank=0;
    pout=0.1;
    this->file=file;
}


void MPLibrary::_infer_insert_size(){
    auto samreader=new SAMR(this->file);
    vector<int> obs_insert;
    int counter=0;
    while(!samreader->sam_has_reads()){
        auto fwd=samreader->get_next_read();
        //we are reading the header of the SAMfile
        if(fwd.pos==-1) {
            continue;
        }
        //auto fwd=parser.parse_entry(lfwd);
        auto rev=samreader->get_next_read();
        //we are in the same contig
        if(rev.contig.compare(fwd.contig) == 0){
            obs_insert.push_back(abs(rev.pos-fwd.pos));
            counter++;
        }
        //we have all the observations that we need
        if(counter >=mobs){
            break;
        }
    }
    //we close the current SAM file
    //readsam.close();
    //we call the destructor
    samreader->~SAMR();
    //we compute the avg insert size distance
    //cout << "Num obs: "<<counter << endl;
    //this->compute_lib_average(obs_insert);
    this->avg_insert_size=compute_average(obs_insert);
    //this->compute_lib_std(obs_insert);
    this->std_insert_size=int(avg_insert_size * 0.1);
}

void MPLibrary::load_library(Contig *a) {
    auto samreader=new SAMR(this->file);
    //boundaries for max and min distance among pairs
    float min = float(-this->std_insert_size * 2.5);
    float max = float(this->avg_insert_size+this->std_insert_size * 2.5);
    //cout << " Min " << min << " Max "<< max<<endl;
    //a tmp storage of links, we dump the content of this every 500 links are stored
    vector<linking> locallinks;
   // unordered_map<string, bool> edges2links; // store the relation between edges and links


    int idl=0;//and id for each stored link
    while(!samreader->sam_has_reads()){
        auto fwd=samreader->get_next_read();
        //we are reading the header of the SAMfile
        if(fwd.pos==-1) {
            continue;
        }
        auto rev=samreader->get_next_read();
        //we are in the same contig
        if(rev.contig.compare(fwd.contig) == 0){
            //we save the contig coverage of the shorter library for repeat identification
            if(this->rank ==0) {
                a->add_ctg_reads(a->ctg2id(fwd.contig), 2);
                //fwd and rev come from the same contig
                //long_reads[fwd.longread]+=1;
                //todo: improve this to be thread safe just save this to a local container and dump the content at the end of the execution
                //a->add_long_reads(a->ctg2id(fwd.contig), fwd.longread);
                //a->add_long_reads(a->ctg2id(fwd.contig), fwd.longread);
            }
            this->ctg2lr[a->ctg2id(fwd.contig)][fwd.lid]+=2;
        }else{
            //we have to create a link between contigs if possible
            //we have to determine if the link is ok
            vector<int> ori=this->get_orientation(&fwd,&rev);
            int rpos1=fwd.pos+fwd.len;
            int rpos2=rev.pos+rev.len;
            int distance=0;
            //contig orientation
            bool oc1=0;
            bool oc2=0;
            //computing the distance and ctg orientation for each putative link
            if(ori[0] == 0 and ori[1] == 0){ // reads  - - => ctg + -
                distance=this->avg_insert_size - (a->getcontig_length(fwd.contig) - fwd.pos) - rpos2;
                oc1=1; oc2=0; // + -
            }else{
                if(ori[0] == 1 and ori[1] == 1){ //reads + + => ctg - +
                    distance=this->avg_insert_size - (a->getcontig_length(rev.contig) - rev.pos) - rpos1;
                    oc1=0; oc2=1; // - +
                }else{
                    if(ori[0] == 0 and ori[1] == 1){ //reads + - => + +
                        distance=this->avg_insert_size - (a->getcontig_length(fwd.contig) - fwd.pos) - (a->getcontig_length(rev.contig) - rev.pos);
                        oc1=1; oc2=1; // + +
                    }else{
                        if(ori[0] == 1 and ori[1] == 0){ //reads - +  =>  - -
                            distance=this->avg_insert_size - rpos1 - rpos2;
                            oc1=0; oc2=0; // - -
                        }
                    }
                }
            }

            //we can determine if the distance is correct to include the edge in the
            if(min < distance and max > distance ){
                //due that is a valid pair we can increase the coverage of both cotigs by one
                if(this->rank ==0) {
                    a->add_ctg_reads(a->ctg2id(fwd.contig),1);
                    a->add_ctg_reads(a->ctg2id(rev.contig),1);
                    //todo: this is not thread safe unless under the current conf
                    //we add the long reads to the contigs to ask after for them
                    //a->add_long_reads(a->ctg2id(fwd.contig),fwd.longread);
                    //a->add_long_reads(a->ctg2id(rev.contig),rev.longread);
                }
                //this sctructure save the long reads to contigs
                this->ctg2lr[a->ctg2id(fwd.contig)][fwd.lid]++;
                this->ctg2lr[a->ctg2id(rev.contig)][rev.lid]++;

                    uint32_t c1,c2;
                    if(oc1){
                        //ctg1+="+";
                        c1=a->get_contig_header(fwd.contig);
                    }else{
                        //ctg1+="-";
                        c1=a->get_contig_tail(fwd.contig);
                    }
                    if(oc2){
                        //ctg2+="+";
                        c2=a->get_contig_header(rev.contig);
                    }else{
                        //ctg2+="-";
                        c2=a->get_contig_tail(rev.contig);
                    }

                    // we create a putative edge connecting both one of the boths ends of the contigs
                    bool cchange=false;
                    if(c1 > c2) {
                        auto t=c2;
                        c2=c1;
                        c1=t;
                        cchange=true;
                        //we have to change the pos also
                        auto tp=fwd.pos;
                        fwd.pos=rev.pos;
                        rev.pos=tp;
                    }

                    //todo: we have to check if is necesary to change the other varibles
                    //we create a link
                   //linking tmp(idl,a->get_contig_id(fwd.contig),a->get_contig_id(rev.contig),fwd.pos,rev.pos,oc1,oc2,bool(ori[1]),bool(ori[2]),distance,e,cchange,fwd.longread);
                   //linking tmp(idl,oc1,oc2,bool(ori[1]),bool(ori[2]),distance,e,cchange,fwd.longread);
                   linking tmp(idl,c1,c2,fwd.pos,rev.pos,oc1,oc2,bool(ori[1]),bool(ori[2]),distance,cchange,fwd.lid,fwd.gaps,fwd.gape);
                //string seq
                    //string el=join("__",vector<string>({e,fwd.longread}));
                    //if this variable is 0, we have to load the tmp seq
                    //if(edges2links.count(el) == 0){
                        //we load the seq, otherwise is not necesary
                        //todo:testing human genomes we dont load the seq for the moment
                        //tmp.addseq(fwd.gap);
                    //}
                    //edges2links[el]=true;//.push_back(idl);
                   locallinks.push_back(tmp);

                    idl++;//we increase the id of link
                //500 hundred is the size of the local buffer
                if(locallinks.size() > 500){
                    for(auto l:locallinks){
                        //we save the links to the current library of links
                        this->links.push_back(l);
                    }
                    locallinks.erase(locallinks.begin(),locallinks.end());
                }
            }
        }

    }

    //we recover the last links added to the local lib storage recover them
    if(!locallinks.empty()){
        for(auto l:locallinks){
            //we save the links to the current library of links
            this->links.push_back(l);
        }
        locallinks.erase(locallinks.begin(),locallinks.end());
    }
    //we clean the links
    //edges2links.erase(edges2links.begin(),edges2links.end());
    //we remove the samreader
    samreader->~SAMR();
}

//funtion to print the links of a given library
void MPLibrary::print_links(){

    for(auto l:this->links){
        /*if(l.has_seq()){
            *//*cout <<"OK "<<this->getRank()<<" "<<l.lid<<" "<<l.edge << " " << l.longread << " " << l.ctg1 << " " << l.ctg2 << " " << l.dist << " "
                 << l.getgapseq() << endl;*//*
            cout <<"OK "<<this->getRank()<<" "<<l.lid << " " << l.longread << " " << l.dist << " "
                 << l.getgapseq() << endl;
        }else{*/
            //cout <<"OK "<<this->getRank()<<" "<<l.lid<<" "<<l.edge << " " << l.longread << " " << l.ctg1 << " " << l.ctg2 << " " << l.dist <<" "<<  endl;
            cout <<"OK "<<this->getRank()<<" "<<l.lid << " " << l.lonread_id << " " << l.dist <<" "<<  endl;
        //}
    }

}


//this function return the orientation of pairs given a FR librarie
vector<int> MPLibrary::get_orientation(shortread *a, shortread *b){

    vector<int> ori(2,0);
    //FROM the SAM reader we got true if the read is mapped in the reverse strand
    //Here we swicth to false and true is used for foward strand
    //fwd read
    if(a->ori){
        ori[0]=0;//means -
    }else{
        ori[0]=1;//means +
    }
    //rev read
    if(b->ori){
        ori[1]=0;//means -
    }else{
        ori[1]=1;//means +
    }

    //all libs are in "FR" orientation
    //we adjust the orientation of the reverse mate
    ori[1]= 1 - ori[1];
    return ori;
}

vector<linking> MPLibrary::getLinks() {
    return this->links;
}


//to store the long reads that traverse the contig
void MPLibrary::load_reads2ctg(Contig *a) {
    for (auto lr: this->ctg2lr) {
        for(auto ll:this->ctg2lr[lr.first]){
            int counter=0;
            //cout << " "<<lr.first<<" "<<ll.first<<" "<<ll.second<<endl;
            a->add_long_reads(lr.first, ll.first,ll.second);
        }
    }
    //we delete the local container
    this->ctg2lr.erase(this->ctg2lr.begin(),this->ctg2lr.end());
}
