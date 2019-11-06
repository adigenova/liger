//
// Created by Alex Di genova on 3/6/18.
//

#include "Contig.h"

//Construct contigs from FASTA file
Contig::Contig(string filename){
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
                tmp.header=ori;
                ori++;
                tmp.tail=ori;
                ori++;
                this->contigs.push_back(tmp);
                this->string2id[tmp.namectg]=id;
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
        tmp.header=ori;
        ori++;
        tmp.tail=ori;
        ori++;
        this->contigs.push_back(tmp);
        this->string2id[tmp.namectg]=id;
    }
}

int Contig::ctg2id(string name) {
    return this->string2id[name];
}

ctg Contig::getcontig(string name) {
    return this->contigs[string2id[name]];
}

ctg Contig::getcontig(int id) {
    return this->contigs[id];
}


int Contig::getcontig_length(string name) {
    return this->contigs[string2id[name]].length;
}

int Contig::getcontig_length(int id) {
    return this->contigs[id].length;
}

string Contig::getcontig_seq(string name) {
    return this->contigs[string2id[name]].seq;
}

string Contig::getcontig_seq(int id) {
    return this->contigs[id].seq;
}

void Contig::print_contigs() {
    for(auto &c:this->contigs){
        cout << c.ctgid<< " " <<c.namectg <<" "<< c.length<<" "<<c.seq<<endl;
    }
}

void Contig::print_contigs_file(string file) {
    ofstream rcon;
    //todo: change the contig name file
    rcon.open(file+".contigs");
    rcon << "#Ctg_id Name length"<<endl;
    for(auto &c:this->contigs){
        rcon << c.ctgid<< " " <<c.namectg <<" "<< c.length<<" "<<c.seq<<endl;
    }
}

vector<ctg> Contig::get_all_contigs() {
    return this->contigs;
}

int Contig::get_contig_id(string name) {
    return this->string2id[name];
}


float Contig::getAvg_ctg_cov() const {
    return avg_ctg_cov;
}

void Contig::setAvg_ctg_cov(float avg_ctg_cov) {
    Contig::avg_ctg_cov = avg_ctg_cov;
}

float Contig::getStd_ctg_cov() const {
    return std_ctg_cov;
}

void Contig::setStd_ctg_cov(float std_ctg_cov) {
    Contig::std_ctg_cov = std_ctg_cov;
}

void Contig::mark_repeats() {
    vector<float> cov;
    for (auto con:this->contigs) {
        float c;
        c = (con.num_reads * 200 / con.length);
            if(c>0) {
                cov.push_back(c);
            }
       this->set_ctg_coverage(con.ctgid,c);
    }
    //AVG coverage and STD
    avg_ctg_cov=compute_average(cov);
    std_ctg_cov=compute_std(cov,this->avg_ctg_cov);
    //todo:improve this when the std_ctg_cov is to large
    //this could be inferrend from the illuumina asssembler
    //float max = avg_ctg_cov + 4*std_ctg_cov;
    auto max = 4.5*avg_ctg_cov;
    int repeats=0;
    cout << "Max allowed coverage "<<max<<endl;
    for (auto con:this->contigs) {
        //cout << "NOREPEAT CTG:"<<con.ctgid<<" "<<con.length<<" "<<con.repeat<<" "<<con.coverage<<endl;
        if(con.coverage > max){
            this->set_repeat(con.ctgid);
            //cout << "REPEAT CTG:"<<con.ctgid<<" "<<con.length<<" "<<con.repeat<<" "<<con.coverage<<endl;
            repeats++;
        }else{
            //cout << "NONREPEAT CTG:"<<con.ctgid<<" "<<con.length<<" "<<con.repeat<<" "<<con.coverage<<endl;
        }
    }
    cout << "AVG "<<avg_ctg_cov <<" STD "<< std_ctg_cov <<" CTG REP "<<repeats<<endl;
}

//method to load the short-read contig coverage from a file [ctgname coverage]
void Contig::mark_repeats(string filename) {
    std::ifstream input(filename);
    if(!input.good()){
        std::cerr << "Error opening '"<<filename<<"'. Bailing out." << std::endl;
        exit(1);
    }
    string line;
    vector<float> cov;
    while(getline( input, line ).good() ){
        auto result=splits(" ",line);
        assert(result.size() == 2);
        //we save the ctg coverage as a float
        auto c=atof(result[1].c_str());
        //we set the coverage of the contig from illumina-reads
        this->set_ctg_coverage(this->get_contig_id(result[0]),c);
        //we save the coverage in tmp vector to compute STD and AVG
        cov.push_back(c);
    }
    //we discard 10% oulayers from boths tails
    sort(cov.begin(),cov.end());
    int a=int(cov.size()*0.05);
    vector<float> ocov(cov.begin()+a,cov.end()-a);
    //AVG coverage and STD
    this->avg_ctg_cov=compute_average(ocov);
    this->std_ctg_cov=compute_std(ocov,this->avg_ctg_cov);

    /*if(this->std_ctg_cov == 1){
        this->std_ctg_cov++;//we use and std min of 2
    }*/
    //I want the upper coverage
    float max = float(avg_ctg_cov+2.5*std_ctg_cov+0.5);
    int repeats=0;
    cout << "Max allowed coverage "<<max<<endl;
    for (auto con:this->contigs) {
        if(con.coverage > max){
            this->set_repeat(con.ctgid);
            //cout << "REPEAT CTG:"<<con.ctgid<<" "<<con.namectg<<" "<<con.length<<" "<<con.repeat<<" "<<con.coverage<<endl;
            repeats++;
        }else{
            //cout << "NONREPEAT CTG:"<<con.ctgid<<" "<<con.namectg<<" "<<con.length<<" "<<con.repeat<<" "<<con.coverage<<endl;
        }
    }
    cout << "AVG "<<avg_ctg_cov <<" STD "<< std_ctg_cov <<" CTG REP "<<repeats<<endl;
    //we clear the local vector used:
    cov.erase(cov.begin(),cov.end());
    ocov.erase(ocov.begin(),ocov.end());
}

//This procedure is similar to the a-statistic of myers but much more simple
void Contig::mark_repeats_astat(string filename, float rf) {
    std::ifstream input(filename);
    if(!input.good()){
        std::cerr << "Error opening '"<<filename<<"'. Bailing out." << std::endl;
        exit(1);
    }
    string line;
    vector< pair<int,float> > cov;
    while(getline( input, line ).good() ){
        auto result=splits(" ",line);
        assert(result.size() == 2);
        //we save the ctg coverage as a float
        auto c=atof(result[1].c_str());
        //we set the coverage of the contig from illumina-reads
        this->set_ctg_coverage(this->get_contig_id(result[0]),c);
        //we save the coverage in tmp vector to compute STD and AVG
        pair<int,float> tmp(this->getcontig_length(this->get_contig_id(result[0])),c);
        //cov.push_back(c);
        cov.push_back(tmp);
    }
    //we use the longest contigs to stimate the average contig  coverage, largest contigs are likely unique
    //This procedure is similar to the a-statistic of myers but much more simple.
    sort(cov.begin(),cov.end(),[](const pair<int,float> a, const  pair<int,float> b){return a.first > b.first;});
    //sort(cov.begin(),cov.end());
    int a=int(cov.size()*0.1);
    vector<float> ocov;
    for(auto j=0; j< a; j++){
            ocov.push_back(cov[j].second);
    }
    //AVG coverage and STD
    this->avg_ctg_cov=compute_average(ocov);
    this->std_ctg_cov=compute_std(ocov,this->avg_ctg_cov);
    //I want the upper coverage
    float max = static_cast<float>(avg_ctg_cov * rf);// repeat any contig with coverage longer 1.5 times the average
    int repeats=0;
    cout << "Max allowed coverage "<<max<<endl;
    for (auto con:this->contigs) {
        if(con.coverage > max){
            this->set_repeat(con.ctgid);
            //cout << "REPEAT CTG:"<<con.ctgid<<" "<<con.namectg<<" "<<con.length<<" "<<con.repeat<<" "<<con.coverage<<endl;
            repeats++;
        }else{
            //cout << "NONREPEAT CTG:"<<con.ctgid<<" "<<con.namectg<<" "<<con.length<<" "<<con.repeat<<" "<<con.coverage<<endl;
        }
    }
    cout << "AVG:"<<avg_ctg_cov <<" STD:"<<std_ctg_cov<<" Limit: "<<max<<" #REP(ctg):"<<repeats<<endl;
    //we clear the local vector used:
    cov.erase(cov.begin(),cov.end());
    ocov.erase(ocov.begin(),ocov.end());
}

uint32_t Contig::get_contig_header(string name) {
    return this->contigs[this->string2id[name]].header;
}

uint32_t Contig::get_contig_tail(string name) {
    return this->contigs[this->string2id[name]].tail;
}

int Contig::check_long_read(int id, const int lr) {
    //return this->contigs[id].long_reads.count(lr) > 0;
    if(this->contigs[id].long_reads.count(lr)){
        return this->contigs[id].long_reads[lr];
    }else {
        return 0;
    }
}

