//
// Created by Alex Digenova on 8/6/18.
//

#include "BSeqDB.h"



BSeqDB::BSeqDB(string fasta,simplehashu &sorder, simplehashu &lcache) {

    //string filename=argv[1];
    prefix=to_string(getpid());
    std::ifstream input(fasta);
    if(!input.good()){
        std::cerr << "Error opening '"<<fasta<<"'. Bailing out." << std::endl;
        exit(1);
    }
    ushort currentfile=0;
    filebuf wb;
    wb.open ("seqtmp."+prefix+"."+to_string(currentfile)+".bin",std::ios::out);
    ostream os(&wb);
    uint32_t id=0;
    std::string line, name, content;
    vector<bseq> bbseq;


    int size=0;

    unordered_map<uint32_t,uint64_t> maxseek;

    while(getline( input, line ).good() ){
        if( line.empty() || line[0] == '>' ){ // Identifier marker
            if( !name.empty()){ // Print out what we read from the last entry

                bseq tmp;
                tmp.id= static_cast<uint32_t>(atoi(name.c_str()));
                //if the sequence is defined in the needed sequences
                if(sorder.count(tmp.id)) {
                    tmp.seq = content;
                    bbseq.push_back(tmp);
                    size+=content.length();
                    //we need to cache this read
                    if(lcache.count(tmp.id)){
                        auto sc=new DnaBitset(tmp.seq.c_str(),tmp.seq.length());
                        //we store the cache reads
                        lrcache[tmp.id]=sc;
                    }
                    totalseqs++;
                    name.clear();
                    if (size > max_file_size) {
                        //we sort the partition
                        sort(bbseq.begin(), bbseq.end(),
                             [&sorder](const bseq a, const bseq b) { return sorder[a.id] < sorder[b.id]; });
                        for (auto s:bbseq) {
                            maxseek[currentfile] = static_cast<uint64_t >(os.tellp());
                            Storeseq(s.id, s.seq, static_cast<uint32_t>(s.seq.length()), os);
                            id++;
                        }
                        wb.close();
                        currentfile++;
                        wb.open("seqtmp." + prefix + "." + to_string(currentfile) + ".bin", std::ios::out);
                        size = 0;
                        bbseq.erase(bbseq.begin(), bbseq.end());
                    }
                }
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
        bseq tmp;
        tmp.id= static_cast<uint32_t>(atoi(name.c_str()));
        //if the sequence is used in the store
        if(sorder.count(tmp.id)) {
            tmp.seq = content;
            bbseq.push_back(tmp);
            totalseqs++;
            //we need to cache this read
            if(lcache.count(tmp.id)){
                auto sc=new DnaBitset(tmp.seq.c_str(),tmp.seq.length());
                //we store the cache reads
                lrcache[tmp.id]=sc;
            }
        }
    }

    //we sort and store the last part of the buffer
    sort(bbseq.begin(),bbseq.end(),[&sorder](const bseq a,const bseq b){return sorder[a.id] < sorder[b.id];});
    for(auto s:bbseq){
        maxseek[currentfile]=static_cast<uint64_t>(os.tellp());
        Storeseq(s.id, s.seq, static_cast<uint32_t>(s.seq.length()), os);
        id++;
    }
    //we clean the buffer and we close the files
    bbseq.erase(bbseq.begin(),bbseq.end());
    wb.close();
    input.close();
    //we prepare the readers of the files, for that we open the files
    numberfiles=currentfile+1;

    //we merge the tmp files using a priority queue
    auto cmp = [&sorder](const bseq a,const bseq b){return sorder[a.id] > sorder[b.id];};
    priority_queue<bseq, vector<bseq>, decltype(cmp) > q2(cmp);

    auto *rb = new filebuf[numberfiles];
    uint32_t seqid=0;
    //we read one seq from the tmp files
    for (int i = 0; i < numberfiles ; i++) {
            rb[i].open("seqtmp."+prefix+"."+to_string(i)+".bin",std::ios::in);
            istream is(&rb[i]);
            is.seekg(0);
            auto s = ReadSeq(is,seqid);
            bseq tmp;
            tmp.id= static_cast<uint32_t>(seqid);
            tmp.seq=s;
            tmp.fid= static_cast<ushort>(i);
            q2.push(tmp);
            //cout << " F1 "<<tmp.id <<" "<<tmp.fid<<" "<<sorder[tmp.id]<<" "<<tmp.seq<<endl;

    }
    //while the queue is not empty we read from the parts
    //cout << "Merging sorting files"<<endl;
    //we create the output file
    //index of the database
    index.open ("index."+prefix+".sorted.bin",std::ios::out);
    ostream oss(&index);
    //int n=0;
    uint32_t current_chunk=0;
    pair<uint64_t ,uint64_t > p(0,0);

    while(!q2.empty()){
        auto lowest = q2.top();
        //cout << lowest.id <<" "<<lowest.fid<<" "<<sorder[lowest.id]<<" "<<lowest.seq<<endl;
        seq2pos[lowest.id]= static_cast<uint64_t>(oss.tellp());
        seq2chunk[lowest.id]=current_chunk;
        Storeseq(lowest.id, lowest.seq, static_cast<uint32_t>(lowest.seq.length()), oss);
        auto tmp_cs=oss.tellp();
            tmp_cs-=p.first;
        if(tmp_cs >= chunk_size){
            p.second= static_cast<uint64_t>(oss.tellp());
            chunks.push_back(p);
            p.first=static_cast<uint64_t>(oss.tellp());
            current_chunk++;
        }
        // remove this record from the queue
        q2.pop();
        // add the next line from the lowest stream (above) to the queue
        // as long as it's not EOF.
        istream is(&rb[lowest.fid]);
        //is.seekg(postfile[lowest.fid]);
        //we read a new record an push it to the end of the queue
        if(is.tellg() <= maxseek[lowest.fid]){
            auto s = ReadSeq(is,seqid);
            bseq tmp;
            tmp.id= static_cast<uint32_t>(seqid);
            tmp.seq=s;
            tmp.fid= lowest.fid;
            q2.push(tmp);
        }else{
            //we delete the tmp file
            string tmpfile="seqtmp."+prefix+"."+to_string(lowest.fid)+".bin";
            rb[lowest.fid].close();
            remove(tmpfile.c_str());
        }
    }
    //we save the last chunk
    p.second=static_cast<uint64_t>(oss.tellp());
    chunks.push_back(p);

    //we end the reading and sorting of the file;
    index.close();
    //freezin memory
    delete[] rb;
    //auto max_chunk_size=0;
    //auto j=0;
    for(auto i:chunks){
       // cout <<"chunks "<<j<<" "<<i.first<<" "<<i.second<<endl;
       // j++;
        auto sizec=i.second - i.first;
        if(sizec > max_chunk_size)
            max_chunk_size=sizec;
    }
    cout << "Allocating max chunk size"<<max_chunk_size<<endl;
    //we allocate the maximum buffer size
    bchunk=new char[max_chunk_size+1];

}

//store a sequence in a binary file
void BSeqDB::Storeseq(uint32_t val,string seq, uint32_t len,ostream& file){

    char buf[4];
    memcpy(buf, &val, 4);
    //we save the readID
    file.write( buf, 4);

    auto eseq=encodeseq(seq);
    memcpy(buf, &len, 4);
    //we write the seq length
    file.write( buf, 4);
    //bytes used
    size_t dna_bytes = (len / 4) + (len % 4 != 0);
    memcpy(buf, &dna_bytes, 4);
    //we write the bytes used to store seq
    file.write( buf, 4);
    //we write the encoded seq
    file.write( eseq, dna_bytes);
    //file.write( seq.c_str(), len);
    //we clear the seqstring object
    //todo:test if this is usefull or not
    seq.clear();
    delete[] eseq;
}


char* BSeqDB::encodeseq(string seq) {

    auto dna_len=seq.length();
    /* number of bytes necessary to store dna_str as a bitset */
    size_t dna_bytes = (dna_len / 4) + (dna_len % 4 != 0);

    auto m_data = new char[dna_bytes];
    //we set the char to caracters
    std::memset(m_data, 0, dna_bytes);

    /* for each base of the DNA sequence */
    for (size_t i = 0; i < dna_len; ++i)
    {
        uint8_t shift = 6 - 2 * (i % 4);

        switch (seq[i])
        {
            case 'A':
                m_data[i / 4] |= BASE_A << shift;
                break;
            case 'C':
                m_data[i / 4] |= BASE_C << shift;
                break;
            case 'G':
                m_data[i / 4] |= BASE_G << shift;
                break;
            case 'T':
                m_data[i / 4] |= BASE_T << shift;
                break;
            default:
                throw std::invalid_argument("invalid DNA base found in DnaBitset class");
        }

        shift = (shift == 0) ? 6 : shift - 2;
    }

    return m_data;
}

char* BSeqDB::decodeseq(char* seq, int len) {
    auto m_len=len;
    char* dna_str = new char[m_len + 1];
    std::memset(dna_str, 0, m_len+1);

    /* for each base of the DNA sequence */
    for (size_t i = 0; i < m_len; ++i)
    {
        uint8_t shift = 6 - 2 * (i % 4);
        uint8_t mask = BASE_MASK << shift;

        /* get the i-th DNA base */
        uint8_t base = (seq[i / 4] & mask) >> shift;

        switch (base)
        {
            case BASE_A:
                dna_str[i] = 'A';
                break;
            case BASE_C:
                dna_str[i] = 'C';
                break;
            case BASE_G:
                dna_str[i] = 'G';
                break;
            case BASE_T:
                dna_str[i] = 'T';
                break;
            default:
                throw std::runtime_error("invalid DNA base");
        }
        //cout << dna_str[i]<<" "<<i<<endl;
    }

    dna_str[m_len] = '\0';
    return dna_str;

}






//read a sequence from a binary file, expensive method
string BSeqDB::ReadSeq(istream& file, uint32_t& id)
{   //read the seqid
    char buf[4];
    file.read( buf, 4 );
    //int32_t val;
    memcpy(&id, buf, 4);
    //read the seq length
    char len[4];
    file.read( len, 4 );
    int32_t slen;
    memcpy(&slen, len, 4);

    //read the bytes used to store the seq
    //char len[4];
    file.read( len, 4 );
    int32_t bytes_len;
    memcpy(&bytes_len, len, 4);
    //read the sequence
    auto bseq = new char[bytes_len];
    file.read( bseq, bytes_len);
    auto dseq=decodeseq(bseq,slen);
    string tmp(dseq);
    //we delete the sequences
    delete[] bseq;
    delete[] dseq;
    return tmp;
}

void BSeqDB::testRead(unordered_map<uint32_t,uint32_t> & lorder) {


    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator
    std::uniform_int_distribution<uint32_t> distr(0, totalseqs-1); // define the range
    //exit(0);
    cout << "performing TESTREAD"<<endl;

    clock_t begin = clock();
    int iter=1;
    //if(totalseqs > )
    //auto iter=int(1000000/totalseqs);
    int round=0,q=0;
    //int seqid=0;
    //we open the index file
    index.open("index."+prefix+".sorted.bin",std::ios::in);


    while(round < iter) {
        //for (uint32_t i = 0; i < totalseqs; i++) {
        for(auto l:lorder){
                //auto s=getSeq(i);
                auto s=ChunkGetSeq(l.first);
                /*if(s.seq.length() > 0)
                cout << s.id <<" "<<s.seq<<endl;
                */
                q++;
        }
        round++;
    }
    index.close();
        clock_t end = clock();
   double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "TotalSeqs : "<<totalseqs<<". Time spent in "<< q<<" queries  from chunk" << elapsed_secs <<endl;

}

void BSeqDB::DBdump() {
    uint32_t seqid=0;
    //we open the master file
    index.open("index."+prefix+".sorted.bin",std::ios::in);
    istream is(&index);
    //we sort by id
    vector<uint32_t > tmp;
    for(auto i:seq2pos){
        tmp.push_back(i.first);
    }
    //we sort according to the pos in the file that should match the given order
    sort(tmp.begin(),tmp.end(),[this](const int32_t a, const int32_t b){return seq2pos[a] < seq2pos[b];});
    for(auto i:tmp){
            is.seekg(seq2pos[i]);
            auto s = ReadSeq(is, seqid);
            //cout << i << " "<<seq2pos[i]<<" "<< seqid << " " << s.length() << " " << s << endl;
        cout <<"LRDB: " <<i<< " "<<seq2pos[i]<<" "<<seq2chunk[i]<<" "<<chunks[seq2chunk[i]].first<<" "<<chunks[seq2chunk[i]].second<<" "<<seqid << " " << s.length() << " " << s << endl;
    }
    index.close();

}

BSeqDB::~BSeqDB() {
    //deleting the index file
    string indexfile="index."+prefix+".sorted.bin";
    remove(indexfile.c_str());
    //deleting the bbuffer
    //if(loaded_chunk>=0)
    delete[] bchunk;
    //deleting containers
    seq2pos.erase(seq2pos.begin(),seq2pos.end());
    seq2chunk.erase(seq2chunk.begin(),seq2chunk.end());
    chunks.erase(chunks.begin(),chunks.end());

}

bseq BSeqDB::getSeq(uint32_t seqid) {
    uint32_t ss=0;
    bseq tmp;
    if(seq2pos.count(seqid)) {
        istream is(&index);
        is.seekg(seq2pos[seqid]);
        auto s = ReadSeq(is, ss);
        tmp.seq=move(s);
        tmp.id=ss;
    }else{
        cerr << "sequence not stored in bseqDB :"<<seqid<<endl;
    }
    return tmp;
}

void BSeqDB::load_chunk(uint32_t cid) {
    //the chunk is in memory
    if(cid == loaded_chunk){
        return;
    }

    //we delete the previous chunk only if the chunk was loaded before cid > 0
    if(cid > 0){
        //we delete the current chunk
     //   delete[]  bchunk;
    }


    auto size=chunks[cid].second - chunks[cid].first;
    //size++;
    //we create the new size
     //bchunk=new char[size];
    //we clear the buffer
     //memcpy(bchunk,"\0",max_chunk_size+1);
     memset(bchunk,'0',max_chunk_size+1);
    index.open("index."+prefix+".sorted.bin",std::ios::in);
    istream is(&index);
    //we go to the chunk start
    is.seekg(chunks[cid].first);
    //we read the chunk lenght
    is.read(bchunk,(int)size);
    if (is){
        cout << "all characters read successfully. "<< is.gcount()<<" "<<cid<<endl;
    }else{
        cout << "error: only " << is.gcount() << " could be read\n";
    }
    //string tmp(bchunk);
   // cout << "size of chunk "<<size<<" "<<strlen(bchunk)<<" "<<tmp.length    ()<<endl;
    index.close();
    loaded_chunk=cid;
}


bseq BSeqDB::ChunkGetSeq(uint32_t seqid) {
    //uint32_t ss=0;
    bseq tmp;
    tmp.id=0;
    tmp.seq="";
    //seq not stored in the database
    if(seq2chunk.count(seqid) == 0){
        cout << "Asking for a sequence not stored in the database"<<endl;
        return tmp;
    }

    //read could be stored in the cache memory
    if(lrcache.count(seqid)){
        auto s=lrcache[seqid];
        tmp.seq=s->to_string();
        tmp.id=seqid;
        return tmp;
    }

    //magic pointer to chunk data
    char* pos;
    //we load the chunk if necessary
    //cout <<"SEQ2CHUNK: "<<seq2chunk[seqid]<<endl;
    load_chunk(seq2chunk[seqid]);
    //cout << "Chunk loaded"<<endl;
    //now we can search the sequence in the bchunk
    //we get the seq start
    auto bstart=seq2pos[seqid]-chunks[seq2chunk[seqid]].first;
    //cout <<seqid<<" "<<bstart<<" "<<seq2pos[seqid]<<" "<<seq2chunk[seqid]<<" "<<chunks[seq2chunk[seqid]].first<<" "<<chunks[seq2chunk[seqid]].second<<endl;
    //bstart is the fisrt bit of the seq
    uint32_t id=0;
    pos=bchunk+bstart;
    memcpy(&id, pos, 4);
    uint32_t len=0;
    pos=bchunk+bstart+4;
    memcpy(&len,pos, 4);
    uint32_t blen=0;
    pos=bchunk+bstart+8;
    memcpy(&blen,pos, 4);
    //auto tseq = new char[blen];
   // memset (tseq,'a',blen);
    //cout << id<<" "<<seqid<<" "<<len<<" "<<blen<<" "<<blen*4<<endl;
    auto seqbits=new char[blen];
    pos=bchunk+bstart+12;
    memcpy(seqbits,pos,blen);
    auto ddseq=decodeseq(seqbits,len);
    string tmpseq(ddseq);
     //we fill the string seq
     tmp.seq=move(tmpseq);
     tmp.id=id;
     //we delete the sequences
     delete[] seqbits;
     delete[] ddseq;
     tmpseq.clear();

    return tmp;
}

uint BSeqDB::getMax_file_size() const {
    return max_file_size;
}

uint BSeqDB::getTotalseqs() const {
    return totalseqs;
}

uint BSeqDB::getNumberfiles() const {
    return numberfiles;
}

void BSeqDB::DumpCache() {
    cout << "DUMPING Cache reads"<<endl;
    for(auto s:lrcache){
        auto seq=lrcache[s.first];
        cout << s.first<<" "<<seq->to_string()<<" "<<seq->to_string().length()<<endl;
    }


}

uint BSeqDB::getNumberLRCahed() const {
    return lrcache.size();
}

