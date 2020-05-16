//
// Created by Alex Di Genova on 25/07/2016.
//

#include "SAMR.h"


SAMR::~SAMR(){
//we delete the buffer and close the opened file
    delete[] buf;
    close(fd);

}
//default constructor
SAMR::SAMR(string file){

    fd = open(file.c_str(), O_RDONLY);
    if(fd == -1)
        handle_error("open");
    // obtain file size
    struct stat sb;
    if (fstat(fd, &sb) == -1)
        handle_error("fstat");
    //if file size is shorter than the BUFFER size
    if(sb.st_size < BUFFER_SIZE){
        //std::cout << sb.st_size<<endl;
        BUFFER_SIZE=sb.st_size-1;
        //cout << "Changing the buffer size to file size : "<< sb.st_size<<" "<<BUFFER_SIZE<<endl;
    }
    //we create the buffer
    buf = new char[BUFFER_SIZE+1];
    //we read an initial portion of the file to fill the buf
    bytes_read = read(fd, buf, BUFFER_SIZE);
    //we dont reads lines from the file
    if(bytes_read == (size_t)-1)
        handle_error("read failed");
    //we check that the file is not empty
    if (!bytes_read){
        //break;
        //there is no read in the file
        is_empty=true;
    }
    //we adjust the pointers to the start of the buffer
    init=buf;
    pos=buf;
    consumed=0;

}

//read data from the file
void SAMR::load_buffer(){
    bytes_read = read(fd, buf, BUFFER_SIZE);
    //we dont reads lines from the file
    if(bytes_read == (size_t)-1)
        handle_error("read failed");
    //we check that the file is not empty
    if (!bytes_read){
        //there is no read in the file
        is_empty=true;
        //means that we read data from the file
    }else{
        is_empty=false;
    }
    //we adjust the pointers to the start of the buffer
    init=buf;
    pos=buf;
    //and the consumed bases also
    consumed=0;

}

shortread SAMR::get_next_read(){
    shortread hit;//we create a hit
    //we look for the newline
    pos = (char*) memchr(pos, '\n', (buf + bytes_read) - pos);
    //we found the "\n" in the characters
    if(pos != NULL){
        //we copy the new line to the buffer
        memcpy(liner,init,pos - init);
        //we add the termination character
        liner[pos - init]='\0';
        //we count the consumed based
        consumed+=pos-init;
        //we create the string seq
        string tmp(liner);
        init=pos+1;
        pos++;
        if(tmp.find("@"))
            hit=parse_entry(tmp);

    }else{
        //we save the remainder of the buffer and we load the buffer again
        if(bytes_read - consumed > 0){
            //there is a remainig in the buffer
            memcpy(remain,init,buf+bytes_read - init);
            remain[buf+bytes_read - init]='\0';
            //cout <<"Remain again "<<bytes_read<<" "<<consumed<<" "<<remain<<endl;
            remainb=true;
        }
        //we load the buffer for the next read
        load_buffer();
        if(remainb && !sam_has_reads()){
            //we complete the current hit
            pos = (char*) memchr(pos, '\n', (buf + bytes_read) - pos);
            memcpy(liner,init,pos - init);
            //we add the termination character
            liner[pos - init]='\0';
            //we count the consumed based
            consumed+=pos-init;
            //we create the string seq
            string tmp(liner);
            //init=pos+1;
            init=pos+1;
            //we add the part of the remain sequence
            //if(remainb){
            string tmp2(remain);
            tmp=tmp2+tmp;
            remainb=false;
            pos++;
            //cout << tmp<<endl;
            if(tmp.find("@"))
                hit=parse_entry(tmp);
        }

    }
//we return the seq anyway
    return hit;

}

//printing errors and other things
void SAMR::handle_error(const char* msg) {
    perror(msg);
    exit(255);
}


//faster line parser with iterators
shortread SAMR::parse_entry(string original){
    shortread rp;
    //it split at 115354560/2m56.626 =~651,720 per sec 1.7x faster
    //if we break the while it took 115354560/2m27.061s=~2X faster
    char separator='\t';
    std::vector<std::string> results;
    std::string::const_iterator start = original.begin();
    std::string::const_iterator end = original.end();
    std::string::const_iterator next = std::find( start, end, separator);
    uint32_t counter=0;
    while ( next != end ) {
        results.push_back( std::string( start, next ) );
        start = next + 1;
        next = std::find( start, end, separator );
        counter++;
        //currently we are using the first 10 columms of the SAM
        if(counter > 10){
            break;
        }
    }

    //I need an assert to check the the SAM was well parsedl
    assert(results.size()>=6);
    //parsing each entry in the SAM
    rp.name=results[0];
    rp.contig=results[2];
    //read orientation for single and paired reads aligments
    rp.ori= (atoi(results[1].c_str()) & 16) ? false : true;
    //rp.len=int(results[9].length());
    rp.len=atoi(results[5].c_str());
    rp.pos=atoi(results[3].c_str());
    rp.lid=atoi(results[0].substr(0, results[0].find("_FG_")).c_str());
    //rp.gap=results[19].substr(5,results[19].length()-5);
    //cout << results[0]<<endl;
    auto lssp=results[0].substr(results[0].find("_FG_")+4,results[0].size());
    //cout << lssp<<endl;
    auto ls=atoi(lssp.substr(0,lssp.find("_")).c_str());
    auto sp=atoi(lssp.substr(lssp.find("_")+1,lssp.size()).c_str());
    //longread names are
    //14660869_FG_6000_13100
    //(lid)_FG_(ls)_(sp)
    //gap=[sp+200,sp+ls-200]
    //gap=[13300,18900]
    rp.gaps=sp+200;
    rp.gape=sp+ls-200;

    results.erase(results.begin(),results.end());
    return rp;
}
