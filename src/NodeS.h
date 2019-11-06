//
// Created by Alex Digenova on 7/4/18.
//

#ifndef LIGER_NODES_H
#define LIGER_NODES_H


#include <utility>
#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>

using namespace std;

class NodeS {
private:
    uint32_t id;
    int ctg_id;
    bool repeat;
    bool cshort;
    int ctg_len;

public:


    NodeS(uint32_t id, int ctg, bool repeat, bool cshort,int ctg_len) : id(id),ctg_id(ctg), repeat(repeat),cshort(cshort),ctg_len(ctg_len){};


    //Getters and Setter
    uint32_t getId() const {
        return id;
    }

    void setId(uint32_t id) {
        NodeS::id = id;
    }


    int getCtg_id() const {
        return ctg_id;
    }

    void setCtg_id(int ctg_id) {
        NodeS::ctg_id = ctg_id;
    }
    bool is_repeat() const{
        return repeat;
    }

    bool is_short() const{
        return cshort;
    }

    int getnode_len() const {
        return ctg_len;
    }

    //copy a node
    NodeS(NodeS* n) {
        //we copy the node not the adj_list
        this->id=n->getId();
        this->ctg_id=n->getCtg_id();
        this->ctg_len=n->getnode_len();
        this->repeat=n->is_repeat();
        this->cshort=n->is_short();
    }
};


#endif //LIGER_NODES_H
