//
// Created by cuper on 2021/3/30.
//

#include "kcell.h"

kcell::kcell() {
    curk=-1;
    objID=-1;
    topk.clear();
    Stau.clear();
}

kcell::~kcell() {
    Stau.clear();
    vector<int>().swap(Stau);
    topk.clear();
    vector<int>().swap(topk);
}

void kcell::WriteToDisk(ofstream& Outfile) {
    Outfile.write((char*) &curk,sizeof(int));
    Outfile.write((char*) &objID,sizeof(int));
    r.WriteToDisk(Outfile);
    int size=topk.size();
    Outfile.write((char*) &size,sizeof(int));
    for (auto it=topk.begin();it!=topk.end();it++){
        Outfile.write((char*) &(*it),sizeof(int));
    }
    size=Stau.size();
    Outfile.write((char*) &size,sizeof(int));
    for (auto it=Stau.begin();it!=Stau.end();it++){
        Outfile.write((char*) &(*it),sizeof(int));
    }
}

void kcell::ReadFromDisk(ifstream& Infile) {
    Infile.read((char*) &curk,sizeof(int));
    Infile.read((char*) &objID,sizeof(int));
    r.ReadFromDisk(Infile);

    
    int size;
    Infile.read((char*) &size,sizeof(int));
    topk.clear();
    int* topk_buffer = new int [size];
    Infile.read((char*) topk_buffer,sizeof(int)*size);
    for (int i=0;i<size;i++){
        topk.push_back(topk_buffer[i]);
    }
    delete[] topk_buffer;

    Infile.read((char*) &size,sizeof(int));
    Stau.clear();
    int* Stau_buffer = new int [size];
    Infile.read((char*) Stau_buffer,sizeof(int)*size);
    for (int i=0;i<size;i++){// for very large k
        Stau.push_back(Stau_buffer[i]);
    }
    delete[] Stau_buffer;

    Get_HashValue();
}

void kcell::TobeRoot(vector<int> &candidates, int dim) {
    curk=0;
    objID=-1;
    topk.clear();
    for (int i=0;i!=candidates.size();i++) Stau.push_back(i);
    r.ToBeRoot(dim);
}

void kcell::Get_HashValue() {
    vector<int> tmp; tmp.clear();
    for (auto it=topk.begin();it!=topk.end();it++) tmp.push_back(*it);
    std::sort(tmp.begin(),tmp.end()); // make it in order to output the unique hash_value for a set
    std::size_t seed = 0;
    for (auto it=tmp.begin();it!=tmp.end();it++){
        boost::hash_combine(seed, *it);
    }

    boost::hash_combine(seed, objID); // objID will be computed twice, to make sure kcell shared the same hash value have the same objID
    hash_value=seed;
    return;
}