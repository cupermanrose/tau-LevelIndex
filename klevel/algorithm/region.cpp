

#include "region.h"

region::region() {
    H.clear();
    V.clear();
    innerPoint.clear();
}

region::~region() {
    H.clear();
    vector<halfspace>().swap(H);
    V.clear();
    vector<vector<float>>().swap(V);
    innerPoint.clear();
    vector<float>().swap(innerPoint);
}

void region::WriteToDisk(ofstream &Outfile, bool OutHS) {
    int size;
    if (OutHS){
        size=H.size();
        Outfile.write((char*) &size, sizeof(int));
        for(auto it=H.begin();it!=H.end();it++){
            int w_size=it->w.size();
            Outfile.write((char*) &w_size,sizeof(int));
            for (auto itt=it->w.begin();itt!=it->w.end();itt++){
                Outfile.write((char*) &(*itt),sizeof(float));
            }
            Outfile.write((char*) &(it->side), sizeof(bool));
        }
    }
    else {// halfspace can be computed from Tk and Ct set
        size=0;
        Outfile.write((char*) &size, sizeof(int));
    }
    size=V.size();
    Outfile.write((char*)&size, sizeof(int));
    for (auto it=V.begin();it!=V.end();it++){
        int w_size=it->size();
        Outfile.write((char*) &w_size,sizeof(int));
        for (auto itt=it->begin();itt!=it->end();itt++){
            Outfile.write((char*) &(*itt),sizeof(float));
        }
    }
    size=innerPoint.size();
    Outfile.write((char*)&size,sizeof(int));
    for (auto it=innerPoint.begin();it!=innerPoint.end();it++){
        Outfile.write((char*) &(*it),sizeof(float));
    }
}

void region::ReadFromDisk(ifstream &Infile) {
    // Read H
    int size;
    Infile.read((char*) &size, sizeof(int));
    H.clear();
    for (int i=0;i<size;i++){
        int w_size;
        Infile.read((char*) &w_size, sizeof(int));
        halfspace HS;
        HS.w.clear();

        float* HS_buffer = new float [w_size];
        Infile.read((char*) HS_buffer,sizeof(float)*w_size);
        for (int j=0;j<w_size;j++) HS.w.push_back(HS_buffer[j]);
        delete[] HS_buffer;
        /*for (int j=0;j<w_size;j++){
            float w_value;
            Infile.read((char*) &w_value,sizeof(float));
            HS.w.push_back(w_value);
        }*/
        Infile.read((char*) &(HS.side),sizeof(bool));
        H.emplace_back(HS);
    }
    // Read V
    Infile.read((char*)&size, sizeof(int));
    V.clear();
    for (int i=0;i<size;i++){
        int w_size;
        Infile.read((char*) &w_size, sizeof(int));
        vector<float> point;
        point.clear();

        float* V_buffer = new float [w_size];
        Infile.read((char*) V_buffer,sizeof(float)*w_size);
        for (int j=0;j<w_size;j++) point.push_back(V_buffer[j]);
        delete[] V_buffer;

        /*for (int j=0;j<w_size;j++){
            float w_value;
            Infile.read((char*) &w_value,sizeof(float));
            point.push_back(w_value);
        }*/
        V.emplace_back(point);
    }
    // Read innerPoint
    Infile.read((char*) &size, sizeof(int));
    innerPoint.clear();

    float* IP_buffer = new float [size];
    Infile.read((char*) IP_buffer,sizeof(float)*size);
    for (int i=0;i<size;i++) innerPoint.push_back(IP_buffer[i]);
    delete[] IP_buffer;
    /*for (int i=0;i<size;i++){
        float w_value;
        Infile.read((char*) &w_value,sizeof(float));
        innerPoint.push_back(w_value);
    }*/
}

void region::ToBeRoot(int dim) {
    H.clear();
    V.clear();
    //Generate vertices for the whole space
    vector<float> origin; origin.clear();
    for (int d = 0; d < dim - 1; d++) origin.push_back(0);
    V.push_back(origin);
    for (int d = 0; d < dim - 1; d++) {
        origin[d] = 1.0;
        V.push_back(origin);
        origin[d] = 0.0;
    }
    innerPoint.clear();
    for (int d = 0; d < dim - 1; d++) innerPoint.push_back(1.0 / (float)dim);
    innerPoint.push_back(0.0);
    return;
}

void region::ComputeHP(vector<float> &w, vector<float> &o1, vector<float> &o2, int dim) {
    w.clear();
    float o1_d = o1[dim - 1];
    float o2_d = o2[dim - 1];
    for (int d = 0; d < dim - 1; d++) {
        w.push_back((o1[d] - o1_d) - (o2[d] - o2_d));
    }
    w.push_back(o2_d - o1_d);
    return;
}

void region::FreeMem() {
    H.clear();
    vector<halfspace>().swap(H);
    V.clear();
    vector<vector<float>>().swap(V);
    innerPoint.clear();
    vector<float>().swap(innerPoint);
}