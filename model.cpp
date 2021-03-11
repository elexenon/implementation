#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "model.h"

Model::Model(const char *filename) : verts_(), faces_() {
    std::ifstream in;
    in.open (filename, std::ifstream::in);
    if (in.fail()) return;
    std::string line;
    while (!in.eof()) {
        std::getline(in, line);
        std::istringstream iss(line.c_str());
        char trash;
        int  itrash;
        std::string strash;
        if (!line.compare(0, 2, "v ")) {
            iss >> trash;
            Vec3f v;
            for (int i=0;i<3;i++) iss >> v[i];
            verts_.push_back(v);
        } else if(!line.compare(0, 2, "vt")){
            iss >> strash;
            Vec2f uv;
            for (int i=0;i<2;i++) iss >> uv[i];
            uvs_.push_back(uv);
        } else if (!line.compare(0, 2, "f ")) {
            std::vector<std::pair<int,int>> f;
            iss >> trash;
            int idx, vt;
            while (iss >> idx >> trash >> vt >> trash >> itrash) {
                idx--; // 
                vt--;
                f.emplace_back(idx, vt);
            }
            faces_.push_back(f);
        }
    }
    std::cerr << "# v# " << verts_.size() << "# vt# " << uvs_.size() << " f# "  << faces_.size() << std::endl;
}

Model::~Model() {
}

int Model::nverts() {
    return (int)verts_.size();
}

int Model::nfaces() {
    return (int)faces_.size();
}

std::vector<std::pair<int,int>> Model::face(int idx) {
    return faces_[idx];
}

Vec3f Model::vert(int i) {
    return verts_[i];
}

Vec2f Model::uv(int i)
{
    return uvs_[i];
}
