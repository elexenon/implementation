#ifndef __MODEL_H__
#define __MODEL_H__

#include <vector>
#include "geometry.h"

class Model {
private:
    // 从上到下存储所有顶点坐标
    std::vector<Vec3f> verts_;
    // 从上到下存储所有uv坐标
    std::vector<Vec2f> uvs_;
    // 存放所有的面，每个面包含一些顶点，每个顶点为一个std::pair
    std::vector<std::vector<std::pair<int,int>> > faces_;
public:
    Model(const char *filename);
    ~Model();
    int nverts();
    int nuvs();
    int nfaces();
    Vec3f vert(int i);
    Vec2f uv(int i);
    std::vector<std::pair<int,int>> face(int idx);
};

#endif //__MODEL_H__
