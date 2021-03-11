#include "tgaimage.h"
#include "geometry.h"
#include "model.h"

#include <cmath>
#include <iostream>

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0,   255, 0,   255);
const TGAColor blue  = TGAColor(0,   0,   255, 255);
const int width  = 800;
const int height = 800;

Vec3f world2screen(Vec3f v) {
    // 向上舍入
    return Vec3f(int((v.x+1.)*width/2.+.5), int((v.y+1.)*height/2.+.5), v.z);
}

std::tuple<int,int> computeUV(float u, float v, TGAImage* textureTga) {
    return {u*textureTga->get_width(), v*textureTga->get_height()};
}

Vec3f barycentric(Vec3f *pts, Vec3f P) {
    Vec3f u = cross(Vec3f(pts[2][0]-pts[0][0], pts[1][0]-pts[0][0], pts[0][0]-P[0]), Vec3f(pts[2][1]-pts[0][1], pts[1][1]-pts[0][1], pts[0][1]-P[1]));
    if (std::abs(u[2])<1) return Vec3f(-1,1,1);
    return {1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z};
}

// Midpoint Bresenham
void line(Vec2i t0, Vec2i t1, TGAImage &image, TGAColor color) { 
    int x0 = t0.x, y0 = t0.y;
    int x1 = t1.x, y1 = t1.y;
    // 确保x为最大位移方向
    bool steep = false; 
    if (std::abs(x0-x1)<std::abs(y0-y1)) { 
        std::swap(x0, y0); 
        std::swap(x1, y1); 
        steep = true; 
    } 
    // 确保计算方向为正向
    if(x0 > x1) {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }

    int dX = std::round(fabs(x1 - x0));
    int dY = std::round(fabs(y1 - y0));
    int delta = dX - 2 * dY;

    int dStepUp = 2 * (dX - dY);
    int dStepDown = -2 * dY;

    int x = x0, y = y0;

    for(int i = x; i <= x1; i++) {
        !steep ? image.set(i, y, color)
               : image.set(y, i, color);
        if(delta < 0) {
            y += (y0<y1?1:-1);
            delta += dStepUp;
        } else {
            delta += dStepDown;
        }
    }
} 
// Sweeping Line
void triangle0(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage& image, TGAColor color) {
    if(t0.y > t1.y) std::swap(t0, t1);
    if(t0.y > t2.y) std::swap(t0, t2);
    if(t1.y > t2.y) std::swap(t1, t2);

    int h1 = t2.y - t0.y;
    int h2 = t1.y - t0.y;
    int h3 = t2.y - t1.y;

    for(int i = 0; i <= h1; i++) {
        bool halfPart = i > h2 || t1.y == t0.y;

        int segmentHeight = halfPart ? h3 : h2;

        float alpha = (float)i / h1;
        float beta  = (float)(i - (halfPart ? h2 : 0)) / segmentHeight;

        Vec2i startPoint = t0 + (t2-t0)*alpha;
        Vec2i endPoint   = halfPart ? t1 + (t2-t1)*beta : t0 + (t1-t0)*beta;

        if(startPoint.x > endPoint.x) std::swap(startPoint, endPoint);
        for(int j = startPoint.x; j <= endPoint.x; j++) {
            image.set(j, t0.y+i, color);
        }
    }
}

// 通过顶点颜色插值
void triangle0(Vec3f *pts, float *zbuffer, TGAImage &image, TGAColor* colors) {
    Vec2f bboxmin( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max());
    Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    Vec2f clamp(image.get_width()-1, image.get_height()-1);
    for (int i=0; i<3; i++) {
        for (int j=0; j<2; j++) {
            bboxmin[j] = std::max(0.f,      std::min(bboxmin[j], pts[i][j]));
            bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
        }
    }
    Vec3f P;
    for (P.x=bboxmin.x; P.x<=bboxmax.x; P.x++) {
        for (P.y=bboxmin.y; P.y<=bboxmax.y; P.y++) {
            Vec3f bc_screen = barycentric(pts, P);
            if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0) continue;
            P.z = 0;
            TGAColor color;
            for (int i=0; i<3; i++) {
                P.z += pts[i][2]*bc_screen[i];
                for(int j = 0;j < 4;j++){
                    color.bgra[j] += colors[i].bgra[j]*bc_screen[i];
                }
            }
            if (zbuffer[int(P.x+P.y*width)]<P.z) {
                zbuffer[int(P.x+P.y*width)] = P.z;
                image.set(P.x, P.y, color);
            }
        }
    }
}

// 通过顶点uv插值
void triangle(std::vector<std::pair<int,int>>& face, Model* model, TGAImage* textureTga,
              float* zbuffer, TGAImage& image, float intensity) {
    Vec2f bboxmin( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max());
    Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    Vec2f clamp(image.get_width()-1, image.get_height()-1);

    Vec3f pts[3];
    for(int i=0; i<3; i++) {
        pts[i] = world2screen(model->vert(face[i].first));
        std::cout << pts[i] << std::endl;
    }

    for (int i=0; i<3; i++) {
        for (int j=0; j<2; j++) {
            bboxmin[j] = std::max(0.f,      std::min(bboxmin[j], pts[i][j]));
            bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
        }
    }

    Vec3f P;
    for (P.x=bboxmin.x; P.x<=bboxmax.x; P.x++) {
        for (P.y=bboxmin.y; P.y<=bboxmax.y; P.y++) {
            Vec3f bc_screen = barycentric(pts, P);
            if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0) continue;

            P.z = 0;
            float u = 0.0,v = 0.0;
            TGAColor color;

            for (int i=0; i<3; i++) {
                P.z += pts[i][2]*bc_screen[i];
                u += model->uv(face[i].second)[0]*bc_screen[i];
                v += model->uv(face[i].second)[1]*bc_screen[i];
            }

            auto[X,Y] = computeUV(u,v,textureTga);

            color = textureTga->get(X, Y);
            color = color * intensity;

            if (zbuffer[int(P.x+P.y*width)]<P.z) {
                zbuffer[int(P.x+P.y*width)] = P.z;
                image.set(P.x, P.y, color);
            }
        }
    }
}

int main(int argc, char** argv) {
    Model* model = nullptr;
    TGAImage* textureTga = nullptr;
    if (2==argc) {
        model = new Model(argv[1]);
    } else if(3==argc){
        model = new Model(argv[1]);
        textureTga = new TGAImage;
        textureTga->read_tga_file(argv[2]);
        textureTga->flip_vertically();
    } else {
        model = new Model("obj/african_head.obj");
    }

    float *zbuffer = new float[width*height];
    for (int i=width*height; i--; zbuffer[i] = -std::numeric_limits<float>::max());

    TGAImage image(width, height, TGAImage::RGB);

    Vec3f light_dir(0, 0, -1);

// //按颜色插值
//    for (int i=0; i<model->nfaces(); i++) {
//        std::vector<std::pair<int,int>> face = model->face(i);
//        Vec3f pts[3];
//        TGAColor colors[3];

//        for (int i=0; i<3; i++) {
//            pts[i] = model->vert(face[i].first);
//            int X = model->uv(face[i].second)[0] * textureTga->get_width();
//            int Y = model->uv(face[i].second)[1] * textureTga->get_height();
//            colors[i] = textureTga->get(X,Y);
//        }

//        Vec3f n = cross(pts[2]-pts[0], pts[1]-pts[0]);

//        for(int i=0; i<3; i++) pts[i] = world2screen(pts[i]);

//        n.normalize();
//        float intensity = n*light_dir;
//        if (intensity>0) {
//            for(int i=0;i<3;i++) colors[i] = colors[i]*intensity;
//            triangle0(pts, zbuffer, image, colors);
//        }
//    }

// 按uv插值
    for (int i=0; i<model->nfaces(); i++) {
        std::vector<std::pair<int,int>> face = model->face(i);
        Vec3f pts[3];

        for (int i=0; i<3; i++) pts[i] = model->vert(face[i].first);

        Vec3f n = cross(pts[2]-pts[0], pts[1]-pts[0]);

        n.normalize();
        float intensity = n*light_dir;
        if (intensity>0) {
            triangle(face, model, textureTga, zbuffer, image, intensity);
        }
    }

    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");
    delete model;
    return 0;
}
