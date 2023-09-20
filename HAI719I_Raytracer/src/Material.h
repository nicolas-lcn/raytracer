#ifndef MATERIAL_H
#define MATERIAL_H

#include "imageLoader.h"
#include "Vec3.h"
#include <cmath>
#include <GL/glew.h>
#include <GL/glut.h>
 //todo use imageloader

enum MaterialType {
    Material_Diffuse_Blinn_Phong ,
    Material_Glass,
    Material_Mirror,
    Material_Texture,
    Material_Transparent
};


struct Material {
    Vec3 ambient_material;
    Vec3 diffuse_material;
    Vec3 specular_material;
    double shininess;
    ppmLoader::ImageRGB texture;

    float index_medium;
    float transparency;

    MaterialType type;

    Material() {
        type = Material_Diffuse_Blinn_Phong;
        transparency = 0.0;
        index_medium = 1.0;
        ambient_material = Vec3(0., 0., 0.);
    }

    void loadTexture2DFromFilePath(const std::string& path) {
       ppmLoader::ImageRGB image;
       ppmLoader::load_ppm(image, path);
       texture = image;
    }

    Vec3 sampleTexture(const float u, const float v) {
        Vec3 color;
        int i = floor(u*this->texture.w);
        int j = floor((1-v)*this->texture.h);
        int index = (i+j*this->texture.w);
        ppmLoader::RGB rgb = this->texture.data[index];
        color = Vec3((float)rgb.r/255.0, (float)rgb.g/255.0, (float)rgb.b/255.0);
        return color;
    }

};



#endif // MATERIAL_H
