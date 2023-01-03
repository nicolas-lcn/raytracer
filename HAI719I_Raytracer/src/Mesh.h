#ifndef MESH_H
#define MESH_H


#include <vector>
#include <string>
#include "Vec3.h"
#include "Ray.h"
#include "Triangle.h"
#include "Material.h"

#include <GL/glut.h>

#include <cfloat>


// -------------------------------------------
// Basic Mesh class
// -------------------------------------------
class Square;
struct MeshVertex {
    inline MeshVertex () {}
    inline MeshVertex (const Vec3 & _p, const Vec3 & _n) : position (_p), normal (_n) , u(0) , v(0) {}
    inline MeshVertex (const MeshVertex & vertex) : position (vertex.position), normal (vertex.normal) , u(vertex.u) , v(vertex.v) {}
    inline virtual ~MeshVertex () {}
    inline MeshVertex & operator = (const MeshVertex & vertex) {
        position = vertex.position;
        normal = vertex.normal;
        u = vertex.u;
        v = vertex.v;
        return (*this);
    }
    // membres :
    Vec3 position; // une position
    Vec3 normal; // une normale
    float u,v; // coordonnees uv
};

struct MeshTriangle {
    inline MeshTriangle () {
        v[0] = v[1] = v[2] = 0;
    }
    inline MeshTriangle (const MeshTriangle & t) {
        v[0] = t.v[0];   v[1] = t.v[1];   v[2] = t.v[2];
    }
    inline MeshTriangle (unsigned int v0, unsigned int v1, unsigned int v2) {
        v[0] = v0;   v[1] = v1;   v[2] = v2;
    }
    unsigned int & operator [] (unsigned int iv) { return v[iv]; }
    unsigned int operator [] (unsigned int iv) const { return v[iv]; }
    inline virtual ~MeshTriangle () {}
    inline MeshTriangle & operator = (const MeshTriangle & t) {
        v[0] = t.v[0];   v[1] = t.v[1];   v[2] = t.v[2];
        return (*this);
    }
    // membres :
    unsigned int v[3];
};




class Mesh {
protected:
    void build_positions_array() {
        positions_array.resize( 3 * vertices.size() );
        for( unsigned int v = 0 ; v < vertices.size() ; ++v ) {
            positions_array[3*v + 0] = vertices[v].position[0];
            positions_array[3*v + 1] = vertices[v].position[1];
            positions_array[3*v + 2] = vertices[v].position[2];
        }
    }
    void build_normals_array() {
        normalsArray.resize( 3 * vertices.size() );
        for( unsigned int v = 0 ; v < vertices.size() ; ++v ) {
            normalsArray[3*v + 0] = vertices[v].normal[0];
            normalsArray[3*v + 1] = vertices[v].normal[1];
            normalsArray[3*v + 2] = vertices[v].normal[2];
        }
    }
    void build_UVs_array() {
        uvs_array.resize( 2 * vertices.size() );
        for( unsigned int vert = 0 ; vert < vertices.size() ; ++vert ) {
            uvs_array[2*vert + 0] = vertices[vert].u;
            uvs_array[2*vert + 1] = vertices[vert].v;
        }
    }
    void build_triangles_array() {
        triangles_array.resize( 3 * triangles.size() );
        for( unsigned int t = 0 ; t < triangles.size() ; ++t ) {
            triangles_array[3*t + 0] = triangles[t].v[0];
            triangles_array[3*t + 1] = triangles[t].v[1];
            triangles_array[3*t + 2] = triangles[t].v[2];
        }
    }
public:
    std::vector<MeshVertex> vertices;
    std::vector<MeshTriangle> triangles;

    std::vector< float > positions_array;
    std::vector< float > normalsArray;
    std::vector< float > uvs_array;
    std::vector< unsigned int > triangles_array;
    std::vector< Vec3 > boardingbox;

    Material material;

    void loadOFF (const std::string & filename);
    void recomputeNormals ();
    void centerAndScaleToUnit ();
    void scaleUnit ();


    virtual
    void build_arrays() {
        recomputeNormals();
        build_positions_array();
        build_normals_array();
        build_UVs_array();
        build_triangles_array();
        boardingbox = boardingBox();
    }


    void translate( Vec3 const & translation ){
        for( unsigned int v = 0 ; v < vertices.size() ; ++v ) {
            vertices[v].position += translation;
        }
    }

    void apply_transformation_matrix( Mat3 transform ){
        for( unsigned int v = 0 ; v < vertices.size() ; ++v ) {
            vertices[v].position = transform*vertices[v].position;
        }

        recomputeNormals();
        build_positions_array();
        build_normals_array();
    }

    void scale( Vec3 const & scale ){
        Mat3 scale_matrix(scale[0], 0., 0.,
                0., scale[1], 0.,
                0., 0., scale[2]); //Matrice de transformation de mise à l'échelle
        apply_transformation_matrix( scale_matrix );
    }

    void rotate_x ( float angle ){
        float x_angle = angle * M_PI / 180.;
        Mat3 x_rotation(1., 0., 0.,
                        0., cos(x_angle), -sin(x_angle),
                        0., sin(x_angle), cos(x_angle));
        apply_transformation_matrix( x_rotation );
    }

    void rotate_y ( float angle ){
        float y_angle = angle * M_PI / 180.;
        Mat3 y_rotation(cos(y_angle), 0., sin(y_angle),
                        0., 1., 0.,
                        -sin(y_angle), 0., cos(y_angle));
        apply_transformation_matrix( y_rotation );
    }

    void rotate_z ( float angle ){
        float z_angle = angle * M_PI / 180.;
        Mat3 z_rotation(cos(z_angle), -sin(z_angle), 0.,
                        sin(z_angle), cos(z_angle), 0.,
                        0., 0., 1.);
        apply_transformation_matrix( z_rotation );
    }


    void draw() const {
        if( triangles_array.size() == 0 ) return;
        GLfloat material_color[4] = {material.diffuse_material[0],
                                     material.diffuse_material[1],
                                     material.diffuse_material[2],
                                     1.0};

        GLfloat material_specular[4] = {material.specular_material[0],
                                        material.specular_material[1],
                                        material.specular_material[2],
                                        1.0};

        GLfloat material_ambient[4] = {material.ambient_material[0],
                                       material.ambient_material[1],
                                       material.ambient_material[2],
                                       1.0};

        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material_specular);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material_color);
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, material_ambient);
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, material.shininess);

        glEnableClientState(GL_VERTEX_ARRAY) ;
        glEnableClientState (GL_NORMAL_ARRAY);
        glNormalPointer (GL_FLOAT, 3*sizeof (float), (GLvoid*)(normalsArray.data()));
        glVertexPointer (3, GL_FLOAT, 3*sizeof (float) , (GLvoid*)(positions_array.data()));
        glDrawElements(GL_TRIANGLES, triangles_array.size(), GL_UNSIGNED_INT, (GLvoid*)(triangles_array.data()));
        
        // glPointSize(20.);
        // glBegin(GL_POINTS);
        // std::vector<Vec3> boardingbox = boardingBox();
        // for(size_t i = 0; i<boardingbox.size(); i++)
        // {
        //     glVertex3f( boardingbox[i][0] , boardingbox[i][1] , boardingbox[i][2] );
        //     i+= boardingbox.size()-1;
        // }
        // glEnd();
    }

    std::vector<Vec3> boardingBox() const
    {
        Vec3 bbmin, bbmax;
        bbmin = Vec3(positions_array[0], positions_array[1], positions_array[2]);
        bbmax = Vec3(positions_array[0], positions_array[1], positions_array[2]);

        for(size_t i = 3; i<positions_array.size(); i+=3)
        {
            for (int j = 0; j < 3; ++j)
            {
                if(bbmin[j] > positions_array[i+j]) bbmin[j] = positions_array[i+j];
                if(bbmax[j] < positions_array[i+j]) bbmax[j] = positions_array[i+j];
            }
        }
        bbmin -= Vec3(0.01, 0.01, 0.01);
        bbmax += Vec3(0.01, 0.01, 0.01);

        std::vector< Vec3 > boardingbox;

        boardingbox.push_back(bbmin);
        boardingbox.push_back(bbmax);

        return boardingbox;
    }

    bool isInBoardingBox(Vec3 origin, Vec3 dir) const
    {
        float tmin, tmax, tYmin, tYmax, tZmin, tZmax;
        tmin = (boardingbox[0][0] - origin[0])/dir[0]; // tXmin;
        tmax = (boardingbox[boardingbox.size()-1][0] - origin[0])/dir[0]; // tXmax;
        tYmin = (boardingbox[0][1] - origin[1])/dir[1]; 
        tYmax = (boardingbox[boardingbox.size()-1][1] - origin[1])/dir[1]; 
        tZmin = (boardingbox[0][2] - origin[2])/dir[2]; 
        tZmax = (boardingbox[boardingbox.size()-1][2] - origin[2])/dir[2]; 

        if(tmin > tmax) std::swap(tmin, tmax); // swap value if tmin is not min;
        if(tYmin > tYmax) std::swap(tYmin, tYmax); //same
        if(tZmin > tZmax) std::swap(tZmin, tZmax); //same
        if(tmin > tYmax || tmax < tYmin) return false;
        if(tmin < tYmin) tmin = tYmin;
        if(tmax > tYmax) tmax = tYmax;
        if(tmin > tZmax || tmax < tZmin) return false;
        return true;
    }

    RayTriangleIntersection intersect( Ray const & ray ) const {
        RayTriangleIntersection closestIntersection;
        closestIntersection.t = FLT_MAX;

        //check intersection with boarding box
        if(!isInBoardingBox(ray.origin(), ray.direction())) return closestIntersection;

        for (int i = 0; i < triangles.size(); ++i)
        {
            MeshTriangle meshTriangle = triangles[i];
            Triangle triangle = Triangle(vertices[triangles[i][0]].position, 
                                         vertices[triangles[i][1]].position, 
                                         vertices[triangles[i][2]].position);
            RayTriangleIntersection intersection = triangle.getIntersection(ray);
            if(intersection.intersectionExists && intersection.t<closestIntersection.t && intersection.t> 0.001f)
            {
                closestIntersection = intersection;

            }

        }
        return closestIntersection;
    }

};




#endif
