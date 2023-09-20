#ifndef SQUARE_H
#define SQUARE_H
#include "Vec3.h"
#include <vector>
#include "Mesh.h"
#include <cmath>

struct RaySquareIntersection{
    bool intersectionExists;
    float t;
    float u,v;
    Vec3 intersection;
    Vec3 normal;
};


class Square : public Mesh {
public:
    Vec3 m_normal;
    Vec3 m_bottom_left;
    Vec3 m_right_vector;
    Vec3 m_up_vector;

    Square() : Mesh() {}
    Square(Vec3 const & bottomLeft , Vec3 const & rightVector , Vec3 const & upVector , float width=1. , float height=1. ,
           float uMin = 0.f , float uMax = 1.f , float vMin = 0.f , float vMax = 1.f) : Mesh() {
        setQuad(bottomLeft, rightVector, upVector, width, height, uMin, uMax, vMin, vMax);
    }

    Vec3 bottomLeft() const { return this->vertices[0].position; }
    Vec3 bottomRight() const { return this->vertices[1].position; }
    Vec3 upRight() const { return this->vertices[2].position; }
    Vec3 upLeft() const { return this->vertices[3].position; }
    Vec3 normal() const { return Vec3::cross((bottomRight() - bottomLeft()) , (upLeft() - bottomLeft())); }

    void setQuad( Vec3 const & bottomLeft , Vec3 const & rightVector , Vec3 const & upVector , float width=1. , float height=1. ,
                  float uMin = 0.f , float uMax = 1.f , float vMin = 0.f , float vMax = 1.f) {
        m_right_vector = rightVector;
        m_up_vector = upVector;
        m_normal = Vec3::cross(rightVector , upVector);
        m_bottom_left = bottomLeft;

        m_normal.normalize();
        m_right_vector.normalize();
        m_up_vector.normalize();

        m_right_vector = m_right_vector*width;
        m_up_vector = m_up_vector*height;

        vertices.clear();
        vertices.resize(4);
        vertices[0].position = bottomLeft;                                      vertices[0].u = uMin; vertices[0].v = vMin;
        vertices[1].position = bottomLeft + m_right_vector;                     vertices[1].u = uMax; vertices[1].v = vMin;
        vertices[2].position = bottomLeft + m_right_vector + m_up_vector;       vertices[2].u = uMax; vertices[2].v = vMax;
        vertices[3].position = bottomLeft + m_up_vector;                        vertices[3].u = uMin; vertices[3].v = vMax;
        vertices[0].normal = vertices[1].normal = vertices[2].normal = vertices[3].normal = m_normal;
        triangles.clear();
        triangles.resize(2);
        triangles[0][0] = 0;
        triangles[0][1] = 1;
        triangles[0][2] = 2;
        triangles[1][0] = 0;
        triangles[1][1] = 2;
        triangles[1][2] = 3;


    }
    bool isInSquare(Vec3 pos) const
    {
        Vec3 bottom = Vec3::cross(bottomRight() - bottomLeft(), pos - bottomLeft());
        Vec3 right = Vec3::cross(upRight() - bottomRight(), pos - bottomRight());
        Vec3 up = Vec3::cross(upLeft() - upRight(), pos - upRight());
        Vec3 left = Vec3::cross(bottomLeft() - upLeft(), pos - upLeft());

        bool direction = (Vec3::dot(bottom, normal()) > 0);

        bool sameDirection = (Vec3::dot(right, normal()) > 0) == direction && 
                             (Vec3::dot(left, normal())  > 0) == direction &&
                             (Vec3::dot(up, normal())    > 0) == direction;
        return sameDirection;
    }

    RaySquareIntersection intersect(const Ray &ray) const {
        RaySquareIntersection intersection;
        Vec3 p = bottomLeft() - ray.origin();
        float t = ( Vec3::dot(p, normal()))/Vec3::dot(ray.direction(), normal());
        Vec3 point = ray.direction() * t + ray.origin();
        float width = (bottomRight() - bottomLeft()).length();
        float height = (upLeft() - bottomLeft()).length();
        if(t>0 && isInSquare(point))
        {
            intersection.intersectionExists = true;
            intersection.t = t;
            intersection.intersection = point;
            intersection.normal = normal();
            intersection.normal.normalize();
            float u,v;
            Vec3 X = bottomRight() - bottomLeft();
            Vec3 Y = upLeft() - bottomLeft();
            Vec3 p_bottom = point - bottomLeft();
            Vec3 projP_X = (Vec3::dot(p_bottom, X)/X.squareLength())*X;
            Vec3 projP_Y = (Vec3::dot(p_bottom, Y)/Y.squareLength())*Y;
            u = projP_X.length()/X.length();
            v = projP_Y.length()/Y.length();
            intersection.u = u;
            intersection.v = v;
            return intersection;
        }
        else
        {
            intersection.intersectionExists = false;
            return intersection;
        }
         
    }
};
#endif // SQUARE_H
