#ifndef TRIANGLE_H
#define TRIANGLE_H
#include "Vec3.h"
#include "Ray.h"
#include "Plane.h"

struct RayTriangleIntersection{
    bool intersectionExists;
    float t;
    float w0,w1,w2;
    unsigned int tIndex;
    Vec3 intersection;
    Vec3 normal;
};

class Triangle {
private:
    Vec3 m_c[3] , m_normal;
    float area;
public:
    Triangle() {}
    Triangle( Vec3 const & c0 , Vec3 const & c1 , Vec3 const & c2 ) {
        m_c[0] = c0;
        m_c[1] = c1;
        m_c[2] = c2;
        updateAreaAndNormal();
    }
    void updateAreaAndNormal() {
        Vec3 nNotNormalized = Vec3::cross( m_c[1] - m_c[0] , m_c[2] - m_c[0] );
        float norm = nNotNormalized.length();
        m_normal = nNotNormalized / norm;
        area = norm / 2.f;
    }
    void setC0( Vec3 const & c0 ) { m_c[0] = c0; } // remember to update the area and normal afterwards!
    void setC1( Vec3 const & c1 ) { m_c[1] = c1; } // remember to update the area and normal afterwards!
    void setC2( Vec3 const & c2 ) { m_c[2] = c2; } // remember to update the area and normal afterwards!

    Vec3 const & normal() const { return m_normal; }
    Vec3 const & c0() const { return m_c[0]; }
    Vec3 const & c1() const { return m_c[1]; }
    Vec3 const & c2() const { return m_c[2]; }

    Vec3 projectOnSupportPlane( Vec3 const & p ) const {
        Vec3 result;
        //TODO completer
        return result;
    }
    float squareDistanceToSupportPlane( Vec3 const & p ) const {
        float result;
        //TODO completer
        return result;
    }
    float distanceToSupportPlane( Vec3 const & p ) const { return sqrt( squareDistanceToSupportPlane(p) ); }
    bool isParallelTo( Line const & L ) const {
        bool result;
        float dirDotNorm = Vec3::dot(L.direction(), normal());
        result = (dirDotNorm == 0);
        return result;
    }
    Vec3 getIntersectionPointWithSupportPlane( Line const & L, float &_t) const {
        Vec3 result;
        if(! isParallelTo(L))
        {
            _t = ( Vec3::dot((c0() - L.origin()), normal()))/Vec3::dot(L.direction(), normal());
            if(_t>0) result = L.direction() * _t + L.origin();
        }
        return result;
    }
    void computeBarycentricCoordinates( Vec3 const & p , 
                                        float & u0 , 
                                        float & u1 , 
                                        float & u2 ,
                                        bool & isInside) const {
        Vec3 c0c1 = (c1() - c0());
        Vec3 c1c2 = (c2() - c1());
        Vec3 c2c0 = (c0() - c2());


        Vec3 c0pc1normal = Vec3::cross(c0c1, (p-c0()));
        Vec3 c1pc2normal = Vec3::cross(c1c2, (p-c1()));
        Vec3 c2pc0normal = Vec3::cross(c2c0, (p-c2()));

        isInside = (Vec3::dot(c1pc2normal, normal()) > 0) && (Vec3::dot(c0pc1normal, normal()) > 0) && (Vec3::dot(c2pc0normal, normal()) > 0);

        // u0 = Vec3::dot(normal(), p0);
        // u2 = Vec3::dot(normal(), p2);
        // u1 = 1 - u0 - u2;

        float area0 = c1pc2normal.length()/2.f;
        float area1 = c2pc0normal.length()/2.f;
        float area2 = c0pc1normal.length()/2.f;

        u0 = area0/this->area;
        u1 = area1/this->area;
        u2 = area2/this->area;
    }

    RayTriangleIntersection getIntersection( Ray const & ray ) const {
        RayTriangleIntersection result;
        // 1) check that the ray is not parallel to the triangle:

        // 2) check that the triangle is "in front of" the ray:

        // 3) check that the intersection point is inside the triangle:
        // CONVENTION: compute u,v such that p = w0*c0 + w1*c1 + w2*c2, check that 0 <= w0,w1,w2 <= 1

        // 4) Finally, if all conditions were met, then there is an intersection! :

        Line rayLine = (Line) ray;
        if(isParallelTo(rayLine))
        {
            result.intersectionExists = false;
            return result;
        }
        float t;
        Vec3 planeIntersection = getIntersectionPointWithSupportPlane(rayLine, t);
        if(t<0)
        {
            result.intersectionExists = false;
            return result;
        }
        float w0,w1,w2;
        bool isInside;
        computeBarycentricCoordinates(planeIntersection, w0,w1,w2, isInside);
        if(( 0.f<=w0 && w0<=1.f) &&
           ( 0.f<=w1 && w1<=1.f)  &&
           ( 0.f<=w2 && w2<=1.f) && isInside)
        {
            result.intersectionExists = true;
            result.t = t;
            result.w0 = w0;
            result.w1 = w1;
            result.w2 = w2;
            result.intersection = Vec3(planeIntersection[0] * w0, planeIntersection[1] * w1, planeIntersection[2] * w2);
            Vec3 norm = Vec3(normal()[0], normal()[1], normal()[2]);
            norm.normalize();
            result.normal = norm;
        }
        else
        {
            result.intersectionExists = false;
            return result;
        }


        return result;
       
    }
};
#endif
