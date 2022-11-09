#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include <string>
#include "Mesh.h"
#include "Sphere.h"
#include "Square.h"


#include <GL/glut.h>


enum ObjectType{
    ObjectType_Sphere,
    ObjectType_Square,
    ObjectType_Mesh
};

enum LightType {
    LightType_Spherical,
    LightType_Quad
};


struct Light {
    Vec3 material;
    bool isInCamSpace;
    LightType type;

    Vec3 pos;
    float radius;

    Mesh quad;

    float powerCorrection;

    Light() : powerCorrection(1.0) {}

};

struct RaySceneIntersection{
    bool intersectionExists;
    unsigned int typeOfIntersectedObject;
    unsigned int objectIndex;
    float t;
    RayTriangleIntersection rayMeshIntersection;
    RaySphereIntersection raySphereIntersection;
    RaySquareIntersection raySquareIntersection;
    RaySceneIntersection() : intersectionExists(false) , t(FLT_MAX) {}
};



class Scene {
    std::vector< Mesh > meshes;
    std::vector< Sphere > spheres;
    std::vector< Square > squares;
    std::vector< Light > lights;

public:


    Scene() {
    }

    void draw() {
        // iterer sur l'ensemble des objets, et faire leur rendu :
        for( unsigned int It = 0 ; It < meshes.size() ; ++It ) {
            Mesh const & mesh = meshes[It];
            mesh.draw();
        }
        for( unsigned int It = 0 ; It < spheres.size() ; ++It ) {
            Sphere const & sphere = spheres[It];
            sphere.draw();
        }
        for( unsigned int It = 0 ; It < squares.size() ; ++It ) {
            Square const & square = squares[It];
            square.draw();
        }
    }




    RaySceneIntersection computeIntersection(Ray const & ray, float minDist) {
        RaySceneIntersection result;
        //TODO calculer les intersections avec les objets de la scene et garder la plus proche
        float t = FLT_MAX;
        for(size_t i = 0; i<spheres.size(); i++)
        {
            RaySphereIntersection intersection = spheres[i].intersect(ray);
            if(intersection.intersectionExists && intersection.t < t && intersection.t > minDist)
            {
                t = intersection.t;
                result.intersectionExists = true;
                result.typeOfIntersectedObject = ObjectType_Sphere;
                result.objectIndex = i;
                result.t = t;
                result.raySphereIntersection = intersection;
            }
        }
        for(size_t i = 0; i<squares.size(); i++)
        {
            RaySquareIntersection intersection = squares[i].intersect(ray);
            if(intersection.intersectionExists && intersection.t < t && intersection.t > minDist)
            {
                t = intersection.t;
                result.intersectionExists = true;
                result.typeOfIntersectedObject = ObjectType_Square;
                result.objectIndex = i;
                result.t = t;
                result.raySquareIntersection = intersection;

            }
        }
        return result;
    }

    RaySceneIntersection computeShadowIntersection(Ray const & ray, float maxDist)
    {
    	RaySceneIntersection result;
    	for(size_t i = 0; i<spheres.size(); ++i)
    	{
    		RaySphereIntersection intersection = spheres[i].intersect(ray);
    		if(intersection.intersectionExists && intersection.t > 0.001f && intersection.t < maxDist)
    		{
    			result.intersectionExists = true;
    			return result;
    		}
    	}
    	for(size_t i = 0; i<squares.size(); ++i)
    	{
    		RaySquareIntersection intersection = squares[i].intersect(ray);
    		if(intersection.intersectionExists && intersection.t > 0.001f && intersection.t < maxDist)
    		{
    			result.intersectionExists = true;
    			return result;
    		}
    	}
    	
    	return result;
    }


    Vec3 phong_sphere (Ray ray, Light light, RaySphereIntersection intersection, int index)
    {

		Vec3 L = light.pos - intersection.intersection;
        L.normalize();
        Vec3 normal = intersection.normal;
        Vec3 reflexion = 2 * (Vec3::dot(normal, L)) * normal - L;
        reflexion.normalize();
        Vec3 rayDir = (-1) * ray.direction();
        rayDir.normalize();
        Vec3 specular = Vec3::compProduct(
        				light.material, 
        				spheres[index].material.specular_material) * 
        				pow(Vec3::dot(reflexion, rayDir), spheres[index].material.shininess);
        Vec3 diffuse=Vec3::compProduct(
        				light.material,
        				spheres[index].material.diffuse_material)*
        			 Vec3::dot(L,normal);
        
        return specular + diffuse;
	
    }

    Vec3 phong_square(Ray ray, Light light, RaySquareIntersection intersection, int index)
    {
    	Vec3 L = light.pos - intersection.intersection;
        L.normalize();
        Vec3 normal = intersection.normal;
        Vec3 reflexion = 2 * (Vec3::dot(normal, L)) * normal - L;
        reflexion.normalize();
        Vec3 rayDir = (-1) * ray.direction();
        rayDir.normalize();
        Vec3 specular = Vec3::compProduct(light.material, 
        								  squares[index].material.specular_material) * 
        								  pow(Vec3::dot(reflexion, rayDir), 
        								  	  squares[index].material.shininess);
        Vec3 diffuse =  Vec3::compProduct(light.material,
        								  squares[index].material.diffuse_material)  *
        								  Vec3::dot(L,normal);

        return diffuse + specular;
    }

    Vec3 rayTraceRecursive( Ray ray , int NRemainingBounces) {

        RaySceneIntersection raySceneIntersection = computeIntersection(ray, 4.9f);
        Vec3 color;
        int light_number = 0;
        if(raySceneIntersection.intersectionExists)
        {
        	
            if(raySceneIntersection.typeOfIntersectedObject == ObjectType_Sphere)
            {
            	for (int i = 0; i < lights.size(); ++i)
            	{
            		Vec3 L = lights[i].pos - raySceneIntersection.raySphereIntersection.intersection;
            		float shadowMinDist = L.length();
            		L.normalize();
            		Ray shadowRay = Ray(raySceneIntersection.raySphereIntersection.intersection, L);
            		RaySceneIntersection shadowIntersection = computeShadowIntersection(shadowRay, shadowMinDist);
            		if(shadowIntersection.intersectionExists)
        			{
        				return Vec3(0,0,0);
        			}
                    else 
                    {
                    	color += phong_sphere(ray, lights[i], raySceneIntersection.raySphereIntersection, raySceneIntersection.objectIndex);
                    	
                    }
            	}
            }
            else if(raySceneIntersection.typeOfIntersectedObject == ObjectType_Square)
            {

            	for (int i = 0; i < lights.size(); ++i)
            	{
            		{
	            		Vec3 L = lights[i].pos - raySceneIntersection.raySquareIntersection.intersection;
	            		float shadowMinDist = L.length();
	            		L.normalize();
	            		Ray shadowRay = Ray(raySceneIntersection.raySquareIntersection.intersection, L);
	            		RaySceneIntersection shadowIntersection = computeShadowIntersection(shadowRay, shadowMinDist);
	            		
	                    if(shadowIntersection.intersectionExists)
	                	{
	                		return Vec3(0,0,0);
	                	}
	                    else 
	                    {
	                    	color += phong_square(ray, lights[i], raySceneIntersection.raySquareIntersection, raySceneIntersection.objectIndex);
	                    	
	                    }
                	}
            	}
            }
        }
        return color;
    }


    Vec3 rayTrace( Ray const & rayStart ) {
        
        Vec3 color;
        color = rayTraceRecursive(rayStart, 1);
        return color;
    }

    void setup_single_sphere() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(-5,5,5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }
        {
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(0. , 0. , 0.);
            s.m_radius = 1.f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,1.,1 );
            s.material.specular_material = Vec3( 0.2,0.2,0.2 );
            s.material.shininess = 20;
        }
    }

    void setup_single_square() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(-5,5,5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }

        {
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.8,0.8,0.8 );
            s.material.specular_material = Vec3( 0.8,0.8,0.8 );
            s.material.shininess = 20;
        }
    }

    void setup_cornell_box(){
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3( 0.0, 1.97, 0.0 );
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }

        { //Back Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.5,0.,0.7 );
            s.material.specular_material = Vec3( 1.,1.,1. );
            s.material.shininess = 16;
        }

        { //Left Wall

            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.rotate_y(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.shininess = 16;
        }

        { //Right Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.0,1.0,0.0 );
            s.material.specular_material = Vec3( 0.0,1.0,0.0 );
            s.material.shininess = 16;
        }

        { //Floor
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.0,1.0,0.0 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }

        { //Ceiling
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.7,0.7,0.7 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }

        { //Front Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(180);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.0,1.0,1.0 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }


        { //GLASS Sphere

            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1.0, -1.25, 0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.shininess = 16;
            s.material.transparency = 1.0;
            s.material.index_medium = 1.4;
        }


        { //MIRRORED Sphere
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1.0, -1.25, -0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Glass;
            s.material.diffuse_material = Vec3( 0.,1.,0. );
            s.material.specular_material = Vec3(  0.,1.,0. );
            s.material.shininess = 16;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;
        }

    }

    void setup_two_spheres()
    {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(-5,5,5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }
        {
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1. , 0. , 0.);
            s.m_radius = 1.f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 0.2,0.2,0.2 );
            s.material.shininess = 20;
        }
        {
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1. , 0. , 0.);
            s.m_radius = 1.f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 0.,1.,0.);
            s.material.specular_material = Vec3( 0.2,0.2,0.2 );
            s.material.shininess = 20;
        }
    }

};



#endif
