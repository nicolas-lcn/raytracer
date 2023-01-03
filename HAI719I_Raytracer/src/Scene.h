#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include <string>
#include "Mesh.h"
#include "Sphere.h"
#include "Square.h"


#include <GL/glut.h>

#define LIGHT_SAMPLES 10.0f


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
    Vec3 direction;
    float radius;

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

    bool displaySoftShadow = true;

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




    RaySceneIntersection computeIntersection(Ray const & ray, float minDist) 
    {
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
        for(size_t i = 0; i<meshes.size(); i++)
        {
            RayTriangleIntersection intersection = meshes[i].intersect(ray);
            if(intersection.intersectionExists && intersection.t < t)
            {
                t = intersection.t;
                result.intersectionExists = true;
                result.typeOfIntersectedObject = ObjectType_Mesh;
                result.objectIndex = i;
                result.t = t;
                result.rayMeshIntersection = intersection;
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
                result.typeOfIntersectedObject = ObjectType_Sphere;
                result.objectIndex = i;
                result.raySphereIntersection = intersection;
    			return result;
    		}
    	}
    	for(size_t i = 0; i<squares.size(); ++i)
    	{
    		RaySquareIntersection intersection = squares[i].intersect(ray);
    		if(intersection.intersectionExists && intersection.t > 0.001f && intersection.t < maxDist)
    		{
    			result.intersectionExists = true;
    			result.typeOfIntersectedObject = ObjectType_Square;
                result.objectIndex = i;
                result.raySquareIntersection = intersection;
    			return result;
    		}
    	}

        for(size_t i = 0; i<meshes.size(); ++i)
        {
            RayTriangleIntersection intersection = meshes[i].intersect(ray);
            if(intersection.intersectionExists && intersection.t > 0.001f && intersection.t < maxDist)
            {
                result.intersectionExists = true;
                result.typeOfIntersectedObject = ObjectType_Mesh;
                result.objectIndex = i;
                result.rayMeshIntersection = intersection;
                return result;
            }
        }
    	
    	return result;
    }

    Vec3 getRefractedRayDir(Vec3 point, Vec3 normal, Vec3 rayDir, double index_medium)
    {
        double cosIncident = Vec3::dot(rayDir, normal);
        double n1, n2;
        if(index_medium == 1.0) return rayDir;
        if(cosIncident > 0.0)
        {
            n1 = index_medium;
            n2 = 1.0;
            cosIncident = -cosIncident;
            normal = -1 * normal;
        }
        else{
            n1 = 1.0;
            n2 = index_medium;
        }
        double n = n1/n2;
        double sinTheta = n *(1 - cosIncident * cosIncident);
        Vec3 refractDir;
        if(sinTheta<1.0)
        {
            double cosTheta = sqrt(1.0 - sinTheta * sinTheta);
            refractDir = n * rayDir + (n * cosIncident - cosTheta) * normal;
            refractDir.normalize();
        }
        return refractDir;
    }


    Vec3 computePhongIllumination(
        Ray ray, 
        Light light, 
        Vec3 normal,
        Vec3 intersection,
        Vec3 ambient_material,
        Vec3 diffuse_material,
        Vec3 specular_material,
        double shininess)
    {
        Vec3 L = light.pos - intersection;
        L.normalize();
        Vec3 reflexion = 2 * (Vec3::dot(L, normal)) * normal - L;
        reflexion.normalize();
        Vec3 rayDir = ray.origin() - intersection;
        rayDir.normalize();
        float alpha = std::max(Vec3::dot(reflexion, rayDir), 0.f);
        float theta = std::max(Vec3::dot(L,normal), 0.f);
        Vec3 ambient = Vec3::compProduct(light.material, ambient_material);
        Vec3 specular = Vec3::compProduct(
                        light.material, 
                        specular_material) * 
                        pow(alpha, shininess);
        Vec3 diffuse=Vec3::compProduct(
                        light.material,
                        diffuse_material)*theta;
        
        return ambient + specular + diffuse;
    }

    Vec3 computeHardShadow(Ray ray, Light light, Vec3 normal,
        Vec3 intersection,
        Vec3 ambient_material,
        Vec3 diffuse_material,
        Vec3 specular_material,
        double shininess)
    {
        Vec3 L = light.pos - intersection;
        float shadowMinDist = L.length();
        L.normalize();
        Ray shadowRay = Ray(intersection, L);
        RaySceneIntersection shadowIntersection = computeShadowIntersection(shadowRay, shadowMinDist);
        if(shadowIntersection.intersectionExists) 
            return Vec3(0,0,0);
        else 
            return computePhongIllumination(ray, light, normal, intersection, ambient_material, diffuse_material, specular_material, shininess);
                
    }

    float computeSoftShadow(Light light, Vec3 intersection)
    {

        float counter = 0;
        for(float k = 0; k<LIGHT_SAMPLES; k++)
        {
        	float x =(float)(rand()/(float)(RAND_MAX / (light.radius)));
            float y =(float)(rand()/(float)(RAND_MAX / (light.radius)));
        	float z =(float)(rand()/(float)(RAND_MAX / (light.radius)));
            Vec3 pos = Vec3(light.pos[0] + x, light.pos[1] + y, light.pos[2] + z);

			Vec3 L = pos - intersection;
			L.normalize();
			Ray shadowRay = Ray(intersection, L);
			RaySceneIntersection shadowIntersection = computeShadowIntersection(shadowRay, 0.9f);
			if(shadowIntersection.intersectionExists) 
            {
                if(shadowIntersection.typeOfIntersectedObject == ObjectType_Sphere)
                {
                    if(spheres[shadowIntersection.objectIndex].material.type != Material_Glass)
                        counter++;
                }
                else if(shadowIntersection.typeOfIntersectedObject == ObjectType_Square)
                {
                    if(squares[shadowIntersection.objectIndex].material.type != Material_Glass)
                        counter++;
                }
                else //Mesh
                {
                    if(meshes[shadowIntersection.objectIndex].material.type != Material_Glass)
                        counter ++;
                }
            }
			
		}
        float shadow = counter/LIGHT_SAMPLES;
		return (1.0f - shadow);
    }


    Vec3 rayTraceRecursive( Ray ray , int NRemainingBounces, float minDist) 
    {
        if(NRemainingBounces == 0) return Vec3();
        RaySceneIntersection raySceneIntersection = computeIntersection(ray, minDist);
        Vec3 color;
        if(!raySceneIntersection.intersectionExists) return color;
        if(raySceneIntersection.typeOfIntersectedObject == ObjectType_Sphere)
        {
            Vec3 normal = raySceneIntersection.raySphereIntersection.normal;
            Vec3 point = raySceneIntersection.raySphereIntersection.intersection;
            float intersected = 0;
            for (int i = 0; i < lights.size(); ++i)
            {
                if(spheres[raySceneIntersection.objectIndex].material.type == Material_Mirror)
                {
                    Vec3 reflectionDir = ray.direction() - (2 * normal * Vec3::dot(ray.direction(), normal));
                    reflectionDir.normalize();
                    Ray reflected = Ray(point, reflectionDir);
                    color += rayTraceRecursive(reflected, NRemainingBounces - 1, 0.0001f);
                }
                else if (spheres[raySceneIntersection.objectIndex].material.type == Material_Glass)
                {
                    Vec3 refractDir = getRefractedRayDir(point, normal, ray.direction(), spheres[raySceneIntersection.objectIndex].material.index_medium);
                    if(refractDir.length()>0)
                    {
                        Ray refracted = Ray(point, refractDir);
                        color += rayTraceRecursive(refracted, NRemainingBounces - 1, 0.0001f);
                    }
                }
                else
                {
                    if(displaySoftShadow)
                    {
                        color = computePhongIllumination(ray, 
                                    lights[i], 
                                    raySceneIntersection.raySphereIntersection.normal, 
                                    raySceneIntersection.raySphereIntersection.intersection,
                                    spheres[raySceneIntersection.objectIndex].material.ambient_material,
                                    spheres[raySceneIntersection.objectIndex].material.diffuse_material,
                                    spheres[raySceneIntersection.objectIndex].material.specular_material,
                                    spheres[raySceneIntersection.objectIndex].material.shininess
                                );
                        float shadow = computeSoftShadow(lights[i], raySceneIntersection.raySphereIntersection.intersection);
                        color *= shadow;
                    }
                    else
                    {
                        color += computeHardShadow(ray, 
                                  lights[i], 
                                  raySceneIntersection.raySphereIntersection.normal, 
                                  raySceneIntersection.raySphereIntersection.intersection,
                                  spheres[raySceneIntersection.objectIndex].material.ambient_material,
                                  spheres[raySceneIntersection.objectIndex].material.diffuse_material,
                                  spheres[raySceneIntersection.objectIndex].material.specular_material,
                                  spheres[raySceneIntersection.objectIndex].material.shininess
                                 );
                    }
                
                }
            }

        }
        else if(raySceneIntersection.typeOfIntersectedObject == ObjectType_Square)
        {
            float intersected = 0;
            
            for (int i = 0; i < lights.size(); ++i)
            {
                if(displaySoftShadow)
                {
                    color = computePhongIllumination(ray, 
                            lights[i], 
                            raySceneIntersection.raySquareIntersection.normal, 
                            raySceneIntersection.raySquareIntersection.intersection,
                            squares[raySceneIntersection.objectIndex].material.ambient_material,
                            squares[raySceneIntersection.objectIndex].material.diffuse_material,
                            squares[raySceneIntersection.objectIndex].material.specular_material,
                            squares[raySceneIntersection.objectIndex].material.shininess);
                    float shadow = computeSoftShadow(lights[i], raySceneIntersection.raySquareIntersection.intersection);
                    color *= shadow;
                }
                else
                {
                    color += computeHardShadow(ray, 
                                      lights[i], 
                                      raySceneIntersection.raySquareIntersection.normal, 
                                      raySceneIntersection.raySquareIntersection.intersection,
                                      squares[raySceneIntersection.objectIndex].material.ambient_material,
                                      squares[raySceneIntersection.objectIndex].material.diffuse_material,
                                      squares[raySceneIntersection.objectIndex].material.specular_material,
                                      squares[raySceneIntersection.objectIndex].material.shininess
                                     );
                }
                
            }
            
        }
        else if(raySceneIntersection.typeOfIntersectedObject == ObjectType_Mesh)
        {
            float intersected = 0;
            
            for (int i = 0; i < lights.size(); ++i)
            {
                color = meshes[raySceneIntersection.objectIndex].material.diffuse_material;
                if(displaySoftShadow)
                {
                    color = computePhongIllumination(ray, 
                            lights[i], 
                            raySceneIntersection.rayMeshIntersection.normal, 
                            raySceneIntersection.rayMeshIntersection.intersection,
                            meshes[raySceneIntersection.objectIndex].material.ambient_material,
                            meshes[raySceneIntersection.objectIndex].material.diffuse_material,
                            meshes[raySceneIntersection.objectIndex].material.specular_material,
                            meshes[raySceneIntersection.objectIndex].material.shininess);
                    intersected += computeSoftShadow(lights[i], raySceneIntersection.rayMeshIntersection.intersection);
                    intersected /= LIGHT_SAMPLES;
                    float shadow = 1.f - intersected;
                    color *= shadow;
                }
                else
                {
                    color += computeHardShadow(ray, 
                                      lights[i], 
                                      raySceneIntersection.rayMeshIntersection.normal, 
                                      raySceneIntersection.rayMeshIntersection.intersection,
                                      meshes[raySceneIntersection.objectIndex].material.ambient_material,
                                      meshes[raySceneIntersection.objectIndex].material.diffuse_material,
                                      meshes[raySceneIntersection.objectIndex].material.specular_material,
                                      meshes[raySceneIntersection.objectIndex].material.shininess
                                     );
                }
                
            }
        }
        return color;

    }

        
    Vec3 rayTrace( Ray const & rayStart ) {
        
        Vec3 color;
        color = rayTraceRecursive(rayStart, 10, 4.9f);
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
            light.pos = Vec3( 0.0, 1.5, 0.0 );
            light.direction = Vec3(0.0, -1.0, 0.0);
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
            s.material.type = Material_Glass;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.shininess = 16;
            s.material.transparency = 1.0;
            s.material.index_medium = 1.0;
        }

       

        { //MIRRORED Sphere
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1.0, -1.25, -0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 0.,1.,0. );
            s.material.specular_material = Vec3(  0.,1.,0. );
            s.material.shininess = 16;
            s.material.transparency = 1.;
            s.material.index_medium = 1.4;
        }

        {
            meshes.resize( meshes.size() + 1);
            Mesh & mesh = meshes[meshes.size() - 1];
            mesh.loadOFF("./data/suzanne.off");
            mesh.centerAndScaleToUnit();
            mesh.build_arrays();
            // mesh.translate(Vec3(1.0, 0.3, 0.5));
            // mesh.build_arrays();
            mesh.material.type = Material_Diffuse_Blinn_Phong;
            mesh.material.diffuse_material = Vec3(0.9,0.5,0.47);
            mesh.material.specular_material = Vec3(0.9,0.5,0.47);
            mesh.material.shininess = 20;
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

    void setup_base_cornell_box(){
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3( 0.0, 1.5, 0.0 );
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
            s.material.diffuse_material = Vec3( 1.,1.,1. );
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
            s.material.diffuse_material = Vec3( 1.0,1.0,1.0 );
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
            s.material.diffuse_material = Vec3( 1.0,1.0,1.0 );
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
            s.material.type = Material_Glass;
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
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,1.,1. );
            s.material.specular_material = Vec3(  1.,1.,1. );
            s.material.shininess = 16;
            s.material.transparency = 1.0;
            s.material.index_medium = 1.4;
        }

    }

    void setup_single_mesh(){
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
            meshes.resize( meshes.size() + 1);
            Mesh & mesh = meshes[meshes.size() - 1];
            mesh.loadOFF("./data/suzanne.off");
            mesh.centerAndScaleToUnit();
            mesh.build_arrays();
            mesh.material.type = Material_Diffuse_Blinn_Phong;
            mesh.material.diffuse_material = Vec3( 0.,1.,0.);
            mesh.material.specular_material = Vec3( 0.2,0.2,0.2 );
            mesh.material.shininess = 20;
        }

    }

};

#endif
