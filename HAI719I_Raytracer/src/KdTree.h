#ifndef KDTREE_H
#define KDTREE_H

#include "Ray.h"
#include "Triangle.h"
#include <algorithm>
#include <cfloat>
#include <stack>
#include <set>


class Node
{
public:
	Node* childLeft, * childRight;
	Vec3 bbmin, bbmax;
	unsigned int checkType; //axis x : 0, axis y : 1, axis z :2;
	int begin, end; // index in Mesh::triangles list	

	Node(Vec3 _bbmin, Vec3 _bbmax, unsigned int _checkType, int _begin, int _end):
	bbmin(_bbmin), bbmax(_bbmax), checkType(_checkType), begin(_begin), end(_end)
	{}
	~Node(){}

	void printNode()
	{
		printf("Node : [%d : %d] axis : %d\nBounding Box \n", begin,end,checkType);
		std::cout<<bbmin<<" | "<<bbmax<<std::endl;
		printf("Triangles : ");
		for (int i = 0; i < triangles.size(); ++i)
		{
			std::cout<<triangles[i]<<" ";
		}
		std::cout<<std::endl;
		
	}

	std::vector<int> triangles;
	void buildNodeTriangles(int nbTriangles)
	{
	    int _end = (end == nbTriangles)? end : end+1;
		for (int i = begin; i < _end; i++)
        {
                triangles.push_back(i);
        }
	}
};

struct RayNodeIntersection
{
	float tmin, tmax;
	bool intersectionExists = false;
	Node* intersected;
};

struct ComparableTriangle
{
	long int id;
	Vec3 centroid;
};

class KdTreeUtils
{
public:
	/*
	* Comparison function to sort triangles depending on the axis and their centroid
	*/
	static bool comp_x(ComparableTriangle &a, ComparableTriangle &b){return a.centroid[0] < b.centroid[0];}
	static bool comp_y(ComparableTriangle &a, ComparableTriangle &b){return a.centroid[1] < b.centroid[1];}
	static bool comp_z(ComparableTriangle &a, ComparableTriangle &b){return a.centroid[2] < b.centroid[2];}

	/*
	* Creates a vector composed by min coordinates of each parameter
	* Used to calculate min bound of bounding box
	*/
	static Vec3 min3(Vec3 &a, Vec3 &b, Vec3 &c)
	{
		Vec3 result;
		for (int i = 0; i < 3; ++i)
		{
			result[i] = std::min(std::min(a[i], b[i]), c[i]);
		}
		return result;
	}

	static Vec3 min2(Vec3 &a, Vec3 &b)
	{
		Vec3 result;
		for (int i = 0; i < 3; ++i)
		{
			result[i] = std::min(a[i], b[i]);
		}
		return result;
	}

	/*
	* Creates a vector composed by max coordinates of each parameter
	* Used to calculate max bound of bounding box
	*/
	static Vec3 max3(Vec3 &a, Vec3 &b, Vec3 &c)
	{
		Vec3 result;
		for (int i = 0; i < 3; ++i)
		{
			result[i] = std::max(std::max(a[i], b[i]), c[i]);
		}
		return result;
	}

	static Vec3 max2(Vec3 &a, Vec3 &b)
	{
		Vec3 result;
		for (int i = 0; i < 3; ++i)
		{
			result[i] = std::max(a[i], b[i]);
		}
		return result;
	}

	static bool intersectNode(Ray const &ray, Node* node)
	{
	    float tmin, tmax, tYmin, tYmax, tZmin, tZmax;
	    Vec3 bbmin = node->bbmin;
	    Vec3 bbmax = node->bbmax;
	    Vec3 origin = ray.origin(); Vec3 dir = ray.direction();
        tmin = (bbmin[0] - origin[0])/dir[0]; // tXmin;
        tmax = (bbmax[0] - origin[0])/dir[0]; // tXmax;
        tYmin = (bbmin[1] - origin[1])/dir[1]; 
        tYmax = (bbmax[1] - origin[1])/dir[1]; 
        tZmin = (bbmin[2] - origin[2])/dir[2]; 
        tZmax = (bbmax[2] - origin[2])/dir[2]; 

        if(tmin > tmax) std::swap(tmin, tmax); // swap value if tmin is not min;
        if(tYmin > tYmax) std::swap(tYmin, tYmax); //same
        if(tZmin > tZmax) std::swap(tZmin, tZmax); //same
        if(tmin > tYmax || tmax < tYmin) return false;
        if(tmin < tYmin) tmin = tYmin;
        if(tmax > tYmax) tmax = tYmax;
        if(tmin > tZmax || tmax < tZmin) return false;
        return true;
	}



};


class KdTree
{
private:

	std::vector<ComparableTriangle> triangles;

public:
	Node* root;
	KdTree(std::vector<Vec3> coordinates, std::vector<unsigned int> indexes, size_t nbTriangles)
	{
		triangles = buildTriangles(coordinates, indexes, nbTriangles);
		root = buildTree(triangles, coordinates, indexes, 0, nbTriangles, 0);
	};
	~KdTree();

	/*
	* Build ComparableTriangle list from meshes info (positions_array, triangles_array) 
	* Theorically it is called only 1 time (in the constructor)
	*/
	std::vector<ComparableTriangle> buildTriangles(std::vector<Vec3> coordinates, std::vector<unsigned int> indexes, size_t nbTriangles)
	{
		std::vector<ComparableTriangle> triangles;
		for (int i = 0; i < nbTriangles; ++i)
		{
			ComparableTriangle triangle;
			triangle.id = i;

			unsigned int v1 = indexes[3*i];
			unsigned int v2 = indexes[3*i + 1];
			unsigned int v3 = indexes[3*i + 2];
			Vec3 c1 = coordinates[v1];
			Vec3 c2 = coordinates[v2];
			Vec3 c3 = coordinates[v3];

			triangle.centroid = (c1+c2+c3)/3.0;

			triangles.push_back(triangle);
		}
		return triangles;
	}

	/*
	* Recursive function to build tree
	* Sort triangles depending on the chosen axis (cycling : 0,1,2) and find median. 
	* The coordinates and indexes lists are only accessed 
	* Begin and end parameters help to recursive aspect of function (instead of give spliced lists)
	*/
	Node * buildTree(std::vector<ComparableTriangle> &triangles, std::vector<Vec3> coordinates, std::vector<unsigned int> indexes, int begin, int end, unsigned int axis)
	{
		if(end-begin < 32) return nullptr;

		switch(axis)
		{
			case 0: //Sort by x-axis
				std::sort(triangles.begin() + begin, triangles.begin() + end, KdTreeUtils::comp_x);
				break;
			case 1: //Sort by y-axis
				std::sort(triangles.begin() + begin, triangles.begin() + end, KdTreeUtils::comp_y);
				break;
			case 2: //Sort by z-axis
				std::sort(triangles.begin() + begin, triangles.begin() + end, KdTreeUtils::comp_z);
				break;
			default:
				printf("Axis error\n");
				break;
		}

		

		int half = (begin + end)/2; // divide triangles
		long int id = triangles[half].id; //median index

		// find BoundingBox 
		unsigned int v1 = indexes[3*id];
		unsigned int v2 = indexes[3*id + 1];
		unsigned int v3 = indexes[3*id + 2];
		Vec3 c1 = coordinates[v1];
		Vec3 c2 = coordinates[v2];
		Vec3 c3 = coordinates[v3];
		Vec3 medianMin = KdTreeUtils::min3(c1, c2, c3);
		Vec3 medianMax = KdTreeUtils::max3(c1, c2, c3);


		Node* median = new Node(medianMin, medianMax, axis, begin, end);
		int nbTriangles = triangles.size();
		median->buildNodeTriangles(nbTriangles);
		

		axis = (axis>=2)? 0 : axis+1; //assure cycle
		median->childLeft  = buildTree(triangles, coordinates, indexes, begin, half, axis); 
		int rightHalf = ((end-begin) == 2) ? half : half+1;
		median->childRight = buildTree(triangles, coordinates, indexes, rightHalf, end, axis);

		// update node bounding box with new bbmin and bbmax 
		updateBoundingBox(median);

		return median;

	}

	void updateBoundingBox(Node* node)
	{
		Vec3 children_bbmin = node->bbmin;
		Vec3 children_bbmax = node->bbmax;
		if(node->childLeft != nullptr && node->childRight != nullptr)
		{
			children_bbmin = KdTreeUtils::min3(node->bbmin, node->childLeft->bbmin, node->childRight->bbmin);
			children_bbmax = KdTreeUtils::max3(node->bbmax, node->childLeft->bbmax, node->childRight->bbmax);
		}
		else 
		{
			if(node->childLeft != nullptr)
			{
				children_bbmin = KdTreeUtils::min2(node->bbmin, node->childLeft->bbmin);
				children_bbmax = KdTreeUtils::max2(node->bbmax, node->childLeft->bbmax);
			}
			else if(node->childRight != nullptr)
			{
				children_bbmin = KdTreeUtils::min2(node->bbmin, node->childRight->bbmin);
				children_bbmax = KdTreeUtils::max2(node->bbmax, node->childRight->bbmax);
			}
		}
		node->bbmin = children_bbmin;
		node->bbmax = children_bbmax;
		
	}


	
	void printTreeFrom(Node* node, int depth)
	{
		if(node == nullptr) 
		{
			return;
		}
		printf("/////////////////// Depth %d /////////////////\n", depth);
		node->printNode();
		if(node->childLeft == nullptr && node->childRight == nullptr) printf("Leaf\n");
		printTreeFrom(node->childLeft, depth+1);
		printTreeFrom(node->childRight, depth+1);
		
	}

	/*
	Iterates over tree and add index of intersected triangle in the given set
	*/
	std::vector<int> intersectRecursive(Ray const& ray, Node* node)
	{
	    // if (node == nullptr)
	    //     return;

	    // if (!KdTreeUtils::intersectNode(ray, node)) return;
	    	

	    // if (node->childLeft == nullptr && node->childRight == nullptr)
	    // {
	    //     // Si c'est une feuille, ajouter les indices qui intersectent le rayon
	    //     int end = (node->end == nbTriangles)? node->end : node->end+1;
	    //     for (int i = node->begin; i < end; ++i)
	    //     {
	    //             indexes.insert(i);
	    //     }
	    // }
	    // else
	    // {

	    //     // Appeler récursivement les fonctions pour les deux enfants
	    //     intersectRecursive(ray, node->childLeft, indexes, nbTriangles);
	    //     intersectRecursive(ray, node->childRight, indexes, nbTriangles);
	    // }
	     //Si il a un fils 1
	    // if (node->childLeft && KdTreeUtils::intersectNode(ray, node->childLeft))
	    // {
	    //     // Si il a aussi un fils 2 
	    //     if (node->childRight && KdTreeUtils::intersectNode(ray, node->childRight))
	    //     {
	    //         // Alors on parcours recursivement encore fils 1 et 2 et return la concaténation
	    //         intersectRecursive(ray, node->childLeft, indexes, nbTriangles);
	    //         intersectRecursive(ray, node->childRight, indexes, nbTriangles);
	    //         return;
	    //     }
	    //     else
	    //     {
	    //     //sinon -> si il n'a pas de fils 2 alors on parcourt recursivement seulement fils 1 
	    //         intersectRecursive(ray, node->childLeft, indexes, nbTriangles);
	    //     	return;
	    //     }
	    // }
	    // // Si il n'a que fils 2 alors on parcourt recursivement fils 2  seulement
	    // else if (node->childRight && KdTreeUtils::intersectNode(ray, node->childRight))
	    // {
	    //     intersectRecursive(ray, node->childRight, indexes, nbTriangles);
	    //     return;
	    // }
	    // // On retourne la liste totale des triangles
	    // int end = (node->end == nbTriangles)? node->end : node->end+1;
     //    for (int i = node->begin; i < end; ++i)
     //    {
     //            indexes.insert(i);
     //    }

		//Si il a un fils 1
	    if (node->childLeft && KdTreeUtils::intersectNode(ray, node->childLeft))
	    {
	        // Si il a aussi un fils 2 
	        if (node->childRight && KdTreeUtils::intersectNode(ray, node->childRight))
	        {
	            // Alors on parcours recursivement encore fils 1 et 2 et return la concaténation
	            std::vector<int> result = intersectRecursive(ray, node->childLeft);
	            std::vector<int> second = intersectRecursive(ray, node->childRight);
	            //Result = result + second;
	            result.insert(result.end(), second.begin(), second.end());
	            return result;
	        }
	        else
	        //sinon -> si il n'a pas de fils 2 alors on parcourt recursivement seulement fils 1 
	            return intersectRecursive(ray, node->childLeft);
	    }
	    // Si il n'a que fils 2 alors on parcourt recursivement fils 2  seulement
	    if (node->childRight && KdTreeUtils::intersectNode(ray, node->childRight))
	    {
	        return intersectRecursive(ray, node->childRight);
	    }
	    // On retourne la liste totale des triangles
	    return node->triangles;

	}

};


#endif