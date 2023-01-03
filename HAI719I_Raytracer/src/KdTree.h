#ifndef KDTREE_H
#define KDTREE_H

#include <Vec3.h>
#include <algorithm>


struct Node
{
	Node* childLeft, childRight;
	Vec3 bbmin, bbmax;
	unsigned int checkType; //axis x : 0, axis y : 1, axis z :2;
	int begin, end; // index in Mesh::triangles_array list
};

struct ComparableTriangle
{
	long int id;
	Vec3 centroid;
};

class KdTreeUtils
{
	/*
	* Comparison function to sort triangles depending on the axis and their centroid
	*/
	static bool comp_x(ComparableTriangle &a, ComparableTriangle &b){return a->centroid[0] < b->centroid[0];}
	static bool comp_y(ComparableTriangle &a, ComparableTriangle &b){return a->centroid[1] < b->centroid[1];}
	static bool comp_z(ComparableTriangle &a, ComparableTriangle &b){return a->centroid[2] < b->centroid[2];}

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

	static Vec3 min2(Vec3 &a, Vec3 &b)
	{
		Vec3 result;
		for (int i = 0; i < 3; ++i)
		{
			result[i] = std::max(a[i], b[i]);
		}
		return result;
	}

};


class KdTree
{
private:
	Node* root;
	std::vector<ComparableTriangle> triangles;

public:
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
		if(end-begin <= 0) return nullptr;

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

		Node* median;
		median->checkType = axis;
		median->bbmin = medianMin;
		median->bbmax = medianMax;
		median->begin = begin;
		median->end = end;

		axis = (axis>=2)? 0 : axis+1; //assure cycle
		median->childLeft  = buildTree(triangles, coordinates, indexes, begin, half, axis);
		median->childRight = buildTree(triangles, coordinates, indexes, half+1, end, axis);

		// update node bounding box with new bbmin and bbmax 

		return median;

	}

	void updateBoundingBox(Node* node)
	{
		Vec3 children_bbmin, children_bbmax;
		if(node->childLeft != nullptr && node->childRight != nullptr)
		{
			children_bbmin = KdTreeUtils::min3(node->bbmin, node->childLeft->bbmin, node->childRight->bbmin);
			children_bbmax = KdTreeUtils::max3(node->bbmax, node->childLeft->bbmax, node->childRight->bbmax);
		}
		else if(node->childLeft != nullptr)
		{
			children_bbmin = KdTreeUtils::min2(node->bbmin, node->childLeft->bbmin);
			children_bbmax = KdTreeUtils::max2(node->bbmax, node->childLeft->bbmax);
		}
		else
		{
			children_bbmin = KdTreeUtils::min2(node->bbmin, node->childRight->bbmin);
			children_bbmax = KdTreeUtils::max2(node->bbmax, node->childRight->bbmax);
		}

		node->bbmin = children_bbmin;
		node->bbmax = children_bbmax;
	}
};


// TODO ITERATOR


#endif