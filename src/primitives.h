//****************************************************
// primitives.h
// Contains:
// - Point
// - Ray
// - Vector
// - Light
// - Color
// - Camera 
// - Triangle
// - AABB Tree?
//****************************************************

#ifndef LIGHTS_H
#define LIGHTS_H

#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <nmmintrin.h>
using namespace std;

#define LARGE_NUM 1000000

/* Type */
//index
#define RAY_ORIGIN_IDX 0
#define RAY_DIRECTION_IDX 1
#define LIGHT_SRC_IDX 0
#define LIGHT_COLOR_IDX 1
#define LIGHT_ISDIR_IDX 2
#define TRIANGLE_PT1_IDX 0
#define TRIANGLE_PT2_IDX 1
#define TRIANGLE_PT3_IDX 2
#define TRIANGLE_NORMAL_IDX 3
#define BRDF_KD_IDX 4
#define BRDF_KS_IDX 5
#define BRDF_KE_IDX 6
#define BRDF_KR_IDX 7
#define BRDF_SP_IDX 8

#define NO_INTERSECTION FLT_MAX
#define X_AXIS 0
#define Y_AXIS 1
#define Z_AXIS 2

//type
typedef __m128 Vector3; // [x, y, z, dummy]
typedef const __m128 constVector3;
typedef __m128 Color; // [r, g, b, dummy]
typedef const __m128 constColor;
typedef __m128 Ray[2]; // [Vector3[origin], Vector3[direction]]
typedef const __m128 constRay[2];
typedef __m128 Light[3]; // [Vector3[source] , Color[color], isDirectional, dummy, dummy, dummy]
typedef const __m128 constLight[3];
//typedef __m128 Brdf[5]; // [Color[kd], Color[ks], Color[ke], kr, sp, dummy, dummy]
//typedef const __m128 constBrdf[5];
typedef __m128 Triangle[9]; // [Vector3[pt1], Vector3[pt2], Vector3[pt3], Vector3[normal], Color[kd], Color[ks], Color[ke], kr, sp]
typedef const __m128 constTriangle[9];

class Brdf {
    public:
        //Color brdf[4];
        Color kd;
        Color ks;
        Color ke;
        float kr, sp;
        Brdf();
        Brdf(Color tkd, Color tks, Color tke, float tkr, float tsp);
};

/* utility function */
//vector
void setVector3(Vector3 &dst, float x, float y, float z);
void vector3Add(Vector3 &dst, constVector3 v1, constVector3 v2);
void vector3Sub(Vector3 &dst, constVector3 v1, constVector3 v2); // v1 - v2
void vector3Scale(Vector3 &dst, constVector3 v1, const float scalar);
float vector3Dot(constVector3 v1, constVector3 v2);
void vector3Cross(Vector3 &dst, constVector3 v1, constVector3 v2); // v1 x v2
void vector3Normalize(Vector3 &dst, constVector3 v1);
void vector3Copy(Vector3 &dst, constVector3 src);
void colorMultiply(Color &dst, constColor c1, constColor c2);
void printVector(Vector3 v);
void printTriangle(Triangle t);
void printRay(Ray r);

//triangle
void setTriangle(Triangle triangle, constVector3 p1, constVector3 p2, constVector3 p3, Brdf brdf);
bool intersectTriangle(constTriangle triangle, constRay ray, float &rr); // check if intersect with the triangle
bool hasIntersect(__m128(*triangleList)[9], const int triangleCount, constRay ray, const float &tmax); // only check if intersect with any triangle
bool nearestIntersect(__m128(*triangleList)[9], const int triangleCount, constRay ray, Vector3 &intersectPt, int &triangleIdx, float &tmax); // check every possible triangle and return the nearest

//ray
void setRayByValue(Ray &ray, float x, float y, float z, float dx, float dy, float dz);
void setRayByPoint(Ray &ray, constVector3 origin, constVector3 p);
void setRayByVector(Ray &ray, constVector3 origin, constVector3 direction);
void getReflection(Vector3 &reflectionDirection, constVector3 rayDirection, constVector3 normal); // get reflection ray direction, normal should be normalized vector

//light
void setLight(Light light, float x, float y, float z, float r, float g, float b, bool isDir);

//brdf
//void brdfCopy(Triangle &dst, Brdf src);

//aabb tree classes

class BoundingBox{
  public:
    float min_x, max_x;
    float min_y, max_y;
    float min_z, max_z; 
    
    BoundingBox(); 
    BoundingBox(constTriangle); 
    BoundingBox(float, float, float, float, float, float); 
    
    int getLongestAxis(); 
    float getMidPoint(int); 
    bool intersect(constRay); 
};
class AABB_Node{
  public:
        BoundingBox bb;
        AABB_Node* children[2];  
        int triCount; 
        __m128 (*containedTriangles)[9];  
        
        AABB_Node(){}; 
        AABB_Node(__m128 (*triList)[9], int, int); 
        float CollisionTest(constRay, Triangle, float*); 
        bool CollisionTest(constRay, float*); 
};






















#endif
