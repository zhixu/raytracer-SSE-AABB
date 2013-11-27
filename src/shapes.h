/*
 * shapes
 * moves to primitives.h
*/
#include <iostream>
#include <math.h>
#include <vector>
#include "primitives.h"

#define X_AXIS 0
#define Y_AXIS 1
#define Z_AXIS 2

#define LARGE_NUM 1000000

/*
class BoundingBox{
  public:
    float min_x, max_x;
    float min_y, max_y;
    float min_z, max_z; 
    
    BoundingBox(){};
    BoundingBox(float, float, float, float, float, float); 
    
    int getLongestAxis(); 
    float getMidPoint(int); 
    bool intersect(Ray); 
};*/

/*
class ShapeList{ //see AggregatePrimitive on website 
    public:
        vector<Shape*> allShapes; 
        bool checkIntersect(Ray, Point*, Shape*&, double);
        bool checkIntersect(Ray, double);
        
        BoundingBox getRootBox(); 
        
        ShapeList() {}; 
        ShapeList(vector<Shape*>); 
};*/
/*
class AABB_Node{
  public:
        BoundingBox bb;
        AABB_Node* children[2];  
        ShapeList containedShapes; 
        
        AABB_Node(ShapeList, int); 
        bool CollisionTest(Ray, Point*, Shape*&); 
};*/