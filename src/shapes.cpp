#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>

#include <cstdlib>

#include <time.h>
#include <math.h>

#include "shapes.h"
/*
BoundingBox::BoundingBox(float minx, float maxx, float miny, float maxy, float minz, float maxz){
    min_x = minx, max_x = maxx, min_y = miny, max_y = maxy, min_z = minz, max_z = maxz; 
}*/
/*
int BoundingBox::getLongestAxis(){
    float xd = max_x - min_x; 
    float yd = max_y - min_y;
    float zd = max_z - min_z;
    float m = max(xd, max(yd, zd));
    if(m == xd)
        return X_AXIS;
    else if (m == yd)
        return Y_AXIS;
    return Z_AXIS; 
    
}*/
/*
float BoundingBox::getMidPoint(int axis){
    if(axis == X_AXIS){
        
        return (max_x - min_x)/2; 
    }
    if(axis == Y_AXIS){
        
        return (max_y - min_y)/2;  
    }
    else{
        
         return (max_z - min_z)/2; 
    }
    
}*/
/*
BoundingBox Triangle::getBB(){
    //don't return bounding box with 0 volume
    
   float xmin = min(min(v1.x, v2.x), v2.x);
   float xmax = max(max(v1.x, v2.x), v2.x);
   float ymin = min(min(v1.y, v2.y), v2.y);
   float ymax = max(max(v1.y, v2.y), v2.y);
   float zmin = min(min(v1.z, v2.z), v2.z);
   float zmax = max(max(v1.z, v2.z), v2.z);
   
   if (xmin == xmax) xmax += 1;
   if (ymin == ymax) ymax += 1;
   if (zmin == zmax) zmax += 1;
    
    return BoundingBox(xmin, xmax,
                       ymin, ymax,
                       zmin, zmax); 
                       
}
*/

/*
ShapeList::ShapeList(vector<Shape*> shapes){
    allShapes = shapes; 
}*/



/*
BoundingBox ShapeList::getRootBox(){
    float minx, miny, minz = LARGE_NUM; 
    float maxx, maxy, maxz = -LARGE_NUM; 
    BoundingBox bb; 
    
	size_t shapeCount = allShapes.size();
    for(size_t i = 0; i < shapeCount; i++){
        bb = allShapes[i]->getBB(); 
        minx = fmin(bb.min_x, minx);
        miny = fmin(bb.min_y, miny);
        minz = fmin(bb.min_z, minz);
        
        maxx = fmax(bb.max_x, maxx);
        maxy = fmax(bb.max_y, maxy);
        maxz = fmax(bb.max_z, maxz);
    }
    
    return BoundingBox(minx, maxx, miny, maxy, minz, maxz); 
    
}*/
/*
bool BoundingBox::intersect(Ray r){ 
    float xmin, xmax, ymin, ymax, zmin, zmax; 
    float xd = r.direction.dx;
    if (xd == 0) xd = 0.00001;
    
    float yd = r.direction.dy;
    if (yd == 0) yd = 0.00001;
    
    float zd = r.direction.dz;
    if (zd == 0) zd = 0.00001;
    
    float a_x = 1 / xd;
    if(a_x >= 0){
        xmin = a_x * (min_x - r.origin.x);
        xmax = a_x * (max_x - r.origin.x); 
    }
    else{
        xmin = a_x * (max_x - r.origin.x);
        xmax = a_x * (min_x - r.origin.x); 
        
    }
    
    float a_y = 1 / yd;
    if(a_y >= 0){
        ymin = a_y * (min_y - r.origin.y);
        ymax = a_y * (max_y - r.origin.y); 
    }
    else{
        ymin = a_y * (max_y - r.origin.y);
        ymax = a_y * (min_y - r.origin.y); 
        
    }
   
   
   float a_z = 1 / zd;
    if(a_x >= 0){
        zmin = a_z * (min_z - r.origin.z);
        zmax = a_z * (max_z - r.origin.z); 
    }
    else{
        zmin = a_z * (max_z - r.origin.z);
        zmax = a_z * (min_z - r.origin.z); 
        
    }
    
    if((xmin > ymax) || (xmin > zmax) || (ymin > xmax) || (ymin > zmax) || (zmin > xmax) || (zmin > ymax))
        return false;
    return true; 
   
   
}*/
/*
//make tree by passing in allshapes to get root
//adapted from http://www.flipcode.com/archives/Dirtypunks_Column-Issue_05_AABB_Trees_Back_To_Playing_With_Blocks.shtml
AABB_Node::AABB_Node(ShapeList sl, int depth){ 
    
    
    bb = sl.getRootBox(); 
    
    if(depth < 2){
        
        int axis = bb.getLongestAxis(); 
        Shape* s; 
        ShapeList shapeBucket[2]; 
        
        float mid = bb.getMidPoint(axis); 
        //divide this bb by midpoint along longest axis, sort objects
		size_t shapeCount = sl.allShapes.size();
		for(size_t i = 0; i < shapeCount; i++){
            if(sl.allShapes[i]->getBB().getMidPoint(axis) < mid)
                shapeBucket[0].allShapes.push_back(sl.allShapes[i]);
            else
                shapeBucket[1].allShapes.push_back(sl.allShapes[i]); 
        }
        
        children[0] = new AABB_Node(shapeBucket[0], depth + 1);
        children[1] = new AABB_Node(shapeBucket[1], depth + 1); 
    }
    else{ //make leaf node 
        children[0] = 0; 
        containedShapes = sl; 
        
        
    }
    
}*/
/*
bool AABB_Node::CollisionTest(Ray ray, Point* p, Shape*& sh){
    if(!bb.intersect(ray)){
        return false; 
    }

        if(children[0]){
            if(children[0]->CollisionTest(ray, p, sh))
                return true;
            if(children[1]->CollisionTest(ray, p, sh))
                return true; 
            return false; 
        
        }
        //do intersection with shapes at this node 
        return containedShapes.checkIntersect(ray, p, sh, LARGE_NUM); 
    
    return false; 
}*/
