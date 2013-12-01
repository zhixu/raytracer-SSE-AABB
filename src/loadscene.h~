#ifndef LOADSCENE_H
#define LOADSCENE_H

#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>

#include <nmmintrin.h>

#include "primitives.h"

class parsedScene{
    public:		
        int width, height, reflectiondepth, triangleCount, lightCount; 
        float fov;        
        float aaJitter; // jitter distance, 0 means none, 16 means approx 16 pixels???? not really, but it looks like it
						// jitter does not work unless aaDepth >= 2.0
		float aaDepth; // how many more iterations for aa, 1 means none, 32 is really strong, shoudl be d^k

        Color ambient;
        Vector3 lookfrom, lookat, up_dir;

        std::string outputFileName;
        
        BoundingBox* bbList; //array of Bounding Boxes corresponding to each triangle
        __m128 (*triangleList)[9];
		//Brdf (*brdfList);
        __m128 (*lightList)[3];
        
        AABB_Node root; 
    
    parsedScene(); 
    
};// class parsedScene

parsedScene loadScene(std::string file);
#endif
