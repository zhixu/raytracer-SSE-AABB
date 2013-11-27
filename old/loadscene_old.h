#ifndef LOADSCENE_H
#define LOADSCENE_H

#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <cstdlib>


#include "primitives.h"
#include "shapes.h"

class parsedScene{
    public:
        int width, height, reflectiondepth; 
        Color ambient;
        Point lookfrom, lookat;
        Vector up_dir;
        double fov;
        
        double aaJitter; // jitter distance, 0 means none, 16 means approx 16 pixels???? not really, but it looks like it
        // jitter does not work unless aaDepth >= 2.0
		double aaDepth; // how many more iterations for aa, 1 means none, 32 is really strong, shoudl be d^k
        
        std::string outputFileName;
        
        ShapeList shapes;
        vector<Light> lights;
    
    parsedScene(); 
    
};


parsedScene loadScene(std::string file);


#endif
