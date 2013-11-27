//****************************************************
// primitives.h
// Contains:
// - Point
// - Ray
// - Vector
// - Light
// - Color
// - Camera 
// - MatrixStack
// - Matrix 
//****************************************************

#ifndef LIGHTS_H
#define LIGHTS_H

#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;


//****************************************************
// CLASS DECLARATIONS AND CONSTRUCTORS 
//****************************************************
class Vector; 
class Point{ 
public:
    Point() {};
    Point(double, double, double); 

    double x, y, z;
    Vector subtract(Point); 
    Point operator + (Vector); //position + direction = position 
    Point operator + (Point); // point addition; shortcut for camera's ray formula
    Point operator * (double); 
    Vector operator - (Point);  
    
    
    
};

    
    
class Vector{
 public:
    Vector() {};  
    Vector(double, double, double);
    Vector(Point); 
    Vector(Point, Point); 
    
    double dx, dy, dz;
    double mag;
    
    Vector normalize();
    double dotProduct(Vector);
    Vector negative();
    Vector crossProduct(Vector);
    
    //allow for scalar * and vector +
    Vector operator * (double);
    Vector operator + (Vector); 
    Vector operator - (Vector); 
    
    
};




class Ray{
public:
    Ray() {}; 
    Ray(double, double, double, double, double, double);
    Ray(Point, Point); 
    Ray(Point, Vector);

    Vector direction; 
    Point origin; 
   
    
};


/*
 *  Color class keeps track of r g b values, which are [0, 1]
 */ 
 class Color{
     public:
        Color() {};
        Color(double, double, double);
        
        double r, g, b;
        
        Color operator + (Color);  
        Color operator * (Color); 
        Color operator * (double); 
        void operator += (Color); 
        
        Color clone();
 };

/*
 *  Light class describes a light source that extends in all directions with rgb value from point x y z
 */ 

class Light {
  public:
    Light() {}; 
    Light (double, double, double, double, double, double, bool); 

    void initPos(double, double, double);
    void initRGB(double, double, double); 
    Color color; 
    Point source; 
    double x, y, z, r, g, b;
	bool directional;

};

class Matrix{
    public:
		double m[4][4]; 
    Matrix();
    
    // Instantiator shortcut for translation, scaling, or rotating (depends upon input char)
    Matrix(char, double, double, double, double); 
    
    Matrix operator * (Matrix); 
    Point operator * (Point);
    Vector operator * (Vector);
    Ray operator * (Ray);
    
    Matrix invert();
    void debug();
    Matrix clone();
    
    Vector vectorTimesM(Vector); 
};

class MatrixStack{
    public:
    MatrixStack(); 
    
    vector<Matrix> stack;
    vector<Matrix> stackT;
    vector<Matrix> stackS;
    vector<Matrix> stackR;
    Matrix product;
    Matrix productT; // product of all the T matrices
    Matrix productS;
    Matrix productR;
    
    void push(); 
    void pop();
    void addTransform(Matrix, char);
     
};







 
#endif
