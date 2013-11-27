#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>

#include <cstdlib>
#include <string.h>
#include <time.h>
#include <math.h>
#include "primitives.h"


//Constructor
Point::Point(double nx, double ny, double nz) {
		x = nx;
		y = ny;
		z = nz;
}

//Constructor
Vector::Vector(double nx, double ny, double nz){
        dx = nx;
        dy =  ny; 
        dz = nz; 
        mag = sqrt(nx*nx + ny*ny + nz*nz);
}
Vector::Vector(Point p){
       dx = p.x;
       dy =  p.y; 
       dz = p.z; 
       mag = sqrt(dx*dx + dy*dy + dz*dz);
}


//Construct vector from two points
//Vector goes from p1 to p2
Vector::Vector(Point p1, Point p2){
    dx = p2.x - p1.x; 
    dy = p2.y - p1.y;
    dz = p2.z - p1.z; 
    mag = sqrt(dx*dx + dy*dy + dz*dz);
}


//Construct ray given starting coordinates and end coordinates
Ray:: Ray(double nx, double ny, double nz, double ndx, double ndy, double ndz) {
		origin = Point (nx, ny, nz);
        direction = Vector(ndx, ndy, ndz).normalize();
}
//Construct ray given an origin and direction vector
Ray::Ray(Point p, Vector dir){
    origin = p;
    direction = dir;     
    
}




//Construct ray given two points
Ray:: Ray(Point startingPoint, Point newPoint) {
		origin = Point(startingPoint.x,
                startingPoint.y,
                startingPoint.z);

		direction = Vector(newPoint.x - startingPoint.x,
                    newPoint.y - startingPoint.y,
                    newPoint.z - startingPoint.z).normalize(); 

}


//Construct a Color
 Color::Color(double red, double green, double blue){
  
    r = red;
    g = green;
    b = blue;
  
 }
 
 
 
 
 //Light constructors
Light::Light(double xp, double yp, double zp, double red, double green, double blue, bool dir){
    x = xp;
    y = yp;
    z = zp; 
    r = red;
    g = green;
    b = blue; 
    source = Point(x, y, z); 
    color = Color(red, green, blue); 
    directional = dir; 
    
    
}   



//****************************************************
// Point functions
//****************************************************



//PointA - PointB = vector from B to A
//return Vector that is result of this - p, which is from p to this
Vector Point:: subtract (Point p){
    return Vector(x - p.x, y - p.y, z - p.z); 
}

Vector Point::operator- (Point p){
    return Vector(x - p.x, y - p.y, z - p.z); 
}

Point Point::operator+ (Vector v){
    return Point(x + v.dx, y + v.dy, z + v.dz); 
}

// implemented only for the sake of camera; be wary of using normally....
Point Point::operator+ (Point p){
    return Point(x + p.x, y + p.y, z + p.z); 
}

Point Point::operator* (double s){
    return Point(x*s, y*s, z*s); 
}



//****************************************************
// Vector functions
//****************************************************

//normalize this vector
Vector Vector::normalize(){
    Vector v = Vector((double)dx / mag, (double)dy / mag, (double)dz / mag);  
    v.mag = 1;
    return v;
}
double Vector::dotProduct(Vector v){
    return dx * v.dx + dy * v.dy + dz * v.dz; 
}

//return negation of this vector
Vector Vector::negative(){
   return Vector(-dx, -dy, -dz); 
}

//multiply this vector by scalar
Vector Vector :: operator*(double a){
   return Vector(a*dx, a*dy, a*dz); 
}

//add to this vector to another vector
Vector Vector::operator+ (Vector v){
    return Vector(dx + v.dx, dy + v.dy, dz + v.dz); 
}
//subtract by another vector
Vector Vector::operator- (Vector v){
    return Vector(dx - v.dx, dy - v.dy, dz - v.dz); 
}


// takes the cross product of this vector and vector V.  Returns 
// in a vector which is perpendicular to both and therefore normal 
// to the plane containing them
Vector Vector::crossProduct(Vector v){
    return Vector(dy*v.dz - dz*v.dy,
				  dz*v.dx - dx*v.dz,
				  dx*v.dy - dy*v.dx);
   
    
}

//****************************************************
// Light functions
//****************************************************
 	
	void Light::initPos(double nx, double ny, double nz) {
		x = nx;
		y = ny;
		z = nz;
	}
	void Light::initRGB(double nr, double ng, double nb) {
		r = nr;
		g = ng;
		b = nb;
	}

//****************************************************
// Color functions
//****************************************************
//add the rgb values of this color to another color 	
Color Color::operator+(Color color){
    return Color(r + color.r, g + color.g, b + color.b); 
}

void Color::operator+=(Color color){
    r += color.r;
    g += color.g;
    b += color.b;
}
    
//multiply the rgb values of this color to another color 	
Color Color::operator*(Color color){
    return Color(r * color.r, g * color.g, b * color.b); 
}
//overload to allow scalar multiplication
Color Color::operator*(double a){
    return Color(r * a, g * a, b * a); 
}

Color Color::clone(){
    return Color(r, g, b); 
}

//****************************************************
//Matrix functions
//****************************************************

//default constructor makes identity matrix
Matrix::Matrix(){
	for(int i = 0; i < 4; i++)
		for(int j = 0; j < 4; j++)
			m[i][j] = 0;

    m[0][0] = 1;   
    m[1][1] = 1;  
    m[2][2] = 1;  
    m[3][3] = 1;  
    
    
}

//Create a 4x4 Matrix. t: translation, r: rotate, s: scale
//Store Matrix in row major form. 
Matrix::Matrix(char type, double x, double y, double z, double angle = 0){
	
	for(int i = 0; i < 4; i++)
		for(int j = 0; j < 4; j++)
			m[i][j] = 0;
    
    if(type == 't'){ 
		m[0][0] = 1;   
		m[1][1] = 1;  
		m[2][2] = 1;  
		m[3][3] = 1;  
		
		m[0][3] = x;   
		m[1][3] = y;  
		m[2][3] = z;  
    }
    //Reference for making rotation matrix: http://www.talisman.org/opengl-1.1/Reference/glRotate.html
    else if(type == 'r'){
        Vector n = Vector(x, y, z).normalize();
        x = n.dx, y = n.dy, z = n.dz; 
        double c = cos(angle); 
        double s = sin(angle); 
        
        /*m[0][0] = x * x * (1 - c) + c;
        m[0][1] = x * y * (1 - c) + z * s;
        m[0][2] = x * z * (1 - c) + y * s; 
    
        m[1][0] = y * z * (1 - c) + z * s;
        m[1][1] = y * y * (1 - c) + c;
        m[1][2] = y * z * (1 - c) - x * s; 
 
        m[2][0] = z * x * (1 - c) - y * s; 
        m[2][1] = z * y * (1 - c) + x * s;
        m[2][2] = z * z * (1 -c ) + c;*/
        
        m[0][0] = x * x * (1 - c) + c;
        m[0][1] = x * y * (1 - c) + z * s;
        m[0][2] = x * z * (1 - c) - y * s; 
    
        m[1][0] = y * x * (1 - c) - z * s;
        m[1][1] = y * y * (1 - c) + c;
        m[1][2] = y * z * (1 - c) + x * s; 
 
        m[2][0] = z * x * (1 - c) + y * s; 
        m[2][1] = z * y * (1 - c) - x * s;
        m[2][2] = z * z * (1 - c) + c;
    
        m[3][3] = 1; 

    } else if(type == 's'){
		m[0][0] = x;   
		m[1][1] = y;  
		m[2][2] = z;  
		m[3][3] = 1;  
    }
    
} 

// Naive Matrix Multiplication;  O(n^3).  Order matters!!! M*B != B*M
Matrix Matrix::operator*(Matrix B){
	Matrix C;
	double sum;
    for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			sum = 0;
			for (int k = 0; k < 4; k++) {
				sum += m[i][k] * B.m[k][j];
			}
			C.m[i][j] = sum;
		}
		
	}
	
	return C;
}

Point Matrix::operator*(Point p){
	double c[4] = {0};
	double pv[4] = {p.x, p.y, p.z, 1};
	
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			c[i] += m[i][j] * pv[j];
		}
	}
	
	return Point(c[0], c[1], c[2]);
}
Vector Matrix::operator*(Vector v){
	double c[4] = {0};
	double pv[4] = {v.dx, v.dy, v.dz, 0};
	
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			c[i] += m[i][j] * pv[j];
		}
	}
	
	return Vector(c[0], c[1], c[2]);//.normalize();
}

Ray Matrix::operator*(Ray r){
	Ray nray;
	nray.origin = (*this) * (r.origin);
	
	nray.direction = (*this) * (r.direction);
	return nray;
}


// vector * matrix.... not the same as M * v
Vector Matrix :: vectorTimesM(Vector v){
   return Point(v.dx*m[0][0]  +  v.dy*m[0][1]  +  v.dz*m[0][2],
				v.dx*m[1][0]  +  v.dy*m[1][1]  +  v.dz*m[1][2],
				v.dx*m[2][0]  +  v.dy*m[2][1]  +  v.dz*m[2][2]);
				//x*M.m[0][0]  +  y*M.m[0][1]  +  z*M.m[0][2] + w,  w term ignored
}

// returns a new matrix that is the inversion of m.
// modified from http://stackoverflow.com/questions/2624422/efficient-4x4-matrix-inverse-affine-transform/7596981#7596981
Matrix Matrix::invert() {
	Matrix inv;
	
	double s0 = m[0][0] * m[1][1] - m[1][0] * m[0][1];
    double s1 = m[0][0] * m[1][2] - m[1][0] * m[0][2];
    double s2 = m[0][0] * m[1][3] - m[1][0] * m[0][3];
    double s3 = m[0][1] * m[1][2] - m[1][1] * m[0][2];
    double s4 = m[0][1] * m[1][3] - m[1][1] * m[0][3];
    double s5 = m[0][2] * m[1][3] - m[1][2] * m[0][3];
    
    double c5 = m[2][2] * m[3][3] - m[3][2] * m[2][3];
    double c4 = m[2][1] * m[3][3] - m[3][1] * m[2][3];
    double c3 = m[2][1] * m[3][2] - m[3][1] * m[2][2];
    double c2 = m[2][0] * m[3][3] - m[3][0] * m[2][3];
    double c1 = m[2][0] * m[3][2] - m[3][0] * m[2][2];
    double c0 = m[2][0] * m[3][1] - m[3][0] * m[2][1];
	
	
    double det = s0*c5 - s1*c4 + s2*c3 + s3*c2 - s4*c1 + s5*c0;
    double invdet = 1.0 / det;

    inv.m[0][0] = ( m[1][1] * c5 - m[1][2] * c4 + m[1][3] * c3) * invdet;
    inv.m[0][1] = (-m[0][1] * c5 + m[0][2] * c4 - m[0][3] * c3) * invdet;
    inv.m[0][2] = ( m[3][1] * s5 - m[3][2] * s4 + m[3][3] * s3) * invdet;
    inv.m[0][3] = (-m[2][1] * s5 + m[2][2] * s4 - m[2][3] * s3) * invdet;

    inv.m[1][0] = (-m[1][0] * c5 + m[1][2] * c2 - m[1][3] * c1) * invdet;
    inv.m[1][1] = ( m[0][0] * c5 - m[0][2] * c2 + m[0][3] * c1) * invdet;
    inv.m[1][2] = (-m[3][0] * s5 + m[3][2] * s2 - m[3][3] * s1) * invdet;
    inv.m[1][3] = ( m[2][0] * s5 - m[2][2] * s2 + m[2][3] * s1) * invdet;

    inv.m[2][0] = ( m[1][0] * c4 - m[1][1] * c2 + m[1][3] * c0) * invdet;
    inv.m[2][1] = (-m[0][0] * c4 + m[0][1] * c2 - m[0][3] * c0) * invdet;
    inv.m[2][2] = ( m[3][0] * s4 - m[3][1] * s2 + m[3][3] * s0) * invdet;
    inv.m[2][3] = (-m[2][0] * s4 + m[2][1] * s2 - m[2][3] * s0) * invdet;

    inv.m[3][0] = (-m[1][0] * c3 + m[1][1] * c1 - m[1][2] * c0) * invdet;
    inv.m[3][1] = ( m[0][0] * c3 - m[0][1] * c1 + m[0][2] * c0) * invdet;
    inv.m[3][2] = (-m[3][0] * s3 + m[3][1] * s1 - m[3][2] * s0) * invdet;
    inv.m[3][3] = ( m[2][0] * s3 - m[2][1] * s1 + m[2][2] * s0) * invdet;

	return inv;
}

Matrix Matrix::clone() {
	Matrix B;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			B.m[i][j] = m[i][j];
		}
		
	}
	return B;
}

void Matrix::debug() {
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << "[" << m[i][j] << "]";
		}
		cout << endl;
	}
}



//****************************************************
// MatrixStack functions
//****************************************************

//default: initialize stack to only have identity matrix; 
MatrixStack::MatrixStack(){
    stack.push_back(Matrix());
    stackT.push_back(Matrix());
    stackR.push_back(Matrix());
    stackS.push_back(Matrix());
    product = Matrix();
    productT = Matrix();
    productS = Matrix();
    productR = Matrix();
}

// push the top to a new matrix product;  inits it to Identity.  Updates
// are added throigh addTransform.
void MatrixStack::push(){
	stack.push_back(Matrix());
    stackT.push_back(Matrix());
    stackR.push_back(Matrix());
    stackS.push_back(Matrix());
    
}

void MatrixStack::pop(){
    stack.pop_back(); 
    stackT.pop_back(); 
    stackR.pop_back(); 
    stackS.pop_back(); 
    
    product = Matrix();
	productT = Matrix();
	productR = Matrix();
	productS = Matrix();

	size_t stackSize = stack.size();
	for (size_t i = 0; i < stackSize; i++) {
		product = product * stack[i];
		productT = productT * stackT[i];
		productR = productR * stackR[i];
		productS = productS * stackS[i];
		
	}
	//product = productT * productR * productS;
}

// mulitiplies in the new transformation matrix B to the top of the stack 
// and updates product to reflect new stack.
void MatrixStack::addTransform(Matrix B, char type) {
	
	if(type == 't'){ 
		stackT[stackT.size() - 1] = stackT[stackT.size() - 1] * B;
		productT = productT * B;
	} else if(type == 'r'){ 
		stackR[stackR.size() - 1] = stackR[stackR.size() - 1] * B;
		productR = productR * B;
	} else if(type == 's'){ 
		stackS[stackS.size() - 1] = stackS[stackS.size() - 1] * B;
		productS = productS * B;
	}
	stack[stack.size() - 1] = stack[stack.size() - 1] * B;
	product = product * B;
	
	
	
	//cout << "new product:" << endl;
	//productR.debug();
}


