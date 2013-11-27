#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <cstdlib>

#include <omp.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#include <cstdio>
#include <cstdlib>
#include <sys/time.h>
#include <time.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

#ifdef OSX
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/glu.h>
#endif

#include <time.h>
#include <math.h>

#include "src/brdf.h"
#include "src/loadscene.h"
#include "src/imgwriter/lodepng.h"

double timestamp()
{
	struct timeval tv;
	gettimeofday(&tv, 0);
	return tv.tv_sec+1e-6*tv.tv_usec;
}
#define PI 3.14159265  // Should be used from mathlib
inline float sqr(float x) { return x*x; }

using namespace std;

//****************************************************
// Some Classes
//****************************************************

class Viewport;

class Viewport {
  public:
    int w, h; // width and height
};



//****************************************************
// Global Variables
//****************************************************
Viewport viewport;
parsedScene Scene;
vector<vector<Color> > imageBuffer;

double cosval = 0;
bool flagDrawToScreen = false;
bool flagDrawToFile = false;

//****************************************************
// reshape viewport if the window is resized
//****************************************************
void myReshape(int w, int h) {
  viewport.w = w;
  viewport.h = h;

  glViewport (0,0,viewport.w,viewport.h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0, viewport.w, 0, viewport.h);

}

void myKeyboardFunc(unsigned char key, int x, int y){
	cout << " ---exiting program" << "\n";

    exit(0);
}

Color trace(Ray ray, int depth, Color baseColor){
    if (depth > Scene.reflectiondepth){
        return Color(0, 0, 0);
    }


    Point inter = Point(0, 0, 0);  //ray from camera intersects a shape at this point
    Shape* intershape = NULL; //which shape is intersected. note you can't instantiate abstract class
    if(!Scene.shapes.checkIntersect(ray, &inter, intershape, LARGE_NUM)) //camera ray hits nothing, return background
        return Color(0, 0, 0);   
          
    BRDF brdf = intershape->brdf; 
    Vector N = intershape->getNormal(inter); 
	size_t lightCount = Scene.lights.size();
    for(size_t i = 0; i < lightCount; i++){
        Ray lightray, L;
        if(Scene.lights[i].directional) {
            lightray = Ray(inter, Vector(Point(0,0,0), Scene.lights[i].source).normalize()); 
            
        } else {
            lightray = Ray(inter, Vector(inter, Scene.lights[i].source).normalize());
             
		}
		
		double dist = sqrt(pow(Scene.lights[i].source.x - inter.x, 2) + 
						pow(Scene.lights[i].source.y - inter.y, 2) + 
						pow(Scene.lights[i].source.z - inter.z, 2));
		
		
        Point bias = inter + lightray.direction * 0.001; //take into account shadow bias
        L = Ray(bias, lightray.direction); 
        if(!Scene.shapes.checkIntersect(L, dist)){ //light ray for this light is not blocked by any shapes 
            baseColor += brdf.kd * Scene.lights[i].color * max(0.0, L.direction.dotProduct(N)); 

			Vector R = getReflection(L.direction, N);
            baseColor += brdf.ks * Scene.lights[i].color * pow(max(0.0, R.dotProduct(Vector(inter, ray.origin).normalize())), brdf.sp); 
        } 
    }
    
    
    baseColor += brdf.ke; // the emission of this object

    baseColor += Scene.ambient; // the overal ambient glow of the scene.
    

    //if (brdf.kr > 0){
    if (brdf.ks.r > 0 || brdf.ks.g > 0 || brdf.ks.b > 0){
		Vector reflectDir = getReflection(ray.direction, N).negative();
		Ray reflectedRay = Ray(inter + reflectDir * 0.1, reflectDir); //bias
		baseColor += trace(reflectedRay, depth + 1, baseColor) * brdf.ks; 
    }
    
    return baseColor; 
    
}



void drawScreen() {
    
    //calculations for the ray from the camera to the screen
    Vector look_vector = Vector(Scene.lookfrom, Scene.lookat).normalize();
    Vector up_dir = Scene.up_dir;
    
    
    Vector right_dir = up_dir.crossProduct(look_vector);
    right_dir = right_dir.normalize();

    up_dir = right_dir.crossProduct(look_vector);
    
    double fov = Scene.fov * PI / 180.0;
    double rat = (double(Scene.width)/double(Scene.height));
    double iph = tan(fov/2);
    double ipw = tan(fov/2);
   
    Vector uv = up_dir*iph;
    Vector rv = right_dir*ipw*rat;
    
    Point imgc = Scene.lookfrom + look_vector;
    Point UL = imgc +  uv     + (rv*-1);
    Point UR = imgc +  uv     +  rv;
    Point LL = imgc + (uv*-1) + (rv*-1);
    Point LR = imgc + (uv*-1) +  rv;
    
    Point point, point1, point2; 
    Ray ray; 
    
    Color allColors = Color(0, 0, 0);
    
    int x, y;
    double u, v; 
	#pragma omp parallel for private(x, y, u, v, point, point1, point2, ray, Scene, fov, rat, iph, ipw, uv, rv, imgc, UL, UR, LL, LR)
        for (x = 0; x < Scene.width ; x += 1) {
            for (y = 0; y < Scene.height; y += 1) {
                u = double(x)/Scene.width;
                v = double(y)/Scene.height;
                
                point1 = (LL*v) + (UL*(1-v));
                point1 = point1 * u;
                
                point2 = (LR*v) + (UR*(1-v));
                point2 = point2 * (1-u);
                point = point1 + point2;
                
                ray = Ray(Scene.lookfrom, point);
                ray.direction = ray.direction.normalize();

                imageBuffer[x][y] = trace(ray, 0, Color(0, 0, 0));
                
        //some printing to keep track of progress                 
        //if((x%int(Scene.width/8)) == 0 && y == 0)  cout << "." << endl; 
        //if(x == int(Scene.width/2) && y == int(Scene.height/2)) cout << "  halfway!" << endl; 
			}//for, y
        }//for, x
        
   
}






//****************************************************
// function for openGl screen rendering.
//***************************************************
bool drawn = false;
void myDisplay() {
	glClear(GL_COLOR_BUFFER_BIT);				// clear the color buffer

	glMatrixMode(GL_MODELVIEW);			        // indicate we are specifying camera transformations
	glLoadIdentity();				        // make sure transformation is "zero'd"

	glBegin(GL_POINTS); 

	Color c;
	for (int x = 0; x < viewport.w; x++) {		
		for (int y = 0; y < viewport.h; y++) {
			c = imageBuffer[x][y];
			//cout << "c.r" << c.r << endl;
			glColor3f(c.r, c.g, c.b);
			glVertex2f(x + 0.5, y + 0.5);  
		}
	}
	
	glEnd(); 

	glFlush();
	glutSwapBuffers();					// swap buffers (we earlier set double buffer)
	
	
}

int main(int argc, char *argv[]) {
    try {
		
		if (argc > 1) {
			std::string fname = std::string(argv[1]);
			cout << " ---Loading from file '";
			cout << fname << "'" << endl;
			Scene = loadScene(fname);
			
			if (argc > 2) { 
				int c = 2;
				while (c < argc) {
					if (!std::string(argv[c]).compare("-screen")) {
						flagDrawToScreen = true;
					} else if (!std::string(argv[c]).compare("-file")) {
						flagDrawToFile = true;
					}
					c++;
				}
			} else {
				throw 6; // no destination
			}
		} else {
			throw 5; // not input file
		}
		
		if (flagDrawToScreen)
			cout << " ---Writing to screen" << endl;
		
		imageBuffer = vector<vector<Color> >(Scene.width, vector<Color>(Scene.height, Color(0,0,0)));
	    
		double time = 0.0;
		time = timestamp();
		drawScreen(); // the main computation hog
		printf("Time %f", timestamp() - time);
		
		if (flagDrawToScreen){
			cout << endl << endl << "::: To Screen :::" << endl << endl;
			glutInit(&argc, argv);

			//This tells glut to use a double-buffered window with red, green, and blue channels 
			glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);

			// Initalize theviewport size
			viewport.w = Scene.width;
			viewport.h = Scene.height;

			//The size and position of the window
			glutInitWindowSize(viewport.w, viewport.h);
			glutInitWindowPosition(0,0);
			glutCreateWindow(argv[0]);

			glutDisplayFunc(myDisplay);				// function to run when its time to draw something
			glutReshapeFunc(myReshape);				// function to run when the window gets resized
			glutKeyboardFunc(myKeyboardFunc);

			glutMainLoop();							// infinite loop that will keep drawing and resizing
			
		}
		
		

		cout << " ---Ray Tracer finished" << endl;

	} catch (int e) {
		if (e == 8) {
			cout << "\n";
			cout << "Unable to open file: '" << std::string(argv[1]) << "'" << endl;
			cout << "\n";
		} else if (e == 6) {
			cout << "\n";
			cout << "Please specify destination(s) of image with these flags:" << endl << endl;
			cout << "-screen" << endl;
			cout << "-file" << endl;
			cout << "\n";
		} else {	
			cout << "\n";
			cout << "Mising input.  Please use ./scene ___.test" << "\n";
			cout << "\n";
		}
	}

    return 0;
}






