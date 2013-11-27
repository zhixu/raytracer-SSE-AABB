#include <cstdlib>


#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

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

#include "primitives.h"


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
Viewport	viewport;



//****************************************************
// Vector operations and functions
//****************************************************


//given a point on the unit sphere, return the normal vector 
Vector getNormal(Point p){
    return Vector(p.x, p.y, p.z).normalize();
}
//given normalized light and normal vectors, return normal reflection vector
Vector getReflection(Vector l, Vector n){
    l = l.normalize();
    n = n.normalize(); 
    double dp = l.dotProduct(n); 
    Vector v =  (n * 2 * dp) - l; 
    return v.normalize(); 
    
    
}
//given point on sphere and location of light, return normal direction vector to light
Vector getLight(Point p, Point l){
    Vector v = l.subtract(p); 
    return v.normalize();
    
}




//****************************************************
// Simple init function
//****************************************************

vector<Light> lights;
double cosval = 0;
bool animated = true;
// TODO add keyboard input for animation



void initScene(){
  int radius = min(viewport.w, viewport.h) / 3.0;
}


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
	cout << "Exiting program" << "\n";

    exit(0);
}

//****************************************************
// A routine to set a pixel by drawing a GL point.  This is not a
// general purpose routine as it assumes a lot of stuff specific to
// this example.
//****************************************************

void setPixel(int x, int y, Color c) {
  glColor3f(c.r, c.g, c.b);
  glVertex2f(x + 0.5, y + 0.5);   // The 0.5 is to target pixel
  // centers 
  // Note: Need to check for gap
  // bug on inst machines.
}

//****************************************************
// Draw a filled circle.  
//****************************************************




std::string termAmbientCoStr ("-ka");
Color ka = Color(0.0, 0.0, 0.0); 
std::string termDiffuseCoStr ("-kd");
Color kd = Color(0.0, 0.0, 0.0); 
std::string termSpecularCoStr ("-ks");
Color ks = Color(0.0, 0.0, 0.0); 
std::string termSpecularExpStr ("-sp");
int sp = 1;
std::string termPointLightStr ("-pl");
std::string termDirLightStr ("-dl");
std::string termAnimated ("-ani");
std::string termMS ("-ms");





void circle(float centerX, float centerY, float radius) {
    // Draw inner circle
    glBegin(GL_POINTS);	

  int i,j;  // Pixel indices

  int minI = max(0,(int)floor(centerX-radius));
  int maxI = min(viewport.w-1,(int)ceil(centerX+radius));

  int minJ = max(0,(int)floor(centerY-radius));
  int maxJ = min(viewport.h-1,(int)ceil(centerY+radius));

  for (i=0;i<viewport.w;i++) {
    for (j=0;j<viewport.h;j++) {

      // Location of the center of pixel relative to center of sphere
      float x = (i+0.5-centerX);
      float y = (j+0.5-centerY);

      float dist = sqrt(sqr(x) + sqr(y));

      if (dist<=radius) {
          // This is the front-facing Z coordinate
          float z = sqrt(radius*radius-dist*dist);
          Point point = Point(x, y, z); 
          Vector  N = getNormal(point); 
        
          Color baseColor = Color(0.0, 0.0, 0.0); 
          Color ambColor = Color(0.0, 0.0, 0.0); 
        
          for (int k = 0; k < lights.size(); k++) {
			  Light light = lights[k];

			  Vector L;
			  if (light.directional) {
                  L = Vector(-light.x*radius, -light.y*radius, -light.z*radius).normalize(); 
			  } else {
                  L = getLight(point, Point(light.x*radius, light.y*radius, light.z*radius)); 
			  }
			  Vector R = getReflection(L, N); 
              /* obj.ambient * lights.ambient */
              ambColor += light.color; 
			  /* obj.diffuse * (L*N) * light.diffuse */
              baseColor += kd * light.color * max(0.0, L.dotProduct(N)); 
              /*specular*/
              baseColor += ks * light.color * pow(max(0.0, R.dz), sp); 
          
		}
         
         baseColor += ka * ambColor; 
         setPixel(i,j, baseColor);
        }
    }
  }
  glEnd();
}
//****************************************************
// function that does the actual drawing of stuff
//***************************************************
void myDisplay() {

  glClear(GL_COLOR_BUFFER_BIT);				// clear the color buffer

  glMatrixMode(GL_MODELVIEW);			        // indicate we are specifying camera transformations
  glLoadIdentity();				        // make sure transformation is "zero'd"

  
  /*if (animated && lights.size() > 0) {
	  cosval += 0.3;
	  lights[0].initPos(2*cos(cosval), 2, 2*sin(cosval));
	  cout << "cosval:" << cosval << "\n";
  }*/
  // Start drawing

		circle(viewport.w / 2.0 , viewport.h / 2.0 , min(viewport.w, viewport.h) / 3.0);
	
  

  glFlush();
  glutSwapBuffers();					// swap buffers (we earlier set double buffer)
}
void myFrameMove() {
  //nothing here for now
#ifdef _WIN32
  Sleep(10);                                   //give ~10ms back to OS (so as not to waste the CPU)
#endif
  glutPostRedisplay(); // forces glut to call the display function (myDisplay())
}




//****************************************************
// main and parser for terminal arguments
//****************************************************


/*   Helper to parse chars out:
*   Returns the the number of extra chars it looked at
* */
int parseInput(int argc, char *argv[], int starting) {

	if (termAmbientCoStr.compare(string(argv[starting])) == 0) {
		// "-ka r g b"
		if (starting + 3 < argc) {
            ka = Color(atof(argv[starting + 1]), atof(argv[starting + 2]), atof(argv[starting + 3])); 
			return 3;
		} else {
			throw 5;
		}
	} else if (termDiffuseCoStr.compare(string(argv[starting])) == 0) {
		if (starting + 3 < argc) {
            kd = Color(atof(argv[starting + 1]), atof(argv[starting + 2]), atof(argv[starting + 3])); 
			return 3;
		} else {
			throw 5;
		}
	} else if (termSpecularCoStr.compare(string(argv[starting])) == 0) {
		if (starting + 3 < argc) {
            ks = Color(atof(argv[starting + 1]), atof(argv[starting + 2]), atof(argv[starting + 3])); 
			return 3;
		} else {
			throw 5;
		}
	} else if (termSpecularExpStr.compare(string(argv[starting])) == 0) {
		// "-sp v"
		if (starting + 1 < argc) {
			sp = atof(argv[starting + 1]);
			return 1;
		} else {
			throw 5;
		}
	} else if (termPointLightStr.compare(string(argv[starting])) == 0 ||
			termDirLightStr.compare(string(argv[starting])) == 0) {
		// "-pl x y z r g b"
		if (starting + 6 < argc) {
			bool directional = (termDirLightStr.compare(string(argv[starting])) == 0);;
            Light l = Light(atof(argv[starting + 1]),
                      atof(argv[starting + 2]),
                      atof(argv[starting + 3]),
                      atof(argv[starting + 4]),
                      atof(argv[starting + 5]),
                      atof(argv[starting + 6]),
                      directional); 
                       
			lights.push_back(l);
			return 6;
		} else {
			throw 5;
		}
	} 
	cout << "Input not supported: " << argv[starting] << "\n";
	return 0;
}



int main(int argc, char *argv[]) {





    try {
		/*  Parse input
		 */ 
		int count = 1;
		while (count < argc) {
			count += 1 + parseInput(argc, argv, count);
		}



		//This initializes glut
		glutInit(&argc, argv);

		//This tells glut to use a double-buffered window with red, green, and blue channels 
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);

		// Initalize theviewport size
		viewport.w = 400;
		viewport.h = 400;

  

		//The size and position of the window
		glutInitWindowSize(viewport.w, viewport.h);
		glutInitWindowPosition(0,0);
		glutCreateWindow(argv[0]);

		initScene();							// quick function to set up scene

		glutDisplayFunc(myDisplay);				// function to run when its time to draw something
		glutReshapeFunc(myReshape);				// function to run when the window gets resized
		//glutIdleFunc(myFrameMove);  			// this makes it animated
		glutKeyboardFunc(myKeyboardFunc);




		glutMainLoop();							// infinite loop that will keep drawing and resizing




	} catch (int e) {
		cout << "\n";
		cout << "\n";
		cout << "MISSING INPUT!" << "\n";
		cout << "\n";

	}

    return 0;
}
