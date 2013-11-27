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
#include "loadscene.h"


#define PI 3.14159265  // Should be used from mathlib
inline float sqr(float x) { return x*x; }

using namespace std;

//****************************************************
// Some Classes
//****************************************************

/*class Viewport;

class Viewport {
  public:
    int w, h; // width and height
};
*/


//****************************************************
// Global Variables
//****************************************************
//Viewport	viewport;



//****************************************************
// Vector operations and functions
//****************************************************

double dotProduct(vector<double> a, vector<double> b){
    if(a.size() != b.size()){
        cout << "Dot product cannot be performed!\n" << endl; 
        exit(1); 
    }
    double result = 0; 
    for(int i = 0; i < a.size() - 1; i++){  //account for last term, which is magnitude 
        result += a[i]*b[i]; 
    }
    return result; 
}

//normalize a given vector by dividing by magnitude
vector<double> normalize(vector<double> a){
    double arr[4]; 
    arr[0] = a[0] / a[3];
    arr[1] = a[1] / a[3]; 
    arr[2] = a[2] / a[3]; 
    arr[3] = 1;
    return vector<double>(arr, arr + 4); 
    
}



//given a point on the unit sphere, return the normal vector 
vector<double> getNormal(vector<double> p, int radius){
    double arr[4] = {p[0], p[1], p[2], 1}; 
    return newNormalize(vector<double>(arr, arr + 4)); 
}

//given normalized light and normal vectors, return normal reflection vector
vector<double> getReflection(vector<double> l, vector<double> n){
    double dp = dotProduct(normalize(l),normalize(n)); 
    double arr[4]; 
    arr[0] = 2 * dp * n[0] - l[0]; 
    arr[1] = 2 * dp * n[1] - l[1]; 
    arr[2] = 2 * dp * n[2] - l[2]; 
    arr[3] = 0; 
    return newNormalize(vector<double>(arr, arr + 4)); 
    
}



/*void myKeyboardFunc(unsigned char key, int x, int y){
	cout << "Exiting program" << "\n";
	
    exit(0);
}*/


//****************************************************
// main and parser for terminal arguments
//****************************************************
/*   Helper to parse chars out:
*   Returns the the number of extra chars it looked at
* */
int parseInput(int argc, char *argv[], int starting) {

	if (string(argv[starting]).compare("-o") == 0) {
		// "-ka r g b"
		cout << "YAY!" << endl;
		if (starting + 1 < argc) {
			/* 
			 *  so String at argv[starting + 1] is the path to the .obj file
			 */
			 //loadScene(std::string file);
			 cout << "obj file path is:" << argv[starting + 1] << endl;
			//ka_r = atof(argv[starting + 1]);
			//ka_g = atof(argv[starting + 2]);
			//ka_b = atof(argv[starting + 3]);
			return 1;
		} else {
			throw 5;
		}
	}
	cout << "Input not supported: " << argv[starting] << "\n";
	return 0;
}



int main(int argc, char *argv[]) {
	
	
	//testfunction();
	

    try {
		/*  Parse input
		 */ 
		int count = 1;
		while (count < argc) {
			count += 1 + parseInput(argc, argv, count);
		}



		//This initializes glut
		//glutInit(&argc, argv);

		//This tells glut to use a double-buffered window with red, green, and blue channels 
		//glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);

		// Initalize theviewport size
		//viewport.w = 400;
		//viewport.h = 400;
		
  

		//The size and position of the window
		//glutInitWindowSize(viewport.w, viewport.h);
		//glutInitWindowPosition(0,0);
		//glutCreateWindow(argv[0]);

		//initScene();							// quick function to set up scene

		//glutDisplayFunc(myDisplay);				// function to run when its time to draw something
		//glutReshapeFunc(myReshape);				// function to run when the window gets resized
		//glutIdleFunc(myFrameMove);  			// this makes it animated
		//glutKeyboardFunc(myKeyboardFunc);
		
		
		
		
		//glutMainLoop();							// infinite loop that will keep drawing and resizing

		
		

	} catch (int e) {
		cout << "\n";
		cout << "\n";
		cout << "MISSING INPUT!" << "\n";
		cout << "\n";

	}

    return 0;
}
