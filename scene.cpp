#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <cstdlib>

#include <nmmintrin.h>
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

#include "src/primitives.h"
#include "src/loadscene.h"
#include "src/imgwriter/lodepng.h"

float timestamp()
{
	struct timeval tv;
	gettimeofday(&tv, 0);
	return 1e-6*tv.tv_usec;
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
float (*imageBuffer)[4];

float cosval = 0;
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

void trace(Ray ray, const int depth, Color &baseColor){
    if (depth > Scene.reflectiondepth){
        return;
    }
    
	//Vector3 inter;  //ray from camera intersects a shape at this point
    Vector3 tempV1, tempV2;
	Color tempC1;
	Color reflectColor; // color from reflected ray
    setVector3(reflectColor, 0, 0, 0);
	Ray lightray, reflectedRay;

    float large = LARGE_NUM; 
    Triangle* intertri = new Triangle[9]; 

    float t = Scene.root.CollisionTest(ray, intertri, &large); 
    if(t == 0 || t == NO_INTERSECTION) return; 

    //calculate point of intersection using t value 
    Vector3 intersect; 
    vector3Scale(intersect, ray[RAY_DIRECTION_IDX], t); 
    vector3Add(intersect, intersect, ray[RAY_ORIGIN_IDX]); 
    
    __m128 N = (*intertri)[TRIANGLE_NORMAL_IDX];
   
	float temp1;
	float dist = LARGE_NUM;
	__m128 *light = NULL;
    
	//color from lights
    for(int i = 0; i < Scene.lightCount; ++i){
		light = &(Scene.lightList[i][0]);

        float isDir[4];
        _mm_storeu_ps(isDir, light[LIGHT_ISDIR_IDX]);
        if(isDir[0] != 0) {
			setRayByVector(lightray, intersect, light[LIGHT_SRC_IDX]);
        } // if
		else { 
			vector3Sub(tempV1, light[LIGHT_SRC_IDX], intersect);
            setRayByVector(lightray, intersect, tempV1);
            dist = sqrt(vector3Dot(tempV1, tempV1));
        }//

		vector3Scale(tempV2, lightray[RAY_DIRECTION_IDX], 0.0001); // shadow bias
		vector3Add(lightray[RAY_ORIGIN_IDX], lightray[RAY_ORIGIN_IDX], tempV2);// shadow bias
        //if(!hasIntersect(Scene.triangleList, Scene.triangleCount, lightray, dist)){ //light ray for this light is not blocked by any shapes. 
        
        if(!Scene.root.CollisionTest(lightray, &dist)){    
			vector3Scale(tempC1, light[LIGHT_COLOR_IDX], max(0.0f, vector3Dot(lightray[RAY_DIRECTION_IDX], N)));
			colorMultiply(tempC1, tempC1, (*intertri)[BRDF_KD_IDX]);
			vector3Add(baseColor, baseColor, tempC1); 
         
            getReflection(tempV1, lightray[RAY_DIRECTION_IDX], N);
			vector3Sub(tempV2, ray[RAY_ORIGIN_IDX], intersect);
			vector3Normalize(tempV2, tempV2);
            
            float sp[4];
            _mm_store_ps(sp, (*intertri)[BRDF_SP_IDX]);
			vector3Scale(tempC1, light[LIGHT_COLOR_IDX], pow(max(0.0f, vector3Dot(tempV1, tempV2)), (int) sp[0]));
			colorMultiply(tempC1, tempC1, (*intertri)[BRDF_KS_IDX]);
			vector3Add(baseColor, baseColor, tempC1); 
        } // if
    } // for   
    
    vector3Add(baseColor, baseColor, (*intertri)[BRDF_KE_IDX]); // the emission of this object
    vector3Add(baseColor, baseColor, Scene.ambient); // the overal ambient glow of the scene.

    float ks[4];
    _mm_store_ps(ks, (*intertri)[BRDF_KS_IDX]);
    //printf("specular\t");
    //printVector((*intertri)[BRDF_KS_IDX]);
    if (ks[0] > 0 || ks[1] > 0 || ks[2] > 0){
		getReflection(tempV1, ray[RAY_DIRECTION_IDX], N);
		vector3Scale(tempV1, tempV1, -1);
		vector3Scale(tempV2, tempV1, 0.001);
		vector3Add(tempV2, tempV2, intersect);//bias

		setRayByVector(reflectedRay, tempV2, tempV1);
        trace(reflectedRay, depth + 1, reflectColor);
		colorMultiply(reflectColor, reflectColor, (*intertri)[BRDF_KS_IDX]);
		vector3Add(baseColor, baseColor, reflectColor);
    }
    // if  
} // trace

void drawScreen() {
    
    //calculations for the ray from the camera to the screen
	Vector3 look_vector, up_dir, right_dir, uv, rv, imgc, UL, UR, LL, LR, point;
	Vector3 tempV1, tempV2;	
	Ray ray; 

	vector3Sub(look_vector, Scene.lookat, Scene.lookfrom);
	vector3Normalize(look_vector, look_vector);
	vector3Copy(up_dir, Scene.up_dir);
    
	vector3Cross(right_dir, up_dir, look_vector);
	vector3Normalize(right_dir, right_dir);

	vector3Cross(up_dir, right_dir, look_vector);
    
    float fov = Scene.fov * PI / 180.0;
    float rat = (float(Scene.width)/float(Scene.height));
    float iph = tan(fov/2);
    float ipw = tan(fov/2);
   
	vector3Scale(uv, up_dir, iph);
	vector3Scale(rv, right_dir, ipw * rat);
    
	vector3Add(imgc, Scene.lookfrom, look_vector);
	
	vector3Add(UL, imgc, uv);
	vector3Sub(UL, UL, rv);
	vector3Add(UR, imgc, uv);
	vector3Add(UR, UR, rv);
	vector3Sub(LL, imgc, uv);
	vector3Sub(LL, LL, rv);
	vector3Sub(LR, imgc, uv);
	vector3Add(LR, LR, rv);	   

    int x, y;
    float u, v; 
	#pragma omp parallel for private(x, y, u, v, point, tempV1, tempV2, ray)
        for (x = 0; x < Scene.width ; x += 1) {
            for (y = 0; y < Scene.height; y += 1) {
                u = float(x)/Scene.width;
                v = float(y)/Scene.height;
                
				vector3Scale(tempV2, LL, v * u);
				vector3Scale(tempV1, UL, (1 - v) * u);
				vector3Add(point, tempV2, tempV1);
                
				vector3Scale(tempV2, LR, v * (1 - u));
				vector3Scale(tempV1, UR, (1 - v) * (1 - u));
				vector3Add(point, point, tempV1);                
				vector3Add(point, point, tempV2);
                
                setRayByPoint(ray, Scene.lookfrom, point);

				Color color;
                setVector3(color, 0, 0, 0);
				trace(ray, 0, color);
                float colorv[3];
                _mm_store_ps(colorv, color);
                
                imageBuffer[x * Scene.height + y][0] = colorv[0];
                imageBuffer[x * Scene.height + y][1] = colorv[1];
                imageBuffer[x * Scene.height + y][2] = colorv[2];
                
			}//for, y
        }//for, x       
   
}//draw screen

//****************************************************
// function for openGl screen rendering.
//***************************************************
bool drawn = false;
void myDisplay() {
	glClear(GL_COLOR_BUFFER_BIT);				// clear the color buffer

	glMatrixMode(GL_MODELVIEW);			        // indicate we are specifying camera transformations
	glLoadIdentity();				        // make sure transformation is "zero'd"
        /*float lookfrom[4]; 
 	_mm_storeu_ps(lookfrom, Scene.lookfrom); 
	lookfrom[0] = lookfrom[0] + 0.5; 
	lookfrom[2] = lookfrom[2] + 0.5; 
	Scene.lookfrom = _mm_load_ps(lookfrom); */

	drawScreen(); 

	glBegin(GL_POINTS); 

	float *c;
	for (int x = 0; x < viewport.w; x++) {		
		for (int y = 0; y < viewport.h; y++) {
			c = &(imageBuffer[x * viewport.h + y][0]);
			glColor3f(c[0], c[1], c[2]);
			glVertex2f(x + 0.5, y + 0.5);  
		}
	}
	
	glEnd(); 

	glFlush();
	glutSwapBuffers();					// swap buffers (we earlier set float buffer)	
}
void myFrameMove(){
	glutPostRedisplay(); 

}

void specialKeys(int key, int x, int y){
	float lookfrom[4]; 
 	_mm_storeu_ps(lookfrom, Scene.lookfrom);  
	Scene.lookfrom = _mm_load_ps(lookfrom);
   
    switch(key){
	case GLUT_KEY_RIGHT:lookfrom[0] = lookfrom[0] + 0.5; 
	break;

	case GLUT_KEY_LEFT:lookfrom[0] = lookfrom[0] - 0.5; 
	break;

}

Scene.lookfrom = _mm_load_ps(lookfrom);

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
		
		
        
        
        
        
        
        if (flagDrawToFile) {
                        cout << endl << endl << "::: To img.png File :::" << endl;
                        imageBuffer = new float[Scene.width * Scene.height][4];
	    
                        double time = 0.0;
                        time = timestamp();
                        drawScreen(); // the main computation hog
                        printf("Time %f\n", timestamp() - time);
                        
                        float r, g, b;
                        std::vector<unsigned char> img;
                        img.resize(Scene.width * Scene.height * 4);
                        for(unsigned y = 0; y < Scene.height; y++)
                                for(unsigned x = 0; x < Scene.width; x++) {
                                        
                                        r = imageBuffer[y*Scene.width + x][0];
                                        g = imageBuffer[y*Scene.width + x][1];
                                        b = imageBuffer[y*Scene.width + x][2];
                                        
                                        img[4 * Scene.width * y + 4 * x + 0] = min(255, int(r * 255));
                                        img[4 * Scene.width * y + 4 * x + 1] = min(255, int(g * 255));
                                        img[4 * Scene.width * y + 4 * x + 2] = min(255, int(b * 255));
                                        img[4 * Scene.width * y + 4 * x + 3] = 255;
                        }
                        
                        //unsigned error = lodepng::encode( "img.png", img, Scene.width, Scene.height);
                        std::string fileName = Scene.outputFileName; 
                        unsigned error = lodepng::encode(fileName+".png", img, Scene.width, Scene.height);
                    //if there's an error, display it
                    if(error) {
                                cout << "!!! lodepng: encoder error " << error << ": "<< lodepng_error_text(error) << endl;
                        } else {
                                cout << "    " <<  fileName << " file successfully written." << endl;
                                //cout << "    " <<  "img.png " << " file successfully written." << endl;
                        }
                        
                        delete[] imageBuffer;
                } 
    
        
        
        
        
        
        
        
        
        
        
        if (flagDrawToScreen)
			cout << " ---Writing to screen" << endl;
		
		imageBuffer = new float[Scene.width * Scene.height][4];
	    
		/*double time = 0.0;
		time = timestamp();
		drawScreen(); // the main computation hog
		printf("Time %f\n", timestamp() - time);*/
		
		if (flagDrawToScreen){
			cout << endl << endl << "::: To Screen :::" << endl << endl;
			glutInit(&argc, argv);

			//This tells glut to use a float-buffered window with red, green, and blue channels 
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
			glutIdleFunc(myFrameMove); 
			glutSpecialFunc(specialKeys);

			glutMainLoop();							// infinite loop that will keep drawing and resizing
			
		}
		
		//test
		//test
		ofstream fout("output_test.txt"); 
		if(!fout) { 
			cout << "can't write to file\n"; 
			return 1; 
		}
		for(int x = 0; x < Scene.width; x++)
			for(int y = 0; y < Scene.height; y++)
			{
				if(imageBuffer[x * Scene.height + y][0] != 0 || imageBuffer[x * Scene.height + y][1] != 0 || imageBuffer[x * Scene.height + y][2] != 0)
					fout << "[x, y] " << x << ", " << y << " color " << imageBuffer[x * Scene.height + y][0] << ", " << imageBuffer[x * Scene.height + y][1] << ", " << imageBuffer[x * Scene.height + y][2] << endl;
			}//y
		//test

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

    //delete[] Scene.triangleList;
	//delete[] Scene.brdfList;
    //	delete[] Scene.lightList;
    //delete[] imageBuffer;
    return 0;
}
