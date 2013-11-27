#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>

#include <cstdlib>

#include <time.h>
#include <math.h>

#include "loadscene.h"


using namespace std;


parsedScene loadScene(std::string file) {

  //store variables and set stuff at the end
  parsedScene S;
  S.outputFileName = "output.png";
  

  std::ifstream inpfile(file.c_str());
  if(!inpfile.is_open()) {
    //std::cout << "Unable to open file" << std::endl;
	throw 8;
  } else {
    std::string line;



    vector<Point> vertices;
    BRDF brdf(Color(0.2,0.2,0.2), Color(0,0,0), Color(0,0,0), 0, 1);
    MatrixStack stack;
    
    

    while(inpfile.good()) {
      std::vector<std::string> splitline;
      std::string buf;

      std::getline(inpfile,line);
      std::stringstream ss(line);

      while (ss >> buf) {
        splitline.push_back(buf);
      }
      //Ignore blank lines
      if(splitline.size() == 0) {
        continue;
      }

      //Ignore comments
      if(splitline[0][0] == '#') {
        continue;
      }

      //Valid commands:
      //size width height
      //  must be first command of file, controls image size
      else if(!splitline[0].compare("size")) {
        S.width = atoi(splitline[1].c_str());
        S.height = atoi(splitline[2].c_str());
      }
      //maxdepth depth
      //  max # of bounces for ray (default 5)
      else if(!splitline[0].compare("maxdepth")) {
          S.reflectiondepth = atoi(splitline[1].c_str());
      }
      //output filename
      //  output file to write image to 
      else if(!splitline[0].compare("output")) {
          S.outputFileName = splitline[1];
      }

      //camera lookfromx lookfromy lookfromz lookatx lookaty lookatz upx upy upz fov
      //  speciﬁes the camera in the standard way, as in homework 2.
      else if(!splitline[0].compare("camera")) {
		  
		  S.lookfrom = Point(atof(splitline[1].c_str()),
							 atof(splitline[2].c_str()),
							 atof(splitline[3].c_str()));
		  S.lookat = Point(atof(splitline[4].c_str()),
						   atof(splitline[5].c_str()),
					       atof(splitline[6].c_str()));
		  S.up_dir = Vector(atof(splitline[7].c_str()),
						   atof(splitline[8].c_str()),
					       atof(splitline[9].c_str()));
		  S.fov = atof(splitline[10].c_str());		       

      }

      //sphere x y z radius
      //  Deﬁnes a sphere with a given position and radius.
      else if(!splitline[0].compare("sphere")) {
		    Point p = Point(atof(splitline[1].c_str()),
						    atof(splitline[2].c_str()),
						    atof(splitline[3].c_str()));
						    
			/*Matrix pr = stack.product.clone();
			Matrix translation;
			translation.m[0][3] = pr.m[0][3];
			translation.m[1][3] = pr.m[1][3];
			translation.m[2][3] = pr.m[2][3];
			translation.m[3][3] = pr.m[3][3];
			*/
			p = stack.productT * p;  // translate center point
			//pr.m[0][3] = 0;   // undo the transformation for the input matrix
			//pr.m[1][3] = 0;
			//pr.m[2][3] = 0;
			
		    Sphere* sp = new Sphere(/*stack.product * */p,
									atof(splitline[4].c_str()),
									stack.productR,
									stack.productR.invert(),
									stack.productS,
									stack.productS.invert());
		    (*sp).brdf = brdf.clone();
		    S.shapes.allShapes.push_back(sp); //push the pointer onto it

        // Create new sphere:
        //   Store 4 numbers
        //   Store current property values
        //   Store current top of matrix stack
        /*
         *  Matrices here????? Research
         */ 
      }
      //maxverts number
      //  Deﬁnes a maximum number of vertices for later triangle speciﬁcations. 
      //  It must be set before vertices are deﬁned.
      else if(!splitline[0].compare("maxverts")) {
        // Care if you want
        // Here, either declare array size
        // Or you can just use a STL vector, in which case you can ignore this
      }
      //maxvertnorms number
      //  Deﬁnes a maximum number of vertices with normals for later speciﬁcations.
      //  It must be set before vertices with normals are deﬁned.
      else if(!splitline[0].compare("maxvertnorms")) {
        // Care if you want
      }
      //vertex x y z
      //  Deﬁnes a vertex at the given location.
      //  The vertex is put into a pile, starting to be numbered at 0.
      else if(!splitline[0].compare("vertex")) {
        // Create a new vertex with these 3 values, store in some array
        vertices.push_back(Point(atof(splitline[1].c_str()),
								atof(splitline[2].c_str()),
								atof(splitline[3].c_str())));
      }
      //vertexnormal x y z nx ny nz
      //  Similar to the above, but deﬁne a surface normal with each vertex.
      //  The vertex and vertexnormal set of vertices are completely independent
      //  (as are maxverts and maxvertnorms).
      else if(!splitline[0].compare("vertexnormal")) {
        // x: atof(splitline[1].c_str()),
        // y: atof(splitline[2].c_str()),
        // z: atof(splitline[3].c_str()));
        // nx: atof(splitline[4].c_str()),
        // ny: atof(splitline[5].c_str()),
        // nz: atof(splitline[6].c_str()));
        // Create a new vertex+normal with these 6 values, store in some array
      }
      //tri v1 v2 v3
      //  Create a triangle out of the vertices involved (which have previously been speciﬁed with
      //  the vertex command). The vertices are assumed to be speciﬁed in counter-clockwise order. Your code
      //  should internally compute a face normal for this triangle.
      else if(!splitline[0].compare("tri")) {
		  Triangle* tr = new Triangle(stack.product * vertices[atof(splitline[1].c_str())],
									  stack.product * vertices[atof(splitline[2].c_str())],
									  stack.product * vertices[atof(splitline[3].c_str())],
									  stack.productR/* * stack.productS.invert()*/);
		  (*tr).brdf = brdf.clone();
		  S.shapes.allShapes.push_back(tr);
        // v1: atof(splitline[1].c_str())
        // v2: atof(splitline[2].c_str())
        // v3: atof(splitline[3].c_str())
        // Create new triangle:
        //   Store pointer to array of vertices
        //   Store 3 integers to index into array
        //   Store current property values
        //   Store current top of matrix stack
      }
      //trinormal v1 v2 v3
      //  Same as above but for vertices speciﬁed with normals.
      //  In this case, each vertex has an associated normal, 
      //  and when doing shading, you should interpolate the normals 
      //  for intermediate points on the triangle.
      else if(!splitline[0].compare("trinormal")) {
        // v1: atof(splitline[1].c_str())
        // v2: atof(splitline[2].c_str())
        // v3: atof(splitline[3].c_str())
        // Create new triangle:
        //   Store pointer to array of vertices (Different array than above)
        //   Store 3 integers to index into array
        //   Store current property values
        //   Store current top of matrix stack
      }

      //translate x y z
      //  A translation 3-vector
      else if(!splitline[0].compare("translate")) {
          Matrix S('t', atof(splitline[1].c_str()),
						atof(splitline[2].c_str()),
						atof(splitline[3].c_str()),
						0.0);
		  stack.addTransform(S, 't');
		  
        // Update top of matrix stack
      }
      //rotate x y z angle
      //  Rotate by angle (in degrees) about the given axis as in OpenGL.
      else if(!splitline[0].compare("rotate")) {
		  Matrix S('r', atof(splitline[1].c_str()),
						atof(splitline[2].c_str()),
						atof(splitline[3].c_str()),
						double(atof(splitline[4].c_str())) * (-3.14156/180.0));
		
		  stack.addTransform(S, 'r');

        // Update top of matrix stack
      }
      //scale x y z
      //  Scale by the corresponding amount in each axis (a non-uniform scaling).
      else if(!splitline[0].compare("scale")) {
		  Matrix S('s', atof(splitline[1].c_str()),
						atof(splitline[2].c_str()),
						atof(splitline[3].c_str()),
						0.0);
		
		  stack.addTransform(S, 's');
        // Update top of matrix stack
      }
      //pushTransform
      //  Push the current modeling transform on the stack as in OpenGL. 
      //  You might want to do pushTransform immediately after setting 
      //   the camera to preserve the “identity” transformation.
      else if(!splitline[0].compare("pushTransform")) {
          stack.push();
      }
      //popTransform
      //  Pop the current transform from the stack as in OpenGL. 
      //  The sequence of popTransform and pushTransform can be used if 
      //  desired before every primitive to reset the transformation 
      //  (assuming the initial camera transformation is on the stack as 
      //  discussed above).
      else if(!splitline[0].compare("popTransform")) {
		  stack.pop();
      }

      //directional x y z r g b
      //  The direction to the light source, and the color, as in OpenGL.
      else if(!splitline[0].compare("directional")) {
		  Light l =  Light(atof(splitline[1].c_str()),
								atof(splitline[2].c_str()),
								atof(splitline[3].c_str()),
								atof(splitline[4].c_str()),
								atof(splitline[5].c_str()),
								atof(splitline[6].c_str()),
								true);
		  S.lights.push_back(l); //push the actual object to it, since its not abstract like shapes
        // add light to scene...
      }
      //point x y z r g b
      //  The location of a point source and the color, as in OpenGL.
      else if(!splitline[0].compare("point")) {
		  Light l =  Light(atof(splitline[1].c_str()),
								atof(splitline[2].c_str()),
								atof(splitline[3].c_str()),
								atof(splitline[4].c_str()),
								atof(splitline[5].c_str()),
								atof(splitline[6].c_str()),
								false);
		  S.lights.push_back(l);
        // add light to scene...
      }
      //attenuation const linear quadratic
      //  Sets the constant, linear and quadratic attenuations 
      //  (default 1,0,0) as in OpenGL.
      else if(!splitline[0].compare("attenuation")) {
        // const: atof(splitline[1].c_str())
        // linear: atof(splitline[2].c_str())
        // quadratic: atof(splitline[3].c_str())
        //cout << "TODO: attenuation" << endl;
      }
      //ambient r g b
      //  The global ambient color to be added for each object 
      //  (default is .2,.2,.2)
      else if(!splitline[0].compare("ambient")) {
        S.ambient = Color(atof(splitline[1].c_str()),  //r
						  atof(splitline[2].c_str()),  //g
						  atof(splitline[3].c_str())); //b
        /*The global ambient color to be added for each object (default is .2,.2,.2)*/
      }

      //diffuse r g b
      //  speciﬁes the diﬀuse color of the surface.
      else if(!splitline[0].compare("diffuse")) {
        brdf.kd = Color(atof(splitline[1].c_str()),  //r
						atof(splitline[2].c_str()),  //g
						atof(splitline[3].c_str())); //b
        // Update current properties
      }
      //specular r g b 
      //  speciﬁes the specular color of the surface.
      else if(!splitline[0].compare("specular")) {
		  brdf.ks = Color(atof(splitline[1].c_str()),  //r
						atof(splitline[2].c_str()),  //g
						atof(splitline[3].c_str())); //b
      }
      //shininess s
      //  speciﬁes the shininess of the surface.
      else if(!splitline[0].compare("shininess")) {
          brdf.sp = atof(splitline[1].c_str());
      }
      //emission r g b
      //  gives the emissive color of the surface.
      else if(!splitline[0].compare("emission")) {
        brdf.ke = Color(atof(splitline[1].c_str()),  //r
						atof(splitline[2].c_str()),  //g
						atof(splitline[3].c_str())); //b
        // Update current properties
      } else if(!splitline[0].compare("aa")) {
        S.aaDepth = atof(splitline[1].c_str());
        S.aaJitter = atof(splitline[2].c_str());
        
        cout << " ---Anti-aliasing flag turned on:" << endl;
						
        // Update current properties
      } else {
        std::cerr << "Unknown command: " << splitline[0] << std::endl;
      }
    }

    inpfile.close();
    

    //cout << "light x" << S.lights[0].x << ", " <<  S.lights[1].x << endl;
    
    return S;
  }
  throw 5;

}


parsedScene::parsedScene(){
    width = 400;
    height = 300;
    reflectiondepth = 5;
    
    aaJitter = 0.0; // 1.0
    aaDepth = 1.0; // 4.0
    // aa off by default
    
    lookfrom = Point(0,0,1); 
    lookat = Point(0,0,0); 
    up_dir = Vector(0,1,0); 
    fov = 33; // in degrees
    
    ambient = Color(0.2, 0.2, 0.2);
    
    //ShapeList shapes;
    
}


