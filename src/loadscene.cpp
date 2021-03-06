#include "loadscene.h"
#include "primitives.h"

#include <cstdio>
using namespace std;


parsedScene loadScene(std::string file) {

  //store variables and set stuff at the end
  parsedScene S;  
  S.outputFileName = "output.png";  
  
    float minX = -1.0;
    float maxX = 1.0;
    
    float minY = -1.0;
    float maxY = 1.0;
    
    float minZ = -1.0;
    float maxZ = 1.0;

  std::ifstream inpfile(file.c_str());
  if(!inpfile.is_open()) {
    //std::cout << "Unable to open file" << std::endl;
	throw 8;
  } else {
    std::string line;

	int lightIdx = 0;
	int triangleIdx = 0;
	int vertexIdx = 0;
	int vertexCount = 0;

    Vector3 (*vertices);
	Brdf brdf_t = Brdf();
    
    Color tkd, tks, tke;
    setVector3(tkd, 0.2, 0.2, 0.2);
    setVector3(tks, 0, 0, 0);
    setVector3(tke, 0, 0, 0);
    
    brdf_t.kd = tkd;
    brdf_t.ks = tks;
    brdf_t.ke = tke;
    brdf_t.kr = 0;
    brdf_t.sp = 0;
    
    /*{0.2, 0.2, 0.2, 0, //kd
				 0, 0, 0, 0,//ks
				 0, 0, 0, 0,//ke
				 0, 0, 0, 0,//kr
				 1, 0, 0, 0};//sp*/

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
		  setVector3(S.lookfrom, 
					 atof(splitline[1].c_str()), //x
					 atof(splitline[2].c_str()), //y
					 atof(splitline[3].c_str())); //z
		  setVector3(S.lookat,
					 atof(splitline[4].c_str()), //x
					 atof(splitline[5].c_str()), //y
					 atof(splitline[6].c_str())); //z
		  setVector3(S.up_dir,
					 atof(splitline[7].c_str()), //x
					 atof(splitline[8].c_str()), //y
					 atof(splitline[9].c_str())); //z
		  S.fov = atof(splitline[10].c_str());		       

      }

      //sphere x y z radius
      //  Deﬁnes a sphere with a given position and radius.
      else if(!splitline[0].compare("sphere")) {
		    //dont care
      }
      //maxverts number
      //  Deﬁnes a maximum number of vertices for later triangle speciﬁcations. 
      //  It must be set before vertices are deﬁned.
      else if(!splitline[0].compare("maxverts")) {
		  vertexCount = atoi(splitline[1].c_str());
		  vertices = new Vector3[vertexCount];
      }
	  //maxtri number
      //  Deﬁnes a maximum number of triangles for later triangle speciﬁcations. 
      //  It must be set before vertices are deﬁned.
      else if(!splitline[0].compare("maxtri")) {
		  S.triangleCount = atoi(splitline[1].c_str());
		  S.triangleList = new __m128[S.triangleCount][9];
          S.bbList = new BoundingBox[S.triangleCount]; 
		  //S.brdfList = new Brdf[S.triangleCount];
      }
	  //maxlight number
      //  Deﬁnes a maximum number of lights for later triangle speciﬁcations. 
      //  It must be set before vertices are deﬁned.
      else if(!splitline[0].compare("maxlight")) {
		  S.lightCount = atoi(splitline[1].c_str());
		  S.lightList = new __m128[S.lightCount][3];
      }
      //maxvertnorms number
      //  Deﬁnes a maximum number of vertices with normals for later speciﬁcations.
      //  It must be set before vertices with normals are deﬁned.
      else if(!splitline[0].compare("maxvertnorms")) {
			//dont care
      }
      //vertex x y z
      //  Deﬁnes a vertex at the given location.
      //  The vertex is put into a pile, starting to be numbered at 0.
      else if(!splitline[0].compare("vertex")) {
        // Create a new vertex with these 3 values, store in some array
		  if(vertexIdx < vertexCount){
              float u = atof(splitline[1].c_str());
              float v = atof(splitline[2].c_str());
              float w = atof(splitline[3].c_str());
              
              if (u < minX) { minX = u; }
              if (u > maxX) { maxX = u; }
              
              if (v < minY) { minY = v; }
              if (v > maxY) { maxY = v; }
              
              if (w < minZ) { minZ = w; }
                if (w > maxZ) { maxZ = w; }
              
			setVector3( vertices[vertexIdx],
						u, //x
						v, //y
						w); //z
			++vertexIdx;
		  }//if
		  else{
			  cout << "more verteices than specified" << endl;
			  exit(1);
		  }//else
      }
      //vertexnormal x y z nx ny nz
      //  Similar to the above, but deﬁne a surface normal with each vertex.
      //  The vertex and vertexnormal set of vertices are completely independent
      //  (as are maxverts and maxvertnorms).
      else if(!splitline[0].compare("vertexnormal")) {
			//dont care
      }
      //tri v1 v2 v3
      //  Create a triangle out of the vertices involved (which have previously been speciﬁed with
      //  the vertex command). The vertices are assumed to be speciﬁed in counter-clockwise order. Your code
      //  should internally compute a face normal for this triangle.
      else if(!splitline[0].compare("tri")) {
		  if(triangleIdx < S.triangleCount){
            constVector3 p1 = vertices[atoi(splitline[1].c_str())];
            constVector3 p2 = vertices[atoi(splitline[2].c_str())];
            constVector3 p3 = vertices[atoi(splitline[3].c_str())];
            
			setTriangle( S.triangleList[triangleIdx],
						 p1, //p1
						 p2, //p2
						 p3, brdf_t); //p3
			//brdfCopy(S.brdfList[triangleIdx], brdf_t);
            S.bbList[triangleIdx] = BoundingBox(S.triangleList[triangleIdx]); 
			++triangleIdx;
		  }//if
		  else{
			  cout << "more triangles than specified" << endl;
			  exit(1);
		  }//else		  
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
			//dont care
      }

      //translate x y z
      //  A translation 3-vector
      else if(!splitline[0].compare("translate")) {
          //dont care
      }
      //rotate x y z angle
      //  Rotate by angle (in degrees) about the given axis as in OpenGL.
      else if(!splitline[0].compare("rotate")) {
		  //dont care
      }
      //scale x y z
      //  Scale by the corresponding amount in each axis (a non-uniform scaling).
      else if(!splitline[0].compare("scale")) {
		  //dont care
      }
      //pushTransform
      //  Push the current modeling transform on the stack as in OpenGL. 
      //  You might want to do pushTransform immediately after setting 
      //   the camera to preserve the “identity” transformation.
      else if(!splitline[0].compare("pushTransform")) {
          //dont care
      }
      //popTransform
      //  Pop the current transform from the stack as in OpenGL. 
      //  The sequence of popTransform and pushTransform can be used if 
      //  desired before every primitive to reset the transformation 
      //  (assuming the initial camera transformation is on the stack as 
      //  discussed above).
      else if(!splitline[0].compare("popTransform")) {
		  //dont care
      }

      //directional x y z r g b
      //  The direction to the light source, and the color, as in OpenGL.
      else if(!splitline[0].compare("directional")) {
		  if(lightIdx < S.lightCount){
			setLight(&(S.lightList[lightIdx][0]),
					 atof(splitline[1].c_str()), //x
					 atof(splitline[2].c_str()), //y
					 atof(splitline[3].c_str()), //z
					 atof(splitline[4].c_str()), //r
					 atof(splitline[5].c_str()), //g
					 atof(splitline[6].c_str()), //b
					 true); //isDir
			++lightIdx;
		  }//if
		  else{
			  cout << "more lights than specified" << endl;
			  exit(1);
		  }//else  
        // add light to scene...
      }
      //point x y z r g b
      //  The location of a point source and the color, as in OpenGL.
      else if(!splitline[0].compare("point")) {
		  if(lightIdx < S.lightCount){
			setLight(&(S.lightList[lightIdx][0]),
					 atof(splitline[1].c_str()), //x
					 atof(splitline[2].c_str()), //y
					 atof(splitline[3].c_str()), //z
					 atof(splitline[4].c_str()), //r
					 atof(splitline[5].c_str()), //g
					 atof(splitline[6].c_str()), //b
					 false); //isDir
			++lightIdx;
		  }//if
		  else{
			  cout << "more lights than specified" << endl;
			  exit(1);
		  }//else  
        // add light to scene...
      }
      //attenuation const linear quadratic
      //  Sets the constant, linear and quadratic attenuations 
      //  (default 1,0,0) as in OpenGL.
      else if(!splitline[0].compare("attenuation")) {
			//dont care
      }
      //ambient r g b
      //  The global ambient color to be added for each object 
      //  (default is .2,.2,.2)
      else if(!splitline[0].compare("ambient")) {
		setVector3( S.ambient,
					atof(splitline[1].c_str()),  //r
					atof(splitline[2].c_str()),  //g
					atof(splitline[3].c_str())); //b
        /*The global ambient color to be added for each object (default is .2,.2,.2)*/
      }

      //diffuse r g b
      //  speciﬁes the diﬀuse color of the surface.
      else if(!splitline[0].compare("diffuse")) {
          Color kd;
		setVector3( kd,
					atof(splitline[1].c_str()),  //r
					atof(splitline[2].c_str()),  //g
					atof(splitline[3].c_str())); //b     
        brdf_t.kd = kd;
        // Update current properties
      }
      //specular r g b 
      //  speciﬁes the specular color of the surface.
      else if(!splitline[0].compare("specular")) {
          Color ks;
		  setVector3( ks,
					atof(splitline[1].c_str()),  //r
					atof(splitline[2].c_str()),  //g
					atof(splitline[3].c_str())); //b   
            brdf_t.ks = ks;
        // Update current properties
      }
      //shininess s
      //  speciﬁes the shininess of the surface.
      else if(!splitline[0].compare("shininess")) {
          brdf_t.sp = atof(splitline[1].c_str());
          //setVector3(brdf_t[BRDF_SP_IDX], atof(splitline[1].c_str()), 0, 0);
      }
      //emission r g b
      //  gives the emissive color of the surface.
      else if(!splitline[0].compare("emission")) {
          Color ke;
        setVector3( ke,
					atof(splitline[1].c_str()),  //r
					atof(splitline[2].c_str()),  //g
					atof(splitline[3].c_str())); //b       
        brdf_t.ke = ke;
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
    S.root = AABB_Node(S.triangleList, 0, S.triangleCount); 
    
	printf("minx: %f maxx: %f miny: %f maxy: %f minz: %f maxz %f\n", minX, maxX, minY, maxY, minZ, maxZ);
    
    
	//delete obj
	delete[] vertices;

    //cout << "light x" << S.lights[0].x << ", " <<  S.lights[1].x << endl;
    
    return S;
  }
  throw 5;

} // load scene

parsedScene::parsedScene(){
    width = 400;
    height = 300;
    reflectiondepth = 5;
	triangleCount = 0;
	lightCount = 0;
    
    aaJitter = 0.0; // 1.0
    aaDepth = 1.0; // 4.0
    // aa off by default
    
	setVector3(lookfrom, 0, 0, 1);
	setVector3(lookat, 0, 0, 0);
	setVector3(up_dir, 0, 1, 0);
    fov = 33; // in degrees
    
	setVector3(ambient, 0.2, 0.2, 0.2);
    
    triangleList = NULL;
	lightList = NULL;    
} //parsedScene()


