#include "primitives.h"
#include <algorithm>
#include <cstdio>

/* utility function */
//vector
void setVector3(Vector3 &dst, float x, float y, float z){
    float temp[4];
    temp[0] = x;
    temp[1] = y;
    temp[2] = z;
    temp[3] = 0;
    
    dst = _mm_load_ps(temp);
    
} // set vector

void vector3Add(Vector3 &dst, constVector3 v1, constVector3 v2){	
	dst = _mm_add_ps(v1, v2);
} // vector add

// v1 - v2
void vector3Sub(Vector3 &dst, constVector3 v1, constVector3 v2){
	dst = _mm_sub_ps(v1, v2);
} // vector subtract

void vector3Scale(Vector3 &dst, constVector3 v1, const float scalar){
    const __m128 sclr = _mm_set1_ps(scalar);
    dst = _mm_mul_ps(v1, sclr);
} // vector scale

float vector3Dot(constVector3 v1, constVector3 v2){
    __m128 temp = _mm_mul_ps(v1, v2);
    temp = _mm_hadd_ps(temp, temp);
    temp = _mm_hadd_ps(temp, temp);
    
    float p[4];
    _mm_storeu_ps(p, temp);
    
	return p[0]; 
} // dot product

// v1 x v2
// using http://fastcpp.blogspot.com/2011/04/vector-cross-product-using-sse-code.html
void vector3Cross(Vector3 &dst, constVector3 v1, constVector3 v2){
	Vector3 temp;
    
    temp = _mm_sub_ps(
        _mm_mul_ps(v2, _mm_shuffle_ps(v1, v1, _MM_SHUFFLE(3, 0, 2, 1))),
        _mm_mul_ps(v1, _mm_shuffle_ps(v2, v2, _MM_SHUFFLE(3, 0, 2, 1))));
    
    vector3Scale(dst, _mm_shuffle_ps(temp, temp, _MM_SHUFFLE(3, 0, 2, 1)), -1);
    
} // cross product

// using http://fastcpp.blogspot.com/2012/02/calculating-length-of-3d-vector-using.html
void vector3Normalize(Vector3 &dst, constVector3 v1){
	float p[4];

    _mm_store_ps(p, v1);
    float magnitude = sqrt(vector3Dot(v1, v1));
    p[0] = p[0]/magnitude;
    p[1] = p[1]/magnitude;
    p[2] = p[2]/magnitude;
    p[3] = 0; 
    
    dst = _mm_load_ps(p);
} // normalize

void vector3Copy(Vector3 &dst, constVector3 src){
	dst = src;
} // copy




void colorAdd(Color &dst, Color c1, Color c2) {
    dst = _mm_add_ps(c1, c2);
}

void colorScale(Color &dst, Color c1, const float scalar) {
    const __m128 sclr = _mm_set1_ps(scalar);
    dst = _mm_mul_ps(c1, sclr);
}

void colorMultiply(Color &dst, constColor c1, constColor c2){
	dst = _mm_mul_ps(c1, c2);
} // color multiply


void printVector(Vector3 v) {
    float p[4];

    _mm_store_ps(p, v);
    //if ( p[0] != p[0] || p[1] != p[1] || p[2] != p[2] || p[3] != p[3] ) {
    //if ( p[0] != 0 || p[1] != 0 || p[2] != 0 || p[3] != 0) {
        printf("VECTOR\t");
        printf("x: %f y: %f z: %f extra: %f\n", p[0], p[1], p[2], p[3]);
    //}
}

void printTriangle(Triangle t) {
    printf("--------------------------TRIANGLE------------------------------------\n");
    printf("kd\t");
    printVector(t[BRDF_KD_IDX]);
    printf("ks\t");
    printVector(t[BRDF_KS_IDX]);
    printf("ke\t");
    printVector(t[BRDF_KE_IDX]);
    printf("kr\t");
    printVector(t[BRDF_KR_IDX]);
    printf("sp\t");
    printVector(t[BRDF_SP_IDX]);
    printf("triangle points though\n");
    printVector(t[TRIANGLE_PT1_IDX]);
    printVector(t[TRIANGLE_PT2_IDX]);
    printVector(t[TRIANGLE_PT3_IDX]);
    printVector(t[TRIANGLE_NORMAL_IDX]);
}

void printRay(Ray r) {
    printf("RAY ORIGIN\t");
    printVector(r[RAY_ORIGIN_IDX]);
    printf("RAY DIRECTION\t");
    printVector(r[RAY_DIRECTION_IDX]);
}

//triangle
void setTriangle(Triangle triangle, constVector3 p1, constVector3 p2, constVector3 p3, Brdf brdf){
	Vector3 v1;
	Vector3 v2;
    //Triangle temp;
    
    Vector3 kr, sp;
    
    setVector3(kr, brdf.kr, 0, 0);
    setVector3(sp, brdf.sp, 0, 0);
	
	vector3Sub(v1, p1, p3); // v1 = p1 - p3
	vector3Sub(v2, p2, p3); // v2 = p2 - p3
	vector3Cross(triangle[TRIANGLE_NORMAL_IDX], v1, v2); // normal = v1 x v2, normalize
	vector3Normalize(triangle[TRIANGLE_NORMAL_IDX], triangle[TRIANGLE_NORMAL_IDX]);

	vector3Copy(triangle[TRIANGLE_PT1_IDX], p1);
	vector3Copy(triangle[TRIANGLE_PT2_IDX], p2);
	vector3Copy(triangle[TRIANGLE_PT3_IDX], p3);
    
    vector3Copy(triangle[BRDF_KD_IDX], brdf.kd);
    vector3Copy(triangle[BRDF_KS_IDX], brdf.ks);
    vector3Copy(triangle[BRDF_KE_IDX], brdf.ke);
    vector3Copy(triangle[BRDF_KR_IDX], kr);
    vector3Copy(triangle[BRDF_SP_IDX], sp);
} // set triangle

// check if intersect with the triangle
// http://geomalgorithms.com/a06-_intersect-2.html
bool intersectTriangle(constTriangle triangle, constRay ray, float* rr){

	Vector3 u, v, w, w0, intersect;
	float nDotV = vector3Dot(triangle[TRIANGLE_NORMAL_IDX], ray[RAY_DIRECTION_IDX]);

	if(nDotV == 0) return false; // ray parallel to triangle

	vector3Sub(u, triangle[TRIANGLE_PT2_IDX], triangle[TRIANGLE_PT1_IDX]);
	vector3Sub(v, triangle[TRIANGLE_PT3_IDX], triangle[TRIANGLE_PT1_IDX]);
	vector3Sub(w0, ray[RAY_ORIGIN_IDX], triangle[TRIANGLE_PT1_IDX]);

	float nDotW0 = -vector3Dot(triangle[TRIANGLE_NORMAL_IDX], w0);
	float r = nDotW0 / nDotV;

	if(r < 0) return false; // ray goes away from triangle

	// intersect = ray.origin + r * ray.direction
	vector3Scale(intersect, ray[RAY_DIRECTION_IDX], r);
	vector3Add(intersect, intersect, ray[RAY_ORIGIN_IDX]);

	// check intersect in triangle
	vector3Sub(w, intersect, triangle[TRIANGLE_PT1_IDX]);
	float uu = vector3Dot(u,u);
    float uv = vector3Dot(u,v);
    float vv = vector3Dot(v,v);
    float wu = vector3Dot(w,u);
    float wv = vector3Dot(w,v);
    float D = uv * uv - uu * vv;
	// get and test parametric coords
    float s = (uv * wv - vv * wu) / D;
    float t = (uv * wu - uu * wv) / D;
    if (s < 0.0 || t < 0.0 || (s + t) > 1.0) return false; // intersect is outside triangle        

	*rr = r;
    return true;                       // intersect is in triangle
} // triangle intersect test

// only check if intersect with any triangle
bool hasIntersect(__m128(*triangleList)[9] , const int triangleCount, constRay ray, float* tmax){
	float t = -1;
	bool intersect;

    for(int i = 0; i < triangleCount; i++){
       intersect = intersectTriangle(&(triangleList[i][0]), ray, &t);
       if (intersect && t < *tmax) {
           *tmax = t; 
           return true; //there is an intersection    
       }          
    } // for
    
    return false;
} // has intersect{
    



//bool nearestIntersect(__m128(*triangleList)[9], const int triangleCount, constRay ray, Vector3 &intersectPt, int &triangleIdx, float &tmax); // check every possible triangle and return the nearest
// check every possible triangle and return the nearest
float nearestIntersect(__m128(*triangleList)[9] , int triangleCount, constRay ray, Triangle* tri, float* tmax){
	float t = -1;
    float min = *tmax;
	bool intersect;

	for(int i = 0; i < triangleCount; i++){
        intersect = intersectTriangle(triangleList[i], ray, &t);
        if (intersect && t < min){
            min = t;
            for (int j = 0; j < 9; j++) {
                (*tri)[j] = triangleList[i][j]; 
            }
        } // if
    } // for

    if(min == LARGE_NUM || min > *tmax) return NO_INTERSECTION; //no intersection
    
    *tmax = min; 
    return min; 
} // nearest intersect

//ray
void setRayByValue(Ray &ray, float x, float y, float z, float dx, float dy, float dz){
	Vector3 origin, direction;
    
    setVector3(origin, x, y, z);
    ray[RAY_ORIGIN_IDX] = origin;

    setVector3(direction, dx, dy, dz);
	ray[RAY_DIRECTION_IDX] = direction;
	vector3Normalize(ray[RAY_DIRECTION_IDX], ray[RAY_DIRECTION_IDX]);
} // set ray by value

void setRayByPoint(Ray &ray, constVector3 origin, constVector3 p){
	Vector3 dir;
	vector3Sub(dir, p, origin);
	vector3Copy(ray[RAY_ORIGIN_IDX], origin);
	vector3Normalize(ray[RAY_DIRECTION_IDX], dir);
} // set ray by two point

void setRayByVector(Ray &ray, constVector3 origin, constVector3 direction){
	vector3Copy(ray[RAY_ORIGIN_IDX], origin);
	vector3Normalize(ray[RAY_DIRECTION_IDX], direction);
} // set ray by origin and direction

void getReflection(Vector3 &reflectionDirection, constVector3 rayDirection, constVector3 normal){
	float dp = vector3Dot(rayDirection, normal); 

	// reflection = (normal * 2dp) - rayDirection
	vector3Scale(reflectionDirection, normal, 2 * dp);
	vector3Sub(reflectionDirection, reflectionDirection, rayDirection);
	vector3Normalize(reflectionDirection, reflectionDirection);
    
    //printf("reflection direction\t");
    //printVector(reflectionDirection);
    
} // get reflection ray direction

//light
void setLight(Light light, float x, float y, float z, float r, float g, float b, bool isDir){
    Vector3 src, color, isDirV;
    float f = isDir ? 1 : 0;
    
    setVector3(src, x, y, z);
    setVector3(color, r, g, b);
    setVector3(isDirV, f, 0, 0);
    
	light[LIGHT_SRC_IDX] = src;
	light[LIGHT_COLOR_IDX] = color;
	light[LIGHT_ISDIR_IDX] = isDirV;
} //set light

//brdf
Brdf::Brdf() { }

Brdf::Brdf(Color tkd, Color tks, Color tke, float tkr, float tsp) {
    kd = tkd;
    ks = tks;
    ke = tke;
    kr = tkr;
    sp = tsp;
}














//*******************************************************************
//aabb tree and bounding box methods
//******************************************************************
BoundingBox::BoundingBox(){
      min_x = 0, max_x = 0, min_y = 0, max_y = 0, min_z = 0, max_z = 0; 
}

BoundingBox::BoundingBox(float minx, float maxx, float miny, float maxy, float minz, float maxz){
    min_x = minx, max_x = maxx, min_y = miny, max_y = maxy, min_z = minz, max_z = maxz; 
}

BoundingBox::BoundingBox(constTriangle tri){
    float p1f[9], p2f[9], p3f[9]; 
    __m128 p1 = tri[0];
    __m128 p2 = tri[1];
    __m128 p3 = tri[2];
    
    _mm_storeu_ps(p1f, p1);
    _mm_storeu_ps(p2f, p2);
    _mm_storeu_ps(p3f, p3);
    
    
    min_x = min(p3f[0], min(p1f[0], p2f[0])); 
    min_y = min(p3f[1], min(p1f[1], p2f[1])); 
    min_z = min(p3f[2], min(p1f[2], p2f[2])); 
    
    max_x = max(p3f[0], max(p1f[0], p2f[0])); 
    max_y = max(p3f[1], max(p1f[1], p2f[1])); 
    max_z = max(p3f[2], max(p1f[2], p2f[2])); 
}

int BoundingBox::getLongestAxis(){
    float xd = abs(max_x - min_x); 
    float yd = abs(max_y - min_y);
    float zd = abs(max_z - min_z);
    float m = max(xd, max(yd, zd));
    if(m == xd)
        return X_AXIS;
    else if (m == yd)
        return Y_AXIS;
    return Z_AXIS; 
    
}

float BoundingBox::getMidPoint(int axis){
    if(axis == X_AXIS){
        
        return (max_x + min_x)/2; 
    }
    if(axis == Y_AXIS){
        
        return (max_y + min_y)/2;  
    }
    else{
        
         return (max_z + min_z)/2; 
    }
    
}

bool BoundingBox::intersect(constRay r){ 
    
    if (!this) return false;  //how are you getting a null BB???
    float rDir[4], rOri[4]; 
    _mm_storeu_ps(rDir, r[RAY_DIRECTION_IDX]);
    _mm_storeu_ps(rOri, r[RAY_ORIGIN_IDX]); 
    
    float xmin = 0, xmax = 0, ymin = 0, ymax = 0, zmin = 0, zmax = 0;
    float xd = rDir[0]; 
    float yd = rDir[1]; 
    float zd = rDir[2]; 
    float ox = rOri[0]; 
    float oy = rOri[1]; 
    float oz = rOri[2]; 
    if (xd == 0) xd = 0.00001;
    if (yd == 0) yd = 0.00001;
    if (zd == 0) zd = 0.00001;
    
    float a_x = 1 / xd;
    if(a_x >= 0){
        xmin = a_x * (min_x - ox);
        xmax = a_x * (max_x - ox);  
    }
    else{
        xmin = a_x * (max_x - ox);
        xmax = a_x * (min_x - ox); 
        
    }
    
    float a_y = 1 / yd;
    if(a_y >= 0){
        ymin = a_y * (min_y - oy);
        ymax = a_y * (max_y - oy); 
    }
    else{
        ymin = a_y * (max_y - oy);
        ymax = a_y * (min_y - oy); 
        
    }
   
    if((xmin > ymax) || (ymin > xmax))
        return false;
    if(ymin > xmin)
        xmin = ymin;
    if(ymax < xmax)
        xmax = ymax;
      
    float a_z = 1 / zd;
    if(a_z >= 0){
        zmin = a_z * (min_z - oz);
        zmax = a_z * (max_z - oz); 
    }
    else{
        zmin = a_z * (max_z - oz);
        zmax = a_z * (min_z - oz); 
    }
    if((xmin > zmax) || (zmin > xmax))
        return false;
    if(zmin > xmin)
        xmin = zmin;
    if(zmax < xmax)
        xmax = zmax; 
    return (xmax > 0); 
   
   
}


BoundingBox getTriBB(constTriangle tri){
    return BoundingBox(tri); 
}

//given a list of triangles, get the BB enclosing all the triangles
BoundingBox getBB(__m128(*tris)[9], int size){
    float minx, miny, minz, maxx, maxy, maxz;  
    minx  = miny = minz = LARGE_NUM; 
    maxx = maxy = maxz = -LARGE_NUM; 
    BoundingBox bb = BoundingBox(); 
    
    for(int i = 0; i < size; i++){
        bb = BoundingBox(tris[i]); 
        minx = fmin(bb.min_x, minx);
        miny = fmin(bb.min_y, miny);
        minz = fmin(bb.min_z, minz);
        
        maxx = fmax(bb.max_x, maxx);
        maxy = fmax(bb.max_y, maxy); 

        maxz = fmax(bb.max_z, maxz);
    }
    
    return BoundingBox(minx, maxx, miny, maxy, minz, maxz); 
    
}

//make tree by passing in allshapes to get root
//adapted from http://www.flipcode.com/archives/Dirtypunks_Column-Issue_05_AABB_Trees_Back_To_Playing_With_Blocks.shtml
AABB_Node::AABB_Node(__m128 (*triList)[9], int depth, int size){ 

    bb = getBB(triList, size); 
     if(size >  3 && depth < 2){ 
        int axis = bb.getLongestAxis(); 
        //vector<__m128> left = vector<__m128>(); 
        //vector<__m128> right = vector<__m128>(); 
        __m128 (*left)[9] = new __m128[size][9];
        __m128 (*right)[9] =  new __m128[size][9]; 
        int lc = 0; 
        int rc = 0;
     
        float mid = bb.getMidPoint(axis); 
        //divide this bb by midpoint along longest axis, sort objects
        for(int i = 0; i < size; i++){
            /*if(tris->at(i).getBB().getMidPoint(axis) < mid)
                left.push_back(tris->at(i));
            else
                right.push_back(tris->at(i)); */
            if(getTriBB(triList[i]).getMidPoint(axis) < mid){
                for (int j = 0; j < 9; j++) {
                    left[lc][j] = triList[i][j]; 
                }
                lc++; 
            }
            else{
                for (int j = 0; j < 9; j++) {
                    right[rc][j] = triList[i][j]; 
                }
                rc++; 
            }
        }
        if(lc > 0){
           if(lc == size) children[0] = new AABB_Node(left, depth + 1, lc); 
           else children[0] = new AABB_Node(left, depth, lc);   
        } 
        else children[0] = 0; 
        if(rc > 0){
            if(rc == size) children[1] = new AABB_Node(right, depth + 1, rc); 
            else children[1] = new AABB_Node(right, depth, rc); 
        } 
        else children[1] = 0; 
        
    }
    else{ //make leaf node 
        children[0] = 0; 
        children[1] = 0; 
    }
    containedTriangles = triList; 
    triCount = size; 
    
}
//simpler test, don't need to know where ray intersects with shapes
bool AABB_Node::CollisionTest(constRay ray, float* limit){
     if(!bb.intersect(ray)){
        return false; 
    }

       if(children[0] || children[1]){
            if(children[0] && children[0]->CollisionTest(ray, limit))
                return true;
            if(children[1] && children[1]->CollisionTest(ray, limit))
                return true; 
            return false; 
        }
    
        
        //do intersection with shapes at this node 
    //return checkTriIntersect(&containedTriangles, ray, limit); 
    return hasIntersect(containedTriangles, triCount, ray, limit); 

}

//return t value for ray, calculate point of intersection with t value
// if returns 0, no intersection
float AABB_Node::CollisionTest(constRay ray, Triangle* tri, float* limit){
    if(!bb.intersect(ray)){
        return NO_INTERSECTION; 
    }
    if(!children[0] && children[1]) return children[1]->CollisionTest(ray, tri, limit);
    if(!children[1] && children[0]) return children[0]->CollisionTest(ray, tri, limit); 
    if(children[0] && children[1])  return min(children[0]->CollisionTest(ray, tri, limit), children[1]->CollisionTest(ray, tri, limit)); 
        
   // return checkTriIntersect(&containedTriangles, ray, tri, limit); 
   return nearestIntersect(containedTriangles, triCount, ray, tri, limit); 
}
















