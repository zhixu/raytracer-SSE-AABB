#include "primitives.h"

//****************************************************
// brdf.h. get the color of a certain pixel with given lights
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

Color getBRDF(Point point, float radius, Color kd, Color ks, Color ke, int sp, vector<Light> lights) {
    
	  Vector  N = getNormal(point); 
	
	  Color baseColor = Color(0.0, 0.0, 0.0); 
	  //Color ambColor = Color(0.0, 0.0, 0.0); 
	

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
		  //ambColor += light.color; 
		  /* obj.diffuse * (L*N) * light.diffuse */
		  baseColor += kd * light.color * max(0.0, L.dotProduct(N)); 
		  /*specular*/
		  baseColor += ks * light.color * pow(max(0.0, R.dz), sp); 
	 
		  baseColor += ke;  // * ambColor; 
		  
		  throw 7;
		  
		  //baseColor += global ambient term
		  cout << "need to addd global ambient term" << endl;
	 }
        
}
