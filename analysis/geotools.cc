#include "geotools.h"
#include <iostream>
#include <cmath>

namespace geotools {

  int RayCylinderIntersection(float vertex[], float dir[], const GeoDef& geom, Intersection& SKinter) {

    // Intersection of Ray and Cylindrical Detector (adapted from use by SK)
    // Gives both the forward intersection point and the backward intersection point (if exist)
    // Assumes the SK Detector Geometry: Cylinder with axis along Z and origin at (0,0)
    // (1) Determine boundaries: wallradius and capabsz
    // (2) Determine points of intersection
    // Output
    //  intersections gives intersection points on SK boundary
    //  dist gives distance between intersection and vertex.  
    //    the sign of dist indiciates if the point is in the direction of the ray (+), or backwards (-)
    //  entries in intersections[] are ordered by forward, then dist

    SKinter.reset();

    // (0) input check
    if (dir[0]==0. && dir[1]==0. && dir[2]==0.) {
      SKinter.intersectionType = kNoIntersect;
      return (int)SKEP_NOINTERSECT;
    }

    // (1) Determine Boundaries
    float capabsz = 0.5*fabs(geom.height);
    float wallradius = geom.radius;

    if (geotool_verbose)
      std::cout << "SKENTRYPOINT: vertex=(" << vertex[0] << ", " << vertex[1] << ", " << vertex[2] << ")"
		<< " dir=(" << dir[0] << ", " << dir[1] << ", " << dir[2] << ")" << std::endl;

    // (2) Get Intersections on cylinder
    float intersections[2][3];
    int status = IntersectionLineWithCylinder( vertex, dir, wallradius, capabsz, intersections);
    if (status!=SKEP_OK) {
      SKinter.intersectionType = kNoIntersect;
      return SKEP_NOINTERSECT;
    }

    // (3) Determine distance to each intersection
    float dist[2] = {0., 0.};
    int numforwards = 0;
    int numbackwards = 0;
    for (int pt=0; pt<2; pt++) {
      dist[pt] = 0.0;
      for (int v=0; v<3; v++)
	dist[pt] += ( intersections[pt][v] - vertex[v] )*( intersections[pt][v] - vertex[v] );
      dist[pt] = sqrt(dist[pt]);
      
      // determine if forward or backward by looking at nonzero change between vertex and intersection point
      int index=0; 
      for (int v=0; v<3; v++) {
	double dd = intersections[pt][v] - vertex[v];
	if ( fabs(dir[v])>0.001 && fabs(dd)>1e-5) {
	  index = v;
	  break;
	}
      }
      
      if ( (intersections[pt][index] - vertex[index])/dir[index] < 0.0 ) {
	numbackwards++;
	dist[pt] = -1*dist[pt];
      }
      else {
	numforwards++;
      }
      if (geotool_verbose)
	std::cout << "  pt " << pt << ": ("  << intersections[pt][0] << ", " << intersections[pt][1] << ", " << intersections[pt][2] << ") " 
		  << "  r=" << sqrt( intersections[pt][0]*intersections[pt][0] + intersections[pt][1]*intersections[pt][1] ) 
		  << "  z=" << intersections[pt][2]
		  << dist[pt] << std::endl;
    }
    
    // (4) Arrange answer in struct
    if (numforwards==2) {
      // both intersections are in the forward direction of ray
      SKinter.intersectionType = k2Forward;
      int nearindex = 0;
      int farindex = 1;
      if (dist[0]>dist[1]) {
	farindex = 0;
	nearindex = 1;
      }
      SKinter.distToNearFor = dist[nearindex];
      SKinter.distToFarFor = dist[farindex];
      for (int v=0; v<3; v++) {
	SKinter.intersectionNearFor[v] = intersections[nearindex][v];
	SKinter.intersectionFarFor[v] = intersections[farindex][v];
      }
    }
    else if (numbackwards==2) {
      SKinter.intersectionType = k2Backward;
      int nearindex = 0;
      int farindex = 1;
      if (fabs(dist[0])>fabs(dist[1])) {
	farindex = 0;
	nearindex = 1;
      }
      SKinter.distToNearBack = fabs(dist[nearindex]);
      SKinter.distToFarBack = fabs(dist[farindex]);
      for (int v=0; v<3; v++) {
	SKinter.intersectionNearBack[v] = intersections[nearindex][v];
	SKinter.intersectionFarBack[v] = intersections[farindex][v];
      }
    }
    else if (numforwards==1 && numbackwards==1) {
      SKinter.intersectionType = kForBack;
      int forindex = 0;
      int backindex = 1;
      if (dist[1]>0) {
	backindex = 0;
	forindex = 1;
      }
      SKinter.distToNearFor = fabs(dist[forindex]);
      SKinter.distToNearBack = fabs(dist[backindex]);
      for (int v=0; v<3; v++) {
	SKinter.intersectionNearFor[v] = intersections[forindex][v];
	SKinter.intersectionNearBack[v] = intersections[backindex][v];
      }    
    }
    else {
      // TANGENT
      SKinter.intersectionType = kTangent;
      SKinter.distToNearFor = fabs(dist[0]);
      SKinter.distToNearBack = fabs(dist[0]);
      for (int v=0; v<3; v++) {
	SKinter.intersectionNearFor[v] = intersections[0][v];
      SKinter.intersectionNearBack[v] = intersections[0][v];
      }
    }//end of filling structure
    
    return SKEP_OK;

  }

//   // ===== End of SKEntryPoint( SKIntersect )  =====================

//   int SKEntryPoint(float vertex[], float dir[], int detboundary, float intersections[][3], float dist[] ) {

//   // Intersection of Ray and SK Detector
//   // Gives both the forward intersection point and the backward intersection point (if exist)
//   // Assumes the SK Detector Geometry: Cylinder with axis along Z and origin at (0,0)
//   // (1) Determine boundaries: wallradius and capabsz
//   // (2) Determine points of intersection
//   // Output
//   //  intersections gives intersection points on SK boundary
//   //  dist gives distance between intersection and vertex.  
//   //    the sign of dist indiciates if the point is in the direction of the ray (+), or backwards (-)
//   //  entries in intersections[] are ordered by forward, then dist

//   // Clear arrays
//   for (int pt=0; pt<2; pt++) {
//     for (int v=0; v<3; v++)
//       intersections[pt][v] = 0.0;
//     dist[pt] = 0;
//   }

//   // (1) Determine Boundaries
//   float capabsz = 0;
//   float wallradius = 0;
//   int status = GetBoundaries( detboundary, wallradius, capabsz );
//   if (status!=SKEP_OK)
//     return SKEP_NOINTERSECT;

//   if (odgeotool_verbose>0)
//     std::cout << "SKENTRYPOINT: vertex=(" << vertex[0] << ", " << vertex[1] << ", " << vertex[2] << ")"
// 	      << " dir=(" << dir[0] << ", " << dir[1] << ", " << dir[2] << ")" << std::endl;

//   // (2) Get Intersections on cylinder
//   status = IntersectionLineWithCylinder( vertex, dir, wallradius, capabsz, intersections);
//   if (status!=SKEP_OK)
//     return SKEP_NOINTERSECT;

//   // (3) Determine distance
//   for (int pt=0; pt<2; pt++) {
//     dist[pt] = 0.0;
//     for (int v=0; v<3; v++)
//       dist[pt] += ( intersections[pt][v] - vertex[v] )*( intersections[pt][v] - vertex[v] );
//     dist[pt] = sqrt(dist[pt]);

//     int index=0;
//     for (int v=0; v<3; v++)
//       if ( dir[v]!=0 ) {
// 	index = v;
// 	break;
//       }
//     if ( (intersections[pt][index] - vertex[index])/dir[index] < 0.0 )
//       dist[pt] = -1*dist[pt];
//   }

//   // (4) Order answer

//   float orderedpts[2][3];
//   float ordereddist[2];
//   if ( dist[0]*dist[1] > 0 ) {
//     // Same sign
//     if ( fabs(dist[0])<fabs(dist[1]) ) {
//       for (int v=0; v<3; v++) {
// 	orderedpts[0][v] = intersections[0][v];
// 	orderedpts[1][v] = intersections[1][v];
//       }     
//       ordereddist[0] = dist[0];
//       ordereddist[1] = dist[1];
//     }
//     else {
//       for (int v=0; v<3; v++) {
// 	orderedpts[0][v] = intersections[1][v];
// 	orderedpts[1][v] = intersections[0][v];
//       }
//       ordereddist[0] = dist[1];
//       ordereddist[1] = dist[0];
//     }
//   }
//   else {
//     // opposite sign
//     if ( dist[0]>0 ) {
//       for (int v=0; v<3; v++) {
// 	orderedpts[0][v] = intersections[0][v];
// 	orderedpts[1][v] = intersections[1][v];
//       }
//       ordereddist[0] = dist[0];
//       ordereddist[1] = dist[1];
//     }
//     else {
//       for (int v=0; v<3; v++) {
// 	orderedpts[0][v] = intersections[1][v];
// 	orderedpts[1][v] = intersections[0][v];
//       }
//       ordereddist[0] = dist[1];
//       ordereddist[1] = dist[0];
//     }
//   }

//   for (int pt=0; pt<2; pt++) {
//     for (int v=0; v<3; v++)
//       intersections[pt][v] = orderedpts[pt][v];
//     dist[pt] = ordereddist[pt];
//   }

//   return SKEP_OK;

// }




  // ----------------------------------------------------------------------------------------
  // UTILITY FUNCTIONS

  int IntersectionLineWithCircle(float m[], float b[], float o[], float r, float points[][3]) {

    // o = origin of circle
    // line defined as y = mx + b
    
    // returns -1 is no intersection
    // returns 0 if point 0 is closest to b
    // returns 1 if point 1 is closest to b

    float x1 = b[0]-o[0];
    float y1 = b[1]-o[1];
    float x2 = m[0]+b[0]-o[0];
    float y2 = m[1]+b[1]-o[1];
    
    float dx = m[0];
    float dy = m[1];
    float dr = sqrt(dx*dx+dy*dy);
    float D = x1*y2-x2*y1;
    
    float discr = r*r*dr*dr - D*D;
    if (geotool_verbose>0) 
      std::cout << "  Discriminant: " << discr << " dr " << dr << std::endl;
    if (discr<0 || (dr==0)) {
      if (geotool_verbose>0)
	std::cout << "  No Intersection!" << std::endl;
      return SKEP_NOINTERSECT;
    }
    
    float signdy = 1;
    if (dy<0)
      signdy = -1;

    //x
    points[0][0] = (D*dy + signdy*dx*sqrt(discr))/(dr*dr) + o[0]; // = X1
    points[1][0] = (D*dy - signdy*dx*sqrt(discr))/(dr*dr) + o[0]; // = X2
    //y
    points[0][1] = (-1.0*D*dx + fabs(dy)*sqrt(discr))/(dr*dr) + o[1]; // = Y1
    points[1][1] = (-1.0*D*dx - fabs(dy)*sqrt(discr))/(dr*dr) + o[1]; // = Y2
    
    float dist[2] = { 0, 0};
    for (int p=0; p<2; p++) {
      for (int v=0; v<2; v++) {
	dist[p] += (points[p][v]-b[v])*(points[p][v]-b[v]);
      }
      dist[p] = sqrt(dist[p]);
    }

    return SKEP_OK;
//     if (dist[0]<dist[1])
//       return 0;
//     else
//       return 1;

  }


// float ClosestPointonLineToPoint(float m[], float b[], float point[], int dim){
//   // dim default is 3
//   float mnorm = 0;
//   for (int v=0; v<dim; v++)
//     mnorm += m[v]*m[v];
//   mnorm = sqrt(mnorm);

//   float s = 0;
//   for (int v=0; v<dim; v++)
//     s += point[v]*m[v] - b[v]*m[v];
//   s = s/(mnorm*mnorm);

//   return s;

// }


// int SKEntryPoint(float vertex[], float dir[], int detboundary, float intersections[][3], float dist[] ) {

//   // Intersection of Ray and SK Detector
//   // Gives both the forward intersection point and the backward intersection point (if exist)
//   // Assumes the SK Detector Geometry: Cylinder with axis along Z and origin at (0,0)
//   // (1) Determine boundaries: wallradius and capabsz
//   // (2) Determine points of intersection
//   // Output
//   //  intersections gives intersection points on SK boundary
//   //  dist gives distance between intersection and vertex.  
//   //    the sign of dist indiciates if the point is in the direction of the ray (+), or backwards (-)
//   //  entries in intersections[] are ordered by forward, then dist

//   // Clear arrays
//   for (int pt=0; pt<2; pt++) {
//     for (int v=0; v<3; v++)
//       intersections[pt][v] = 0.0;
//     dist[pt] = 0;
//   }

//   // (1) Determine Boundaries
//   float capabsz = 0;
//   float wallradius = 0;
//   int status = GetBoundaries( detboundary, wallradius, capabsz );
//   if (status!=SKEP_OK)
//     return SKEP_NOINTERSECT;

//   if (odgeotool_verbose>0)
//     std::cout << "SKENTRYPOINT: vertex=(" << vertex[0] << ", " << vertex[1] << ", " << vertex[2] << ")"
// 	      << " dir=(" << dir[0] << ", " << dir[1] << ", " << dir[2] << ")" << std::endl;

//   // (2) Get Intersections on cylinder
//   status = IntersectionLineWithCylinder( vertex, dir, wallradius, capabsz, intersections);
//   if (status!=SKEP_OK)
//     return SKEP_NOINTERSECT;

//   // (3) Determine distance
//   for (int pt=0; pt<2; pt++) {
//     dist[pt] = 0.0;
//     for (int v=0; v<3; v++)
//       dist[pt] += ( intersections[pt][v] - vertex[v] )*( intersections[pt][v] - vertex[v] );
//     dist[pt] = sqrt(dist[pt]);

//     int index=0;
//     for (int v=0; v<3; v++)
//       if ( dir[v]!=0 ) {
// 	index = v;
// 	break;
//       }
//     if ( (intersections[pt][index] - vertex[index])/dir[index] < 0.0 )
//       dist[pt] = -1*dist[pt];
//   }

//   // (4) Order answer

//   float orderedpts[2][3];
//   float ordereddist[2];
//   if ( dist[0]*dist[1] > 0 ) {
//     // Same sign
//     if ( fabs(dist[0])<fabs(dist[1]) ) {
//       for (int v=0; v<3; v++) {
// 	orderedpts[0][v] = intersections[0][v];
// 	orderedpts[1][v] = intersections[1][v];
//       }     
//       ordereddist[0] = dist[0];
//       ordereddist[1] = dist[1];
//     }
//     else {
//       for (int v=0; v<3; v++) {
// 	orderedpts[0][v] = intersections[1][v];
// 	orderedpts[1][v] = intersections[0][v];
//       }
//       ordereddist[0] = dist[1];
//       ordereddist[1] = dist[0];
//     }
//   }
//   else {
//     // opposite sign
//     if ( dist[0]>0 ) {
//       for (int v=0; v<3; v++) {
// 	orderedpts[0][v] = intersections[0][v];
// 	orderedpts[1][v] = intersections[1][v];
//       }
//       ordereddist[0] = dist[0];
//       ordereddist[1] = dist[1];
//     }
//     else {
//       for (int v=0; v<3; v++) {
// 	orderedpts[0][v] = intersections[1][v];
// 	orderedpts[1][v] = intersections[0][v];
//       }
//       ordereddist[0] = dist[1];
//       ordereddist[1] = dist[0];
//     }
//   }

//   for (int pt=0; pt<2; pt++) {
//     for (int v=0; v<3; v++)
//       intersections[pt][v] = orderedpts[pt][v];
//     dist[pt] = ordereddist[pt];
//   }

//   return SKEP_OK;

// }

// // ===== End of SKEntryPoint( Arrays )  =====================


// int GetBoundaries( int detboundary, float& wallradius, float& capabsz ) {

//   switch (detboundary) {
//   case SKEP_OUTEROD:
//     capabsz = HITKTK*0.5;
//     wallradius = RTKTK;
//     break;
//   case SKEP_INNEROD:
//     capabsz = ZPINTK+ZMED;
//     wallradius = RINTK+RMED;
//     break;
//   case SKEP_ID:
//     capabsz = ZPINTK;
//     wallradius = RINTK;
//     break;
//   default:
//     capabsz = -1;
//     wallradius = -1;
//     return SKEP_BAD_BOUNDARY;
//     break;
//   }

//   return SKEP_OK;

// }

  int IntersectionLineWithCylinder(float vertex[], float dir[], float radius, float absz, float intersections[][3]) {
    
    // Line is given in the form of a ray, which we extend
    // Approach:
    // (1) Determine intersection with circle
    // (2) Determine intersection with cone planes
    // (3) Determine which of the intersections (at most 4) occur
    // Return intersection points via intersections.  

    // INTERSECTIONS WITH WALLS
    bool isvertical = false;
    float wallpts[2][3] = { {0, 0, 0}, {0, 0, 0} };
    if ( dir[0]==0 && dir[1]==0 ) {
      isvertical = true;
    }
    else {
      float origin[3] = {0, 0, 0};
      int closestpt =  IntersectionLineWithCircle(dir, vertex, origin, radius, wallpts);
      if (closestpt<0) {
	return SKEP_NOINTERSECT;
      }
      // Set Z
      int index = 0;
      if ( fabs(dir[0])>fabs(dir[1]) )
	index = 0;
      else
	index = 1;
      for (int pt=0; pt<2; pt++) {
	double s = (wallpts[pt][index] - vertex[index])/dir[index];
	wallpts[pt][2] = vertex[2] + s*dir[2];
      }
    }
    
    
    // INTERSECTIONS WITH PLANES
    bool isparallel = false;
    float cappts[2][3] = { {0, 0, 0}, {0, 0, 0} };
    
    if ( dir[2] == 0.0 ) {
      isparallel = true;
    }
    else {
      for (int plane=-1; plane<2; plane = plane+2 ) {
	double s = (plane*absz - vertex[2])/dir[2];
	int index = (plane + 1)/2; // either 0 or 1
	for (int v=0; v<3; v++)
	  cappts[index][v] = vertex[v] + dir[v]*s;
      }
    }
    
    // PICK THE INTERSECTION POINTS
    // (1) There can only be 2 intersection points
    // (2) Intersection points are within the cylinder
    // (2) Merely loop through intersection points, once I found 2, I am done.
    
    int inters=0;
    for (int pt=0; pt<2; pt++) {
      if ( fabs(wallpts[pt][2]) < absz && isvertical==false) {
	for (int v=0; v<3; v++)
	  intersections[inters][v] = wallpts[pt][v];
	inters++;
      }
    }
    
    if (geotool_verbose)
      std::cout << "Number of intersections after wall test: " << inters << std::endl;
    
    if (inters==2)
      return SKEP_OK;
    
    for (int pt=0; pt<2; pt++) {
      double capr = sqrt( cappts[pt][0]*cappts[pt][0] + cappts[pt][1]*cappts[pt][1] );
      if (capr <= radius && isparallel==false) {
	for (int v=0; v<3; v++)
	  intersections[inters][v] = cappts[pt][v];
	inters++;
      }
    }
    
    if (geotool_verbose)
      std::cout << "Number of intersections after cap test: " << inters << std::endl;
    
    if (inters==2)
      return SKEP_OK;
    else
      return SKEP_NOINTERSECT;
    
  }

// void clearSKIntersection( SKIntersection& skinter ) {

//   skinter.boundary = -1;
//   skinter.intersectionType = SKEP_NOINTERSECT;
//   skinter.distToNearFor = -1;
//   skinter.distToFarFor = -1;
//   skinter.distToNearBack = -1;
//   skinter.distToFarBack = -1;

//   for (int v=0; v<3; v++) {
//     skinter.intersectionNearFor[v] = 0;
//     skinter.intersectionFarFor[v] = 0;
//     skinter.intersectionNearBack[v] = 0;
//     skinter.intersectionFarBack[v] = 0;
//   }

// }

// std::string getIntersectionTypeName( int intersectType ) {

//   std::string intersectName;
//   switch (intersectType) {
//   case SKEP_TWOFORWARD:
//     intersectName =  "SKEP_TWOFORWARD";
//     break;
//   case SKEP_TWOBACKWARD:
//     intersectName =  "SKEP_TWOBACKWARD";
//     break;
//   case SKEP_FORBACK:
//     intersectName =  "SKEP_FORBACK";
//     break;
//   case SKEP_TANGENT:
//     intersectName = "SKEP_TANGENT";
//     break;
//   case SKEP_NOINTERSECT:
//     intersectName = "SKEP_NOINTERSECT";
//     break;
//   default:
//     intersectName = "ERR: Intersection Type Not Defined";
//     break;
//   }

//   return intersectName;
// }

  void getIntersectionPoints( Intersection& skintersect, float intersections[][3] ) {

    // unordered
  
    switch (skintersect.intersectionType) {
    case kForBack:
      for (int v=0; v<3; v++) {
	intersections[0][v] = skintersect.intersectionNearFor[v];
	intersections[1][v] = skintersect.intersectionNearBack[v];
      }
      break;
    case k2Forward:
      for (int v=0; v<3; v++) {
	intersections[0][v] = skintersect.intersectionNearFor[v];
	intersections[1][v] = skintersect.intersectionFarFor[v];
      }
      break;
    case k2Backward:
      for (int v=0; v<3; v++) {
	intersections[0][v] = skintersect.intersectionNearBack[v];
	intersections[1][v] = skintersect.intersectionFarBack[v];
      }
      break;
    case kTangent:
      for (int v=0; v<3; v++) {
	intersections[0][v] = intersections[1][v] = skintersect.intersectionNearFor[v];
      }
      break;
    }//end of switch(intersectionType)
    
  }

// float cwallSK(float pos[]) {

//   float radius = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
//   float absz = fabs(pos[2]);

//   if ( radius<=RINTK && absz <= ZPINTK) {

//     // Easy method for inside the water tank
//     float dr = RINTK-radius;
//     float dz = ZPINTK-absz;
    
//     if (fabs(dr)<fabs(dz))
//       return dr;
//     else
//       return dz;
//   }
//   else {
//     // OUTSIDE THE ID
//     if ( absz <=ZPINTK && radius>RINTK) {
//       // on the sides
//       return RINTK-radius;
//     }
//     else if ( radius<=RINTK && absz>ZPINTK) {
//       // outiside over the caps
//       return ZPINTK-absz;
//     }
//     else if ( radius>RINTK && absz>ZPINTK ) {
//       // outside on the corners
//       float dr = RINTK-radius;
//       float dz = ZPINTK-absz;
//       return -sqrt(dr*dr+dz*dz);
//     }
//   }
  
// }

// void SKSafeVertex( float vertex[], float safevert[] ) {

//   // safe vertex. rounding errors are killing me.
//   // we check if the point is close to an SK boundary.
//   // if it is, we bump it towards the origin
//   double safePhi = atan2( vertex[1], vertex[0] );

//   if ( fabs( fabs(vertex[2])-ZPINTK ) < 0.001 ) {
//     // on the caps                                                                                                                                                                                          
//     if ( vertex[2]>0)
//       safevert[2] = vertex[2] - 0.01;
//     else
//       safevert[2] = vertex[2] + 0.01;
//   }
//   else {
//     safevert[2] = vertex[2];
//   }

//   if (fabs( sqrt(vertex[0]*vertex[0]+vertex[1]*vertex[1])-RINTK ) < 0.001 ) {
//     safevert[0] = vertex[0] - 0.01*cos( safePhi );
//     safevert[1] = vertex[1] - 0.01*sin( safePhi );
//   }
//   else {
//     safevert[0] = vertex[0];
//     safevert[1] = vertex[1];
//   }
//   //std::cout << "safe vertex: " << sqrt( safevert[0]*safevert[0] + safevert[1]*safevert[1] ) << " " << safevert[2] << std::endl;
// }

}//end of geotools namespace
