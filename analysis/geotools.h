#ifndef __geotools_h__
#define __geotools_h__

#include <string>

namespace geotools {

  // store cylinder geometry definition
  struct GeoDef {
    double height;
    double radius;
  };

  // Verbose
  static bool geotool_verbose = false;

  // Intersection Cases
  typedef enum { kUnevaluated, kNoIntersect, k2Forward, k2Backward, kForBack, kTangent } IntersectionCase_t;

  // return status
  enum { SKEP_NOINTERSECT=-1, SKEP_OK };
  
  struct Intersection {
  
    IntersectionCase_t intersectionType;

    float distToNearFor;
    float distToFarFor;
    float distToNearBack;
    float distToFarBack;
    
    float intersectionNearFor[3];
    float intersectionFarFor[3];
    float intersectionNearBack[3];
    float intersectionFarBack[3];
    
    void reset() {
      distToNearFor = distToFarFor = distToNearBack = distToFarBack = -1.0;
      intersectionType = kUnevaluated;
      memset( intersectionNearFor, 0, sizeof(float)*3 );
      memset( intersectionFarFor, 0, sizeof(float)*3 );
      memset( intersectionNearBack, 0, sizeof(float)*3 );
      memset( intersectionFarBack, 0, sizeof(float)*3 );
    }
  };

  // Main function
  int RayCylinderIntersection(float vertex[], float dir[], const GeoDef& geom, Intersection& SKinter);

  // helper functions
  int IntersectionLineWithCylinder(float vertex[], float dir[], float radius, float absz, float intersections[][3]);
  int IntersectionLineWithCircle(float m[], float b[], float o[], float r, float points[][3]);
/*   int RayCylinderIntersection(float vertex[], float dir[], float intersections[][3], float interdist[], GeoDef& geom); */
/*   std::string getIntersectionTypeName( int intersectType ); */
  void getIntersectionPoints( Intersection& skinter, float intersections[][3] );
/*   int GetBoundaries( int detboundary, float& wallradius, float& capabsz ); */
/*   float ClosestPointonLineToPoint(float m[], float b[], float point[], int dim=3); */


}//end of geotools namespace

#endif
