

#ifndef SPK_APP_POTTSWELD_POOL_SHAPE_H
#define SPK_APP_POTTSWELD_POOL_SHAPE_H

#include<vector>

using std::vector;

namespace weld {

namespace pool_shape {

   enum ShapeType {undefined=-1, ellipse=0, teardrop=1};

   class PoolShape {
      public:
         virtual double distance(const double *XYZ) const = 0;
         PoolShape() {}
         virtual ~PoolShape() {}
   };

}

}

#endif
