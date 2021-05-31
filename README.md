# Clothoid
一个计算缓和曲线上任意一点坐标的demo.
   /**         
   *  Calculate segment angle when splitting clothoid. relative to north
   *  a = (Re-Rs)/L         
   *  angle = Rs * L + 1/2 * a * L*L * 180/ PI;
   *  Param a curvature change rate         
   *  Param Rs start curvature         
   *  Param L length of clothoid, uint:m         
   *  Result angle, uint: degree         
   *  double a = (endCurvature - startCurvature) / length;
   *  double angle_diff = startCurvature * length + (a * length * length) / 2;
   *  double angle = startAngle + angle_diff;
   **/
   
   
