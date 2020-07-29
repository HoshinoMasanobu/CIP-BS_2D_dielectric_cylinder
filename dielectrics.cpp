#include "set_variables.h"
#include "dielectrics.h"
#include <iostream>
#include <cmath>

double dielectrics(double x, double y){
 
  x -= (RX[2] + RX[1] + RX[0])/2.0; y -= (RY[2] + RY[1] + RY[0])/2.0;
  double r = std::sqrt(x*x + y*y);
  double theta = atan2(y,x);
  if(  r <= outer_r){
    return eps*eps_r;
  }
  
  return eps;
}

































