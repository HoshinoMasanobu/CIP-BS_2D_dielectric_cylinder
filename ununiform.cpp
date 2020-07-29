#include "set_variables.h"
#include "iostream"
#include "cmath" 
/*
  不均一なノード番号付け
  NX[0] = NX1, NX[1] = NX1 + NX2, NX[2] = NX1 + NX2 + NX3 
  NY[1] = NY1, NY[1] = NY1 + NY2, NY[2] = NY1 + NY2 + NY3 

  
  NX1,NX2,NX3,NY1,NY2,NY3はVariabls.txtで設定
  
  
  NY[2]----------------------------|
     |         |         |         |
  NY[1]----------------------------|
     |         |         |         |
  NY[0]----------------------------|
     |         |         |         |
    0------------------------------|
     0       NX[0]     NX[1]     NX[2]
           　　　模式図
*/

double  position_x(double ix){
  if(ix < NX[0])
    return DX[0]*ix;
  if(ix < NX[1])
    return RX[0] + DX[1]*(ix - NX[0]);

  return RX[0] + RX[1] + DX[2]*(ix - NX[1] );
}
double  position_y(double iy){
  if(iy < NY[0])
    return DY[0]*iy;
  if(iy < NY[1])
    return RY[0] + DY[1]*(iy - NY[0]);

  return RY[0] + RY[1] + DY[2] * (iy - NY[1]);
}

double position(int xy,double i){
  if(xy == 0)return position_x(i);
  if(xy == 1)return position_y(i);
  

  exit(0);
}










