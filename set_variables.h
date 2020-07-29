#ifndef GLOBAL_VARIABLES
#define GLOBAL_VARIABLES
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
 
constexpr double 
c { 3.0e8 },//光速
  mu { 4.0 * M_PI *1.0e-7 },//μ、透磁率
  eps { 1.0 / (mu * c * c) };//ε、誘電率
extern int
K,//CIP-BS法の次数
  NT,//時間ステップ
  NX[3],//不均一分割による境界のノード番号
  NY[3],//不均一分割による境界のノード番号
  pml_L,//PMLの層数

//TFSF境界の座標
  start_x ,
  end_x ,
  start_y ,
  end_y ,
  TFSF_N,
  TFSF_corner[4][2];
extern double
  C,//クーラン数 
  RX[3],//領域の長さ x
  RY[3],//領域の長さ y
  DX[3],//不均一分割
  DY[3],//
  DT,
  FT_x_start,
  FT_x_end,
  FT_y_start,
  FT_y_end;
void set_variables();
#endif





