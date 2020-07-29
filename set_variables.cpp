#include "set_variables.h"
#include "dielectrics.h"


  /*解析空間の離散化数毎の領域の長さをRX[n],RY[n]に格納する

              0     NX[0]       NX[1]    NX[2] 
              | RX[0] |    RX[1]   | RX[2] |   
  NY[2] -----------------------------------
        RY[0] |       |            |       |
  NY[1] -----------------------------------
        RY[1] |       | Dielectric |       |
  NY[0] -----------------------------------
        RY[2] |       |            |       |
  0     -----------------------------------

    RX[n],RY[n]はそれぞれNX[n] - NX[n-1],NY[n] - NY[n-1]個のセルに、
    Variables.txt内のNXn,NYn個のセルに離散化される。
  */

int
K,//CIP-BSの次数
  NT,//総時間ステップ数であるが使わない
  NX[3],
  NY[3],
  
  //PMLの層数
  pml_L,
  
//TFSF境界の頂点の座標
  start_x ,
  end_x ,
  start_y ,
  end_y ,
  
  TFSF_N,
  TFSF_corner[4][2];

 
double
  C,//クーラン数 
  RX[3],//x
  RY[3],//y
  DX[3],
  DY[3],
  DT,

//フーリエ変換の観測位置である長方形の頂点の座標
  FT_x_start,
  FT_x_end,
  FT_y_start,
  FT_y_end;

void set_variables(){
 
  std::ifstream ifs("Variables.txt");
  std::string dummy;
  ifs >> dummy >> K;//CIP-BS法の次数
  ifs >> dummy >> NT;//ステップ数であるが使わない
  ifs >> dummy >> RX[2];//解析領域のX方向の長さ 
  ifs >> dummy >> RY[2];//解析領域のY方向の長さ 
  
  ifs >> dummy >> C ;//クーラン数
  ifs >> dummy >> NX[0];//誘電体の無いx = 0m からx = 2mを何セルに離散化するか
  ifs >> dummy >> NX[1];//誘電体のあるx = 2m からx = 3mを何セルに離散化するか
  ifs >> dummy >> NX[2];//誘電体の無いx = 3m からx = 5mを何セルに離散化するか
  
  C = 0.4 * NX[1] / 30.0;//時間ステップをC=0.4,NX[1] = 30の時に等しくする様にクーラン数を設定する.
  
  //NXと同じ
  ifs >> dummy >> NY[0];
  ifs >> dummy >> NY[1];
  ifs >> dummy >> NY[2];

 
  //フーリエ変換を観測する長方形の頂点の座標を設定
  ifs >> dummy >>  FT_x_start ;
  ifs >> dummy >>  FT_x_end ;
  ifs >> dummy >>  FT_y_start ;
  ifs >> dummy >>  FT_y_end ;

  /*解析空間の離散化数毎の領域の長さをRX[n],RY[n]に格納する
    　
          | RX[0] |   RX[1]   | RX[2] |   
    ----------------------------------
    RY[0] |       |           |       |
    ----------------------------------
    RY[1] |       |Dielectrics|       |
    ----------------------------------
    RY[2] |       |           |       |
    ----------------------------------

    RX[2]に領域の全体の長さが格納されていたが
    以下ではRX[0] + RX[1] + RX[2]が全体の長さとなる。 
  */
  RX[1] = 2.0*outer_r;
  RY[1] = 2.0*outer_r;
  for(int i = 0; i <= 2;i += 2){
    RX[i] = (RX[2] - RX[1])/2.0;
    RY[i] = (RY[2] - RY[1])/2.0;

  }
  
  pml_L = 12;

  for(int i = 0; i<=2 ; i++){
    
    DX[i] = RX[i]/(double)NX[i];
    DY[i] = RY[i]/(double)NY[i];
    if(i != 0){
      NX[i] += NX[i-1];
      NY[i] += NY[i-1];
    }
    std::cout << RX[i] <<" " << NX[i]<< " " << DX[i] <<std::endl;
    std::cout << RY[i] <<" " << NY[i]<< " " <<  DY[i] <<std::endl;
  }
  
  NX[2] += pml_L;
  NY[2] += pml_L;
  
  for(int i = 0; i <= 2; i ++){
    
    NX[i] += pml_L;
    NY[i] += pml_L;
  }
  for(int i = 0; i <= 2; i+= 2){
    RX[i] += DX[i] * pml_L;
    RY[i] += DY[i] * pml_L;
  }
  
  
  DT = C/c/std::sqrt(1.0/DX[1]/DX[1] + 1.0/DY[1]/DY[1]);

  //TFSFは誘電体の直径より1セル大きい辺を持つ長方形で囲む
  start_x = NX[0] - 1;
  end_x =   NX[1] + 1;
  start_y = NY[0] - 1;
  end_y =   NY[1] + 1;

  //TFSFの頂点の座標をiterator用配列に格納
  TFSF_N = 2 * ( end_x - start_x ) + 2 * ( end_y - start_y);
  int tmp_corner[4][2] =
    {{start_x,start_y},{end_x,start_y},
     {end_x,end_y},    {start_x,end_y}};
  for(int corner = 0; corner <=3;corner++){
    for(int xy = 0; xy <= 1; xy++){
      TFSF_corner[corner][xy] = tmp_corner[corner][xy];
    }
  }
}














