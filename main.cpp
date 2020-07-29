#include <iostream>
#include <Eigen/Sparse>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <chrono>
#include "inner_product.h"
#include "set_variables.h"
#include "make_coefficient_matrix.h"
#include "set_array.h"
#include "initialize.h"
#include "gaussian_pulse.h"
#include "pml.h"
#include "tfsf.h"
#include "ununiform.h"
#include "fourie_transform.h"
#include "interpolation.h"
#include "dielectrics.h"

int main(){

    
  std::chrono::system_clock::time_point start_time =
    std::chrono::system_clock::now();

  
  set_variables();//変数設定
  basis_set();//基底関数の各項の係数計算
  set_node_number();//電界や磁界の配列番号設定
  set_node_num_pml();//pmlに関する配列番号設定

  
  //メモリ確保
  sparse_matrix
    **RHS_Coef_Mat,/*
                     係数行列格納
                     RHS_Coef_Mat[Unknouwn itelator 0 = Hx,1 = Hy,2 = Ez][RHS itelator]
                     
                     RHS_Coef_Mat[0 = Hx][0 = Matrix for Hx, 1 = Matrix for Ez, (2 = Matrix for Hy)]
                     RHS_Coef_Mat[1 = Hy][0 = Matrix for Hy,(1 = Matrix for Hx), 1 = Matrix for Ez]
                     RHS_Coef_Mat[2 = Ez][0 = Matrix for Ez, 1 = Matrix for Hy,  2 = Matrix for Hx]
                     ※()は存在しない係数行列
                   */
    **TFSF_Coef;//TFSF様の係数格納 配列の意味は同上
  Solver
    *solver;//solver[0]はHx,solver[1]はHy,solver[2]はEz様のソルバ
  Eigen::VectorXd
    //入射場
    *incidentEz,//[0]に時刻 (n+1)*DT,[1]に時刻 n*DT のEzを格納
    **incidentH,//incidentH[xy][時刻]

    //電磁界格納
    *EM_density,//EM_density[0] = Hx, EM_density[1] = Hy, EM_density[2] = Ez
    *EM_old,//1時刻前の電磁場
    
    //Nearly PMLに関する補助場
    *pmlEz,//pmlEz[0] = Ezx, pmlEz[1] = Ezy
    *pmlH;//pmlH[0] = Hx,pmlH[1] = Hy
  
  //係数
  RHS_Coef_Mat = new sparse_matrix *[3];
  TFSF_Coef = new sparse_matrix *[3];
  solver = new Solver [3];
  
  //電磁界の値
  EM_density = new Eigen::VectorXd [3];
  EM_old = new Eigen::VectorXd [3];

  for(int xyz = 0; xyz <= 2; xyz ++){
    RHS_Coef_Mat[xyz] = new sparse_matrix [3];
    TFSF_Coef[xyz] = new sparse_matrix [3];
  
    EM_density[xyz].resize(row_size[xyz]);
    EM_old[xyz].resize(row_size[xyz]);
    
  }
  
  pmlEz = new Eigen::VectorXd [2];
  pmlH = new Eigen::VectorXd [2];

  incidentEz = new Eigen::VectorXd [2];
  incidentH = new Eigen::VectorXd *[2];
  
  for(int i = 0; i <= 1; i ++){
    int xy = i;
    pmlEz[xy].resize(row_size[2]) ;
    pmlH[xy].resize(row_size[xy]) ;
    incidentH[xy] = new Eigen::VectorXd [2];

    int new_or_old_Ez = i;
    incidentEz[new_or_old_Ez].resize(row_size[2]);

    for(int new_or_old_H = 0; new_or_old_H <=1;new_or_old_H ++){
      incidentH[xy][new_or_old_H].resize(row_size[xy]);
      
    }
  }
  
  //入射波初期値代入
  incident_update(0,0,&incidentH[0][0]);
  incident_update(0,-1,&incidentH[0][1]);
  
  incident_update(1,0,&incidentH[1][0]);
  incident_update(1,-1,&incidentH[1][1]);
  
  incident_update(2,0,&incidentEz[0]);
  incident_update(2,-1,&incidentEz[1]);
  
  
  //行列生成
  Make_Matrix(&solver,&RHS_Coef_Mat);
  set_TFSF(&TFSF_Coef);
  
  initialize(&EM_density);
  

  //fourie transformに係るパラメータ等設定
  set_fourie_trans_parametar_circle();
  set_fourie_trans_parametar();
  
   std::ofstream ofs_tr("Ez_time_response" + std::to_string(pml_L) + ".dat");

   /*
   std::chrono::system_clock::time_point ite_start =
     std::chrono::system_clock::now();
   */
   double maxT = 6.28539e-08; 
   for(int n = 0; DT*n < maxT; n++){
     
     std::cout<<DT*n/maxT*100.0<<"% " << n <<std::endl;
      
     //Fourie Transform
     update_FT_circle(&EM_density,n);
     update_FT(&EM_density,n);

     //Animaion
     /*
     if(n%10 == 0){   
       
       double
         LX = RX[0] + RX[1] + RX[2],
         LY = RY[0] + RY[1] + RY[2];
       
       std::ofstream ofs2( "ANIM/interpolate_Ez"+std::to_string(n)+".dat");
       for(int ix = 0; ix <= 100; ix++){
        for(int iy = 0; iy <= 100; iy++){
          
          double
            x = LX/100.0 * ix,
            y = LY/100.0 * iy;
          
          ofs2 << x << " " << y << " "
               << EM_interpolated(&EM_density, 2,  x,  y) <<std::endl;
        }
        ofs2 << std::endl;
       }
       ofs2.close();
       
     }
     */

     
     //電界計算
     
     int
       xyz = 2,
       R_xyz_plus,R_xyz_minus;
     
     EM_old[xyz] = EM_density[xyz];
    relative_index(xyz,&R_xyz_minus,&R_xyz_plus);
   
    EM_density[xyz] 
      = RHS_Coef_Mat[xyz][0]* ( EM_density[xyz] )
      + RHS_Coef_Mat[xyz][1]* ( pmlH[R_xyz_plus]  )
      + RHS_Coef_Mat[xyz][2]* ( pmlH[R_xyz_minus] )

      
      //TFSFについての項
      + TFSF_Coef[xyz][0] * (incidentEz[1] - incidentEz[0])
      + TFSF_Coef[xyz][1] * incidentH[R_xyz_plus][1]
      + TFSF_Coef[xyz][2] * incidentH[R_xyz_minus][1];
      

    
    EM_density[xyz] = solver[xyz].solve(EM_density[xyz]);//Ez更新
     
    update_E_pml(&EM_density, EM_old, &pmlEz);//補助場pmlEz更新

    //入射波更新
    incidentEz[1] = incidentEz[0];//時刻 n*DT の磁界を過去の磁界として格納
    incident_update(2,n + 1,&incidentEz[0]);//時刻(n + 1)*DTの入射電界を計算
    
    //磁界計算
    
    for( xyz = 0; xyz <= 1; xyz++){
      EM_old[xyz] = EM_density[xyz];
      int  R_xyz_plus,R_xyz_minus;
      relative_index(xyz,&R_xyz_minus,&R_xyz_plus);
      //磁界についての項
      EM_density[xyz]
        = RHS_Coef_Mat[xyz][0]*EM_density[xyz]

        + TFSF_Coef[xyz][0] * (incidentH[xyz][1] - incidentH[xyz][0]);

      
      //Ezについての項
      if(R_xyz_plus == 2){

        EM_density[xyz] -=
          (
           RHS_Coef_Mat[xyz][1]*pmlEz[!xyz]
           
           + TFSF_Coef[xyz][1] * incidentEz[1]
           );
        
      }else if(R_xyz_minus == 2){

        EM_density[xyz] -=
          (
           RHS_Coef_Mat[xyz][2]*pmlEz[!xyz]
           
           + TFSF_Coef[xyz][2] * incidentEz[1]
           );
        
      }
      
      EM_density[xyz] = solver[xyz].solve(EM_density[xyz]);//磁界更新
      
      //入射波更新
      incidentH[xyz][1] = incidentH[xyz][0];//時刻( n + 1/2 )*DT の磁界を過去の磁界として格納
      incident_update(xyz,n + 1,&incidentH[xyz][0]);//時刻(n + 1 + 1/2)*DTの入射磁界を計算
      
    }
        
    update_H_pml( &EM_density, EM_old, &pmlH);//補助場　pmlH　更新
  }
   /*
   std::chrono::system_clock::time_point end_time =
     std::chrono::system_clock::now();
   
   double all_elapsed = std::chrono::duration_cast <std::chrono::milliseconds>
     (end_time - start_time).count(),
    ite_elapsed = std::chrono::duration_cast <std::chrono::milliseconds>
     (end_time - ite_start).count();
   
   std::ofstream ofs_time("calctime.dat",std::ios::app);
   ofs_time<< (NX[0] - pml_L) /2.0<<" "<<all_elapsed<< " "<< ite_elapsed<<std::endl;
  */

   //Fourie transform 結果出力
   output_result_FourieTrans_circle();
   output_result_FourieTrans();
  
  
  return 0;
}

































































































































































































