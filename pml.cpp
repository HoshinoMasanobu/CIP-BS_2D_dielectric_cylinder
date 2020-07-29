#include "set_array.h"
#include "set_variables.h"
#include "ununiform.h"

/*sigmaについて
  PMLでは導電率と導磁率を導入するが
  NPMLではsigmaがそれらにあたるものである.
 */
int
  **pml_start,
  **pml_len;
double sigma_max;
int
M = 4;
double 
R = 1.0e-6;
double Delta[2] = {},N[2]={},
  ***sigma_array;

/*変数x_area,y_areaの意味
  11が解析領域、他はPML
         x_area
           0  1  2
         --------- 
       0| 00 01 02
y_area 1| 10 11 12
       2| 20 21 22
　　

 */
//sigmaの最大値
void set_sigma_max(){
  sigma_max = -  ( M + 1.0)*c/(2.0*pml_L*DX[0])*std::log(std::abs(R));
  
}
//sigmaとその導関数
double sigma(int xy/*成分 x = 0, y = 1*/, int i /*node*/, int diff/*導関数の次数*/){
    
  double P = 1;
  for(int i = 0 ; i < diff; i++){
    P *= (M - i);
  }
  if(i < pml_start[xy][1]){
    return
      sigma_max
      *P*std::pow(-1.0/pml_len[xy][0]*Delta[xy],diff )
      *std::pow(
                
                (pml_len[xy][0]*Delta[xy] - position(xy,i))
                /(pml_len[xy][0]*Delta[xy]),M - diff);
  }
  if(i >= pml_start[xy][2]){
    return
      sigma_max
      *P*std::pow(1.0/pml_len[xy][2]*Delta[xy],diff )
      *std::pow(
                
                ( position(xy,i) - position(xy,N[xy] - pml_len[xy][2]))
                /(pml_len[xy][2]*Delta[xy]),M - diff);
  }
  
  return 0.0;
}

//境界のノード番号を設定
void set_node_num_pml(){
  set_sigma_max();
  Delta[0] = DX[0];
  Delta[1] = DY[0];
  N[0] = NX[2];
  N[1] = NY[2];
  
  pml_start = new int *[2]; 
  pml_len = new int *[2];
  for(int xy = 0; xy <= 1; xy++){
    pml_start[xy] = new int [3];
    pml_len[xy] = new int [3];
    
  }
  pml_start[0][0] = 0;
  pml_start[1][0] = 0;

  pml_len[0][0] = pml_L;
  pml_len[1][0] = pml_L;
  pml_len[0][2] = pml_L;
  pml_len[1][2] = pml_L;
    
  pml_start[0][1] = pml_start[0][0] + pml_len[0][0] + 1;
  pml_start[1][1] = pml_start[1][0] + pml_len[1][0] + 1;
     
  pml_start[0][2] = NX[2] - pml_len[0][2];
  pml_start[1][2] = NY[2] - pml_len[1][2];

  pml_len[0][1] =  pml_start[0][2] - pml_start[0][1] - 1;
  pml_len[1][1] =  pml_start[1][2] - pml_start[1][1] - 1;


  int N[2] ={NX[2],NY[2]};
  //sigmaを毎回計算すると遅くなるため配列に格納しておく。
  sigma_array = new double** [2];
  for(int xy = 0; xy <=1; xy++){
    sigma_array[xy] = new double* [N[xy] + 1];
    for(int i = 0; i <= N[xy]; i++){
      sigma_array[xy][i] = new double [K+1];
      for(int order = 0; order <= K; order ++){
        sigma_array[xy][i][order] = sigma(xy,i,order); 
      }
    }
  }
   
}

void update_E_pml(Eigen::VectorXd **EM_density/*電磁界のポインタ*/,
                  Eigen::VectorXd *EM_past/*１時間ステップだけ過去の電磁界の値*/,
                  Eigen::VectorXd **pmlEz/*電界の補助場*/){
  for(int xy = 0; xy <= 1; xy ++){
    Eigen::VectorXd
        pmlEz_past = (*pmlEz)[xy];
    
    for(int x_area = 0; x_area <= 2; x_area ++){
      for(int y_area = 0; y_area <= 2; y_area ++){  
                
        for(int ix = pml_start[0][x_area]; ix <= pml_start[0][x_area] + pml_len[0][x_area]; ix++){
          for(int iy = pml_start[1][y_area]; iy <= pml_start[1][y_area] + pml_len[1][y_area]; iy++){
     
            int i[2]={ix,iy}; 
            double
              sigma0 = sigma_array[xy][i[xy]][0]/*sigma(xy,i[xy],0)*/,
              sigma1 = sigma_array[xy][i[xy]][1]/*sigma(xy,i[xy],1)*/;
     
            for(int x_order = 0; x_order <= K; x_order++){
              for(int y_order = 0; y_order <= K; y_order++){
                int order[2] ={x_order,y_order};
                
                int node = Node( 2, ix, iy ,x_order ,y_order);
                if(node != -1){
                 
                  (*pmlEz)[xy](node) =
                    (2.0 - DT*sigma0)/(2.0 + DT*sigma0)*pmlEz_past(node)
                    + 2.0/(2.0 + DT*sigma0)*((*EM_density)[2](node) - EM_past[2](node));
                  
                  if(order[xy] == 1){
                    order[xy] = 0;
                    int node_RHS = Node(2, ix, iy, order[0] ,order[1]);
                    
                    if(node_RHS != -1){
                      
                      (*pmlEz)[xy](node ) +=
                        - DT *sigma1/(2.0 + DT*sigma0)
                        * ((*pmlEz)[xy](node_RHS) + pmlEz_past(node_RHS));
                      
                    }
                  }
                }
              }
            }
            

          }
        }
      }
    }

  }

}

void update_H_pml(Eigen::VectorXd **EM_density,Eigen::VectorXd *EM_past,
                  Eigen::VectorXd **pmlH){
  for(int xy = 0; xy <= 1; xy ++){
    Eigen::VectorXd
      pmlH_past = (*pmlH)[xy];
    for(int x_area = 0; x_area <= 2; x_area ++){
      for(int y_area = 0; y_area <= 2; y_area ++){  
       
        for(int ix = pml_start[0][x_area]; ix <= pml_start[0][x_area] + pml_len[0][x_area]; ix++){
          for(int iy = pml_start[1][y_area]; iy <= pml_start[1][y_area] + pml_len[1][y_area]; iy++){
            int i[2]={ix,iy};
            double
              sigma0 = sigma_array[!xy][i[!xy]][0]/*sigma(!xy,i[!xy],0)*/,
              sigma1 = sigma_array[!xy][i[!xy]][1]/*sigma(!xy,i[!xy],1)*/;
            
            for(int x_order = 0; x_order <= K; x_order++){
              for(int y_order = 0; y_order <= K; y_order++){
                int order[xy] ={x_order,y_order};
                
                int node = Node( xy, ix, iy ,x_order ,y_order);
                if(node != -1){
                  (*pmlH)[xy](node) =
                    (2.0 - DT*sigma0)/(2.0 + DT*sigma0)*pmlH_past(node)
                    + 2.0/(2.0 + DT*sigma0)*((*EM_density)[xy](node) - EM_past[xy](node));

                  if(order[!xy] == 1){
                    order[!xy] = 0;
                    int node_RHS = Node(xy, ix, iy, order[0] ,order[1]);
                    if(node_RHS!=-1)
                      (*pmlH)[xy](node) +=
                        - DT *sigma1/(2.0 + DT*sigma0)
                        * ((*pmlH)[xy](node_RHS) + pmlH_past(node_RHS));
                  }
                }
              }
            }
          }
        }
      }
    }
  }
 
}






















