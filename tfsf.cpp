#include "gaussian_pulse.h"
#include "set_array.h"
#include "set_variables.h"
#include "inner_product.h"
#include "make_coefficient_matrix.h"
#include "ununiform.h"
#include "dielectrics.h"
#include <cmath>

const
int
corner_val[4][2] = {{-1,-1}, {1,-1}, {1,1}, {-1,1}},
  update_xy[4][2] = {{1,0},{0,1},{-1,0},{0,-1}};
/*
  関数set_TFSFと関数update_incidentの参照順序
  TFSF境界を内側、外側の順に二周する。
  
  |<----------
  ||<--------A
  ||        A|
  |V        ||
  V-------->||
  ---------->|


  CIP-BS TFSF境界利用上の注意点
  1)急峻な波形は入射しない。
  2)境界上のセルは正方形にする。
*/
class tfsf_vector{
public:
  int
  ix,
    iy,
    x_order,y_order,
    row,
    in_out;
};
tfsf_vector **iterator_tfsf;
void set_plane_wave();

int max_i_circum[3];

void set_TFSF(sparse_matrix ***Coef_Mat){
  
  
  set_plane_wave();
  
  int half[2] = {0,0};

  for(int xyz = 0; xyz <= 2; xyz ++){//E_xyz or H_xyz

    int i_circum = 0;
      
    int E_or_H = 1;

    if(xyz == 2)E_or_H = 0;
    
    int
      RHS_xyz_minus,
      RHS_xyz_plus; 
    relative_index(xyz,&RHS_xyz_minus,&RHS_xyz_plus);
    int RHS_xyz[3] = {xyz,RHS_xyz_plus,RHS_xyz_minus};
      
    for(int list_index = 0; list_index <= 2; list_index++){
     
      int
        target_xyz,
        target_E_or_H = E_or_H;
      if(list_index >= 1)target_E_or_H = !E_or_H;
    
      std::vector <T> list;
    
      if(list_index == -1){
        target_xyz = xyz;
        list.resize(TFSF_N * 3);
      }else{
        
        target_xyz = RHS_xyz[list_index];
        list.resize(TFSF_N * 3);
      }
      if(E_or_H == 1){
        if(list_index > 0 && target_xyz != 2)continue;
      }
  
      
      //試験関数の中心決定
      //TFSF境界の内側0か外側1
      for(int in_out = 0; in_out <= 1; in_out ++){

        //TFSF境界を一周する
        for(int corner = 0; corner <= 3; corner ++){
       
          int
            next_corner = (corner + 1)%4;

          //cornerからcorner+1まで
          int ix = TFSF_corner[corner][0] + in_out*corner_val[corner][0];
          do{
            
            int iy = TFSF_corner[corner][1] + in_out*corner_val[corner][1];
            do{

              //試験関数の導関数の次数 決定
              for(int test_x_order = 0; test_x_order <= K; test_x_order++){
                for(int test_y_order = 0; test_y_order <= K; test_y_order++){
                  
                  int row = Node( xyz, ix, iy,  test_x_order, test_y_order);
                  if(row == -1)continue;
                  if(list_index == 0)i_circum++;
                  //基底関数の中心(試験関数の中心からの距離) 決定
                  for(int x_val = -1; x_val <= 1; x_val++){
                    for(int y_val = -1; y_val <= 1; y_val++){
                      if(x_val == 0 && y_val == 0)continue;//
                      
                      /*試験関数の中心ノードがTF領域なら隣接ノードからSF領域を、
                        中心ノードがSF領域なら隣接ノードからTF領域を検出する。
                       */                      
                      int tfsf = 1;//SF外
                      if(ix + x_val >= start_x && ix + x_val <= end_x
                         &&
                         iy + y_val >= start_y && iy + y_val <= end_y)
                        tfsf = 0;//TF中
                      if(tfsf == in_out)continue;
                      

                      //基底関数の導関数の次数 決定
                      for(int basis_x_order = 0; basis_x_order <= K; basis_x_order++){
                        for(int basis_y_order = 0; basis_y_order <= K; basis_y_order++){
                          
                          int column
                            = Node(target_xyz, ix + x_val, iy + y_val
                                   ,basis_x_order, basis_y_order);
                          
                          if(column != -1){
                            
                            list.push_back(
                                           T(
                                             row,column,
                                             element(
                                                     E_or_H,
                                                     xyz,
                                                     target_E_or_H,
                                                     target_xyz,
                                                     half,
                                                     test_x_order,
                                                     test_y_order,
                                                     basis_x_order,
                                                     basis_y_order,
                                                     ix,iy,
                                                     x_val,
                                                     y_val
                                                     )
                                    
                                             )
                                           );
      
                          }
                        }
                      }
                    }
                  }
                 
                }
              }
         
         
              iy += update_xy[corner][1];
            }while( iy != TFSF_corner[next_corner][1]  + in_out*corner_val[next_corner][1]);
            
            ix += update_xy[corner][0];
          }while(ix != TFSF_corner[next_corner][0] + in_out*corner_val[next_corner][0]);
        }
      }
     
      
      (*Coef_Mat)[xyz][list_index ]
        .resize(row_size[xyz],row_size[target_xyz]);
      
      (*Coef_Mat)[xyz][list_index ]
        .setFromTriplets(
                         list.begin(),
                         list.end() );
      
      (*Coef_Mat)[xyz][list_index ]
        .makeCompressed();
      
    }
    max_i_circum[xyz] = i_circum;
  }
  iterator_tfsf = new tfsf_vector* [3];
  for(int xyz = 0; xyz <= 2; xyz++){
    iterator_tfsf[xyz] = new tfsf_vector[max_i_circum[xyz]];
    int i_circum = 0;
    for(int in_out = 0; in_out <= 1; in_out ++){
      for(int corner = 0; corner <= 3; corner ++){
      
        int
          sign = 1.0,
          next_corner = (corner + 1)%4;
      
        if(in_out == 0)sign = -1.0;
        int ix = TFSF_corner[corner][0] + in_out*corner_val[corner][0];
        do{
          int iy = TFSF_corner[corner][1] + in_out*corner_val[corner][1];
          do{
          
            //試験関数の導関数の次数 決定
            for(int x_order = 0; x_order <= K; x_order++){
              for(int y_order = 0; y_order <= K; y_order++){
              
                int row = Node( xyz, ix, iy,  x_order, y_order);
                if(row == -1)continue;
                  
                iterator_tfsf[xyz][i_circum].ix = ix;
                iterator_tfsf[xyz][i_circum].iy = iy;
                iterator_tfsf[xyz][i_circum].x_order
                  = x_order;
                iterator_tfsf[xyz][i_circum].y_order
                  = y_order;
                iterator_tfsf[xyz][i_circum].row = row;
                iterator_tfsf[xyz][i_circum].in_out = in_out;
                i_circum ++;
                  
              
             
              }
            }
          
            iy += update_xy[corner][1];
          }while( iy != TFSF_corner[next_corner][1]   + in_out*corner_val[next_corner][1]);
          ix += update_xy[corner][0];
        }while(ix != TFSF_corner[next_corner][0] + in_out*corner_val[next_corner][0] );
      }
    }
  
  }
  
}

int
EH[3] ={1,1,0};
double
phi_i = 0.0,
  cosphi = cos(phi_i),//for x
  sinphi = sin(phi_i),//for y
  k_array[3] ={sinphi*std::sqrt(eps/mu), -cosphi*std::sqrt(eps/mu),1.0},
  delta = 1.0/20.0,
  kk = M_PI/outer_r,
  kx,
  ky,
  omega = kk*c;

void set_plane_wave(){
  
  kx = M_PI/32.0/std::sqrt(delta*delta + delta*delta)*cosphi;
  ky = M_PI/32.0/std::sqrt(delta*delta + delta*delta)*sinphi;
  
}

double plane_wave_gaussian
(int xyz, int n, int ix ,int iy, int x_order,int y_order){
  
  return
    k_array[xyz]*
    Gaussian_Pulse(
                   c*DT*(n - 0.5*EH[xyz] ),
                   cosphi*position_x(ix) + sinphi*position_y(iy),
                   32.0*(delta*delta + delta*delta),
                   x_order + y_order )
    *std::pow(cosphi,x_order )
    *std::pow(sinphi,y_order );
}

double pulse_Gauss_sin(int xyz, int n, int ix ,int iy,
                       int gauss_x_order, int gauss_y_order,
                       int sin_x_order,   int sin_y_order,
                       int max_x_order,   int max_y_order){
  
  if( (sin_x_order + gauss_x_order) > max_x_order // || gauss_x_order > max_x_order ||
     || (sin_y_order + gauss_y_order) > max_y_order )// || gauss_y_order > max_y_order )
    return 0.0;
  if((gauss_x_order + gauss_y_order + sin_x_order + sin_y_order)
     == (max_x_order + max_y_order))  {
   
    return
      2.0*k_array[xyz]
      *std::pow(cosphi,gauss_x_order)
      *std::pow(sinphi,gauss_y_order)
      *Gaussian_Pulse(
                      c*DT*(n + 0.5*EH[xyz] ),
                      cosphi *  position_x(ix) 
                      +
                      sinphi *  position_y(iy)
                      - (cosphi * DX[0] *pml_L + sinphi *DY[0] *pml_L)
                      + std::sqrt(2.0)*1.5  - 2.5,
                     
                      64.0*(delta * delta + delta * delta),
                      gauss_x_order + gauss_y_order )
      
      *std::pow( kk, sin_x_order + sin_y_order)
      *std::pow( cosphi,sin_x_order)
      *std::pow( sinphi,sin_y_order)
      *std::sin(
                kk
                *(
                  cosphi * position_x(ix) + sinphi * position_y(iy) 
                  - c*DT*(n + 0.5*EH[xyz])
                  - (cosphi * DX[0] *pml_L + sinphi *DY[0] *pml_L)
                  + std::sqrt(2.0)*1.5  - 2.5
                 
                  )
                + M_PI/2.0*(sin_x_order + sin_y_order)
                );
  }
  
  if(max_x_order > sin_x_order + gauss_x_order){

    return
      pulse_Gauss_sin(xyz, n, ix, iy,
                      gauss_x_order + 1, gauss_y_order,
                      sin_x_order, sin_y_order,
                      max_x_order, max_y_order)
      + 
      pulse_Gauss_sin(xyz, n, ix, iy,
                      gauss_x_order, gauss_y_order,
                      sin_x_order + 1, sin_y_order,
                      max_x_order, max_y_order);
    
  }else if(max_y_order > sin_y_order + gauss_y_order){
    return
      pulse_Gauss_sin(xyz, n, ix, iy,
                      gauss_x_order, gauss_y_order + 1,
                      sin_x_order, sin_y_order,
                      max_x_order, max_y_order)
      + 
      pulse_Gauss_sin(xyz, n, ix, iy,
                      gauss_x_order, gauss_y_order,
                      sin_x_order, sin_y_order + 1,
                      max_x_order, max_y_order);
  }
}

double plane_wave_gauss_sin
(int xyz, int n, int ix ,int iy, int x_order,int y_order){
 
  return
    pulse_Gauss_sin( xyz, n, ix, iy ,
                     0, 0, 0, 0, x_order, y_order);
  
}

void incident_update( int xyz, int n, Eigen::VectorXd *incident){
  double sign[2] = {-1.0,1.0};
  for(int i = 0; i < max_i_circum[xyz]; i++){
    (*incident)(iterator_tfsf[xyz][i].row) =
      sign[iterator_tfsf[xyz][i].in_out]
      *plane_wave_gauss_sin(xyz,n,
                            iterator_tfsf[xyz][i].ix,
                            iterator_tfsf[xyz][i].iy,
                            iterator_tfsf[xyz][i].x_order,
                            iterator_tfsf[xyz][i].y_order) ;
   
   
  }
  
}



















































