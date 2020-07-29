#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
double Gaussian_Pulse(double eta,double xi,double sigma,int order){
  if(order == 0) return std::exp (
                                  -(
                                    std::pow(xi - eta, 2.0)
                                   )/(2.0*std::pow(sigma, 2.0))
                                  );
  if(order == 1)return 
                  -(xi - eta)
                  /(std::pow(sigma, 2.0))
                  *
                  std::exp (
                            -(
                              std::pow(xi - eta, 2.0)
                              )/(2.0*std::pow(sigma, 2.0))
                            );
  if(order == 2) return
                   -1.0/std::pow(sigma, 2.0)
                   *
                   std::exp (
                             -(
                               std::pow(xi - eta, 2.0)
                               )/(2.0*std::pow(sigma, 2.0))
                             )
                   +
                   std::pow(xi - eta, 2.0)
                   /std::pow(sigma, 4.0)
                   *
                   std::exp (
                             -(
                               std::pow(xi - eta, 2.0)
                               )/(2.0*std::pow(sigma, 2.0))
                             );
  if(order == 3)return (xi - eta)
                  /std::pow(sigma, 4.0)
                  *
                  std::exp (
                            -(
                              std::pow(xi - eta, 2.0)
                              )/(2.0*std::pow(sigma, 2.0))
                            )
                  /////////////////1end
                  +
                  2.0 * (xi - eta)
                  /std::pow(sigma, 4.0)
                  *
                  std::exp (
                            -(
                              std::pow(xi - eta, 2.0)
                              )/(2.0*std::pow(sigma, 2.0))
                            )
                  -
                  std::pow(xi - eta, 3.0)
                  /std::pow(sigma, 6.0)
                  *
                  std::exp (
                            -(
                              std::pow(xi - eta, 2.0)
                              )/(2.0*std::pow(sigma, 2.0))
                            )
                  /////////////2end
                  ;
  if(order == 4)return 1.0
                  /std::pow(sigma, 4.0)
                  *
                  std::exp (
                            -(
                              std::pow(xi - eta, 2.0)
                              )/(2.0*std::pow(sigma, 2.0))
                            )
                  -
                  std::pow(xi - eta,2.0)
                  /std::pow(sigma, 6.0)
                  *
                  std::exp (
                            -(
                              std::pow(xi - eta, 2.0)
                              )/(2.0*std::pow(sigma, 2.0))
                            )
                  +
                  ///////////////1end
                  
                  2.0
                  /std::pow(sigma, 4.0)
                  *
                  std::exp (
                            -(
                              std::pow(xi - eta, 2.0)
                              )/(2.0*std::pow(sigma, 2.0))
                            )
                  -
                  2.0 * std::pow(xi - eta,2.0)
                  /std::pow(sigma, 6.0)
                  *
                  std::exp (
                            -(
                              std::pow(xi - eta, 2.0)
                              )/(2.0*std::pow(sigma, 2.0))
                            )
                  /////////////2end
                  -
                  3.0*std::pow(xi - eta, 2.0)
                  /std::pow(sigma, 6.0)
                  *
                  std::exp (
                            -(
                              std::pow(xi - eta, 2.0)
                              )/(2.0*std::pow(sigma, 2.0))
                            )
                  +
                  std::pow(xi - eta, 4.0)
                  /std::pow(sigma, 8.0)
                  *
                  std::exp (
                            -(
                              std::pow(xi - eta, 2.0)
                              )/(2.0*std::pow(sigma, 2.0))
                            )
                  ////////////////3end
                  ;



  if(order == 5)return  -(xi - eta)
                  /std::pow(sigma, 6.0)
                  *
                  std::exp (
                            -(
                              std::pow(xi - eta, 2.0)
                              )/(2.0*std::pow(sigma, 2.0))
                            )
                  
                  
                  /////////////1end
                  -
                  2.0 * (xi - eta)
                  /std::pow(sigma, 6.0)
                  *
                  std::exp (
                            -(
                              std::pow(xi - eta, 2.0)
                              )/(2.0*std::pow(sigma, 2.0))
                            )
                  
                  +
                  std::pow(xi - eta,3.0)
                  /std::pow(sigma, 8.0)
                  *
                  std::exp (
                            -(
                              std::pow(xi - eta, 2.0)
                              )/(2.0*std::pow(sigma, 2.0))
                            )
                  ///////////////2end
                  -
                  
                  2.0*(xi - eta)
                  /std::pow(sigma, 6.0)
                  *
                  std::exp (
                            -(
                              std::pow(xi - eta, 2.0)
                              )/(2.0*std::pow(sigma, 2.0))
                            )
                  
                  
                  
                  /////////////3end
                  -
                  
                  4.0 * (xi - eta)
                  /std::pow(sigma, 6.0)
                  *
                  std::exp (
                            -(
                              std::pow(xi - eta, 2.0)
                              )/(2.0*std::pow(sigma, 2.0))
                            )
                  
                  +
                  
                  2.0 * std::pow(xi - eta,3.0)
                  /std::pow(sigma, 8.0)
                  *
                  std::exp (
                            -(
                              std::pow(xi - eta, 2.0)
                              )/(2.0*std::pow(sigma, 2.0))
                            )
                  ////////////4end
                  -
                  6.0 * (xi - eta)
                  /std::pow(sigma, 6.0)
                  *
                  std::exp (
                            -(
                              std::pow(xi - eta, 2.0)
                              )/(2.0*std::pow(sigma, 2.0))
                            )
                  +
                  
                  3.0*std::pow(xi - eta, 3.0)
                  /std::pow(sigma, 8.0)
                  *
                  std::exp (
                            -(
                              std::pow(xi - eta, 2.0)
                              )/(2.0*std::pow(sigma, 2.0))
                            )
                  //////////////5end
                  + 
                  4.0*std::pow(xi - eta, 3.0)
                  /std::pow(sigma, 8.0)
                  *
                  std::exp (
                            -(
                              std::pow(xi - eta, 2.0)
                              )/(2.0*std::pow(sigma, 2.0))
                      )
                  -
                  std::pow(xi - eta, 5.0)
                  /std::pow(sigma, 10.0)
                  *
                  std::exp (
                            -(
                              std::pow(xi - eta, 2.0)
                              )/(2.0*std::pow(sigma, 2.0))
                            )
                  ;
}

int current(){
  

  
}
