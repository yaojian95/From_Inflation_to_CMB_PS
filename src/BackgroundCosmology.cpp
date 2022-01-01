#include "BackgroundCosmology.h"

//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double OmegaB, 
    double OmegaCDM, 
    double OmegaK,
    double Neff, 
    double TCMB) :
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
  OmegaK(OmegaK),
  Neff(Neff), 
  TCMB(TCMB)
{

  //=============================================================================
  // TODO: Compute OmegaR, OmegaNu, OmegaLambda, H0, ...
  //=============================================================================

    H0 = Constants.H0_over_h*h;
    double pi = 3.14159265358;
    OmegaR = 2*pow(pi, 2)/30*pow(Constants.k_b*TCMB,4)/pow(Constants.hbar,3)/pow(Constants.c,5)*8*pi*Constants.G/3/pow(H0,2);
    OmegaNu = Neff*7.0/8.0*pow(4.0/11.0, 4.0/3.0)*OmegaR;
    OmegaLambda = 1 - OmegaB - OmegaK - OmegaCDM - OmegaR - OmegaNu;
    
}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(){
  Utils::StartTiming("Eta");
    
  //=============================================================================
  // TODO: Set the range of x and the number of points for the splines
  // For this Utils::linspace(x_start, x_end, npts) is useful
  //=============================================================================
  
  // const double x_min = -15.0;
  // const double x_max =  1.0;
  // const int    n_pts =  200;
  
  Vector x_array = Utils::linspace(x_start, x_end, 1000);

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){

    //=============================================================================
    // TODO: Set the rhs of the detadx ODE
    //=============================================================================

    detadx[0] = Constants.c/Hp_of_x(x);

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set the initial condition, set up the ODE system, solve and make
  // the spline eta_of_x_spline 
  //=============================================================================
  
  double etaini = Constants.c/Hp_of_x(x_array[0]);
  Vector eta_ic{etaini};
  
  // Solve the ODE
  ODESolver ode;
  ode.solve(detadx, x_array, eta_ic);
  auto eta_array = ode.get_data_by_component(0);
  
  // Spline eta_of_x_spline; //already defined in .h file
 
  eta_of_x_spline.create(x_array,eta_array,"eta_of_x"); 

  Utils::EndTiming("Eta");
}

//====================================================
// Get methods
//====================================================

double BackgroundCosmology::H_of_x(double x) const{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...
  double a = exp(x);
  double H;
  
  H = H0*sqrt((OmegaB + OmegaCDM)*pow(a, -3) + (OmegaR + OmegaNu)*pow(a, -4) + OmegaK*pow(a, -2) + OmegaLambda);

    return H;
}

double BackgroundCosmology::Hp_of_x(double x) const{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...
  double a = exp(x);
  return H_of_x(x)*a;
}

double BackgroundCosmology::dHpdx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  double a = exp(x);
  double f_a = (OmegaB + OmegaCDM)*pow(a, -3) + (OmegaR + OmegaNu)*pow(a, -4) + OmegaK*pow(a, -2) + OmegaLambda;
  double df_a = -3.0*(OmegaB + OmegaCDM)*pow(a, -4) - 4.0*(OmegaR + OmegaNu)*pow(a, -5) - 2.0*OmegaK*pow(a, -3);
  
  double dHp_a = H_of_x(x) + a*H0*0.5*pow(f_a, -0.5)*df_a;

  return dHp_a*a;
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  double a = exp(x);
  double under_root = (OmegaB+OmegaCDM)*pow(a,-1) + (OmegaR+OmegaNu)*pow(a,-2) + OmegaLambda*pow(a, 2);
  double dunder_root_dx = -(OmegaB+OmegaCDM)*pow(a,-1) -2.0*(OmegaR+OmegaNu)*pow(a,-2) + 2.0*OmegaLambda*pow(a, 2);
  double ddunder_root_dx2 = (OmegaB+OmegaCDM)*pow(a,-1) +4.0*(OmegaR+OmegaNu)*pow(a,-2) + 4.0*OmegaLambda*pow(a, 2);
  double ddHpddx = H0*0.5*(-0.5)*pow(under_root, -1.5)*pow(dunder_root_dx,2) + H0*0.5*pow(under_root, -0.5)*ddunder_root_dx2;
  return ddHpddx;

}

double BackgroundCosmology::get_OmegaB(double x) const{ 
  if(x == 0.0) return OmegaB;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...
  double a = exp(x);
  double OmegaBB = OmegaB/(pow(a, 3)*pow(H_of_x(x), 2)/pow(H0, 2));
  return OmegaBB;
}

double BackgroundCosmology::get_OmegaR(double x) const{ 
  if(x == 0.0) return OmegaR;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================

  double a = exp(x);
  double OmegaRR = OmegaR/(pow(a, 4)*pow(H_of_x(x), 2)/pow(H0, 2));

  return OmegaRR;
}

double BackgroundCosmology::get_OmegaNu(double x) const{ 
  if(x == 0.0) return OmegaNu;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  double a = exp(x);
  double OmegaNuNu = OmegaNu/(pow(a, 4)*pow(H_of_x(x), 2)/pow(H0, 2));

  return OmegaNuNu;
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  if(x == 0.0) return OmegaCDM;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  double a = exp(x);
  double OmegaCDMCDM = OmegaCDM/(pow(a, 3)*pow(H_of_x(x), 2)/pow(H0, 2));

  return OmegaCDMCDM;
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  if(x == 0.0) return OmegaLambda;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  double a = exp(x);
  double OmegaLambdaL = OmegaLambda/(pow(H_of_x(x), 2)/pow(H0, 2));

  return OmegaLambdaL;
}

double BackgroundCosmology::get_OmegaK(double x) const{ 
  if(x == 0.0) return OmegaK;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  double a = exp(x);
  double OmegaKK = OmegaK/(pow(a, 2)*pow(H_of_x(x), 2)*pow(H0, 2));

  return OmegaKK;
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const{
  return H0; 
}

double BackgroundCosmology::get_h() const{ 
  return h; 
}

double BackgroundCosmology::get_Neff() const{ 
  return Neff; 
}

double BackgroundCosmology::get_TCMB(double x) const{ 
  if(x == 0.0) return TCMB;
  return TCMB * exp(-x); 
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:      " << OmegaB      << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM    << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda << "\n";
  std::cout << "OmegaK:      " << OmegaK      << "\n";
  std::cout << "OmegaNu:     " << OmegaNu     << "\n";
  std::cout << "OmegaR:      " << OmegaR      << "\n";
  std::cout << "Neff:        " << Neff        << "\n";
  std::cout << "h:           " << h           << "\n";
  std::cout << "H0:           " << H0           << "\n";
  std::cout << "TCMB:        " << TCMB        << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = -15.0;
  const double x_max =  1.0;
  const int    n_pts =  200;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << exp(x)              << " ";
    
    fp << H_of_x(x) /H0        << " ";
    
    // fp << x                  << " ";
    fp << eta_of_x(x)/Constants.Mpc << " ";
    // fp << Hp_of_x(x)         << " ";
    fp << get_OmegaB(x)      << " ";
    fp << get_OmegaCDM(x)    << " ";
    fp << get_OmegaLambda(x) << " ";
    fp << get_OmegaR(x)      << " ";
    fp << get_OmegaNu(x)     << " ";
    fp << get_OmegaK(x)      << " ";
    fp << dHpdx_of_x(x)      << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

