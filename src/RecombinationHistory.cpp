#include"RecombinationHistory.h"

//====================================================
// Constructors
//====================================================
   
RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo, 
    bool Reionization,
    double z_reion,
    double delta_z_reion,
    double z_Hereion,
    double delta_z_Hereion,
    double Yp) :
  cosmo(cosmo),
  Reionization(Reionization), 
  z_reion(z_reion),
  delta_z_reion(delta_z_reion),
  z_Hereion(z_Hereion),
  delta_z_Hereion(delta_z_Hereion),
  Yp(Yp)
{}

//====================================================
// Do all the solving we need to do
//====================================================

void RecombinationHistory::solve(){
    
  // Compute and spline Xe, ne
  solve_number_density_electrons();
   
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");
  
  //=============================================================================
  // TODO: Set up x-array and make arrays to store X_e(x) and n_e(x) on
  //=============================================================================
  const double G           = Constants.G;
  const double m_H         = Constants.m_H;
  const double pi          = 3.14159265358;  
  
  double OmegaB0      = cosmo->get_OmegaB(0.0);
  double H0           = cosmo->get_H0();
  double rho_c        = 3.0*pow(H0, 2)/8.0/pi/G; 
  
  Vector x_array = Utils::linspace(Constants.x_start, Constants.x_end, npts_rec_arrays); // npt defined in the header

  Vector Xe_arr(npts_rec_arrays);
  Vector ne_arr(npts_rec_arrays);

  // Calculate recombination history
  bool saha_regime = true;
  for(int i = 0; i < npts_rec_arrays; i++){

    //==============================================================
    // TODO: Get X_e from solving the Saha equation so
    // implement the function electron_fraction_from_saha_equation
    //==============================================================
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit)
      saha_regime = false;

    if(saha_regime){
      
      //=============================================================================
      // TODO: Store the result we got from the Saha equation
      //=============================================================================
      //...
      //...
      Xe_arr[i] = Xe_ne_data.first;
      ne_arr[i] = Xe_ne_data.second;

    } else {

      //==============================================================
      // TODO: Compute X_e from current time til today by solving 
      // the Peebles equation (NB: if you solve all in one go remember to
      // exit the for-loop!)
      // Implement rhs_peebles_ode
      //==============================================================
      //...
      //...
      
      
      // The Peebles ODE equation
      ODESolver peebles_Xe_ode;
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        return rhs_peebles_ode(x, Xe, dXedx);
      };
      
      //=============================================================================
      // TODO: Set up IC, solve the ODE and fetch the result 
      //=============================================================================
      //...
      //...
      
      double Xeini = Xe_arr[i-1];
      Vector Xe_ic{Xeini};

      Vector x_array_Peebles;
      x_array_Peebles = {x_array[i-1],x_array[i]}; //solve the ode at present and previous x;
      
      // Solve the ODE
      peebles_Xe_ode.solve(dXedx, x_array_Peebles, Xe_ic);
      auto Xe_i = peebles_Xe_ode.get_data_by_component(0);
      
      if (Reionization)
      {
        double fHe = Yp/(4.0*(1.0-Yp));
        double y = exp(-3.0*x_array[i]/2.0);
        double y_reion = pow(1.0 + z_reion, 1.5);
        double delta_y_reion = 3.0/2.0*sqrt(1.0 + z_reion)*delta_z_reion;
        
        double new1 = (1 + fHe)/2.0*(1 + tanh((y_reion - y)/delta_y_reion));
        double new2 = fHe/2.0*(1.0 + tanh((z_Hereion -(1/exp(x_array[i]) - 1))/delta_z_Hereion));
        Xe_arr[i] = Xe_i[1] + new1 + new2;
      }
      else
      {
        Xe_arr[i] = Xe_i[1];
      }
      
        
      const double a           = exp(x_array[i]);
      double nb = OmegaB0*rho_c/m_H/pow(a, 3);
      double nH = (1 - Yp)*nb;
      
      ne_arr[i] = Xe_arr[i]*nH;
    
    }
  }

  //=============================================================================
  // TODO: Spline the result. Implement and make sure the Xe_of_x, ne_of_x 
  // functions are working
  //=============================================================================
  //...
  //...
  Vector log_Xe_arr(npts_rec_arrays);
  Vector log_ne_arr(npts_rec_arrays);
  
  for(int i = 0; i < npts_rec_arrays; i++){
    log_Xe_arr[i] = log(Xe_arr[i]);
    log_ne_arr[i] =  log(ne_arr[i]);
  }

  log_Xe_of_x_spline.create(x_array, log_Xe_arr,"Xe"); 
  log_ne_of_x_spline.create(x_array, log_ne_arr, "ne");

  Utils::EndTiming("Xe");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  const double a           = exp(x);
 
  // Physical constants
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double epsilon_0   = Constants.epsilon_0;
  const double H0_over_h   = Constants.H0_over_h;
  const double pi          = 3.14159265358;
  
  const double xhi0        = Constants.xhi0;              // Ionization energy for neutral Helium
  const double xhi1        = Constants.xhi1;             // Ionization energy for singly ionized Helium
  // Fetch cosmological parameters
  
  double OmegaB0      = cosmo->get_OmegaB(0.0);
  double H0           = cosmo->get_H0();
  double rho_c        = 3.0*pow(H0, 2)/8.0/pi/G;
  
  double TCMB_a         = cosmo->get_TCMB(x);
  double Tb_a           = TCMB_a;

  // Electron fraction and number density
  double Xe = 0.0;
  double ne = 0.0;
  
  //=============================================================================
  // TODO: Compute Xe and ne from the Saha equation
  //=============================================================================
  //...
  //...
  double nb = OmegaB0*rho_c/m_H/pow(a, 3);
  double nH = (1 - Yp)*nb;
  
  // rhs of 3 equations for xHe+, xHe++, xH+
  double r1 = 2*(pow(m_e*pow(hbar, -2)*k_b*Tb_a/2/pi, 3.0/2.0))*exp(-1*xhi0/Tb_a/k_b);
  double r2 = 4*(pow(m_e*pow(hbar, -2)*k_b*Tb_a/2/pi, 3.0/2.0))*exp(-1*xhi1/Tb_a/k_b);
  double r3 = 1*(pow(m_e*pow(hbar, -2)*k_b*Tb_a/2/pi, 3.0/2.0))*exp(-1*epsilon_0/Tb_a/k_b);
    
  double fe = 1;
  double fe_old = 0;  
  double xHe2, xHe1, xH1; 
  
   // iteratively solve the Saha equation with helium considered
  // while(abs(fe - fe_old) > 10e-10)
  int i;
  while(i < 5)

  {
  fe_old = fe;
    
  ne = fe*nb;
  
  xHe1 = r1/(ne + r1 + r1*r2/ne);
  xHe2 = r2*xHe1/ne;
  xH1  = r3/(ne + r3);
    
  fe = (2*xHe2 + xHe1)*Yp/4.0 + xH1*(1-Yp);
        // std::cout << "fe:          " << abs(fe - fe_old)          << "\n";
  i++; 
  }
  
  ne = fe*nb;
  Xe = ne/nH;

  return std::pair<double,double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of a and X_e
  const double X_e         = Xe[0];
  const double a           = exp(x);

  // Physical constants in SI units
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double c           = Constants.c;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double sigma_T     = Constants.sigma_T;
  const double lambda_2s1s = Constants.lambda_2s1s;
  const double epsilon_0   = Constants.epsilon_0;
  const double pi          = 3.14159265358;

  // Cosmological parameters
  // const double OmegaB      = cosmo->get_OmegaB();
  double OmegaB0      = cosmo->get_OmegaB(0.0);
  double H0           = cosmo->get_H0();
  double H            = cosmo->H_of_x(x);
  double rho_c        = 3.0*pow(H0, 2)/8.0/pi/G;
  
  double TCMB_a         = cosmo->get_TCMB(x);
  double Tb_a           = TCMB_a;

  //=============================================================================
  // TODO: Write the expression for dXedx
  //=============================================================================
  double nb = OmegaB0*rho_c/m_H/pow(a, 3);
  double nH = (1 - Yp)*nb;
  
  double phi2 = 0.448*log(epsilon_0/Tb_a/k_b);
  double alpha = 1.0/137.0359992;
  double alpha2 = 64*pi/sqrt(27*pi)*pow(alpha, 2)*pow(m_e, -2)*pow(hbar, 2)/c*sqrt(epsilon_0/Tb_a/k_b)*phi2;
  
  double beta = alpha2*(pow(m_e*pow(hbar, -2)*k_b*Tb_a/2/pi, 3.0/2.0))*exp(-1*epsilon_0/Tb_a/k_b);
  double beta2;
  
  if(3.0*epsilon_0/(4.0*k_b*Tb_a) > 200){
     beta2 = 0.0;
  } 
  else {
  beta2 = beta*exp(3.0*epsilon_0/(4.0*k_b*Tb_a));
  }
  
  double n1s = (1 - X_e)*nH;
  double L_alpha = H*pow(3*epsilon_0, 3)/pow(8*pi, 2)/n1s/pow(hbar, 3)/pow(c, 3);
  double Cr = (lambda_2s1s + L_alpha)/(lambda_2s1s + L_alpha + beta2);    
  
  dXedx[0] = Cr/H*(beta*(1-X_e) - nH*alpha2*pow(X_e, 2));

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");
  const double sigma_T     = Constants.sigma_T;
  const double c           = Constants.c;

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation
  const int npts = 1000;
  Vector x_array = Utils::linspace(x_start, x_end, npts);
  Vector g_tilde_array(npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){

    //=============================================================================
    // TODO: Write the expression for dtaudx
    //=============================================================================
    //...
    //...
    double H            = cosmo->H_of_x(x);
    double ne           = ne_of_x(x);

    // Set the derivative for photon optical depth
    dtaudx[0] = -1*c*ne*sigma_T/H;

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set up and solve the ODE and make tau splines
  //=============================================================================
    
  double tauini = 1.0;
  Vector tau_ic{tauini};
  
  // Solve the ODE
  ODESolver ode;
  ode.solve(dtaudx, x_array, tau_ic);
  auto tau_array = ode.get_data_by_component(0);
  
  for (int i = 0; i< npts; i++)
  {
  tau_array[i] = tau_array[i] - tau_array[npts-1];
  }
  tau_of_x_spline.create(x_array, tau_array, "tau");

  //=============================================================================
  // TODO: Compute visibility functions and spline everything
  //=============================================================================
  
  for (int i = 0; i< npts; i++)
  {
    g_tilde_array[i] = -1*dtaudx_of_x(x_array[i])*exp(-1*tau_array[i]);
  }

  g_tilde_of_x_spline.create(x_array, g_tilde_array, "g");
  
  Utils::EndTiming("opticaldepth");
}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement. Either from the tau-spline tau_of_x_spline.deriv_(x) or 
  // from a separate spline if you choose to do this
  //=============================================================================
  //...
  //...

  return tau_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return tau_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return g_tilde_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return g_tilde_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::Xe_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return exp(log_Xe_of_x_spline(x));
}

double RecombinationHistory::ne_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return exp(log_ne_of_x_spline(x));
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";
  // auto xe = electron_fraction_from_saha_equation(-9.12736);
  // std::cout << "xe:"            << xe.first<< "\n";
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = 1000;
  const double x_min   = x_start;
  const double x_max   = x_end;
  
  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                    << " ";
    fp << Xe_of_x(x)           << " ";
    fp << electron_fraction_from_saha_equation(-9.12736).first <<" ";
    fp << ne_of_x(x)           << " ";
    fp << tau_of_x(x)          << " ";
    fp << dtaudx_of_x(x)       << " ";
    fp << ddtauddx_of_x(x)     << " ";
    fp << g_tilde_of_x(x)      << " ";
    fp << dgdx_tilde_of_x(x)   << " ";
    fp << ddgddx_tilde_of_x(x) << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

