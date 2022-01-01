#include"PowerSpectrum.h"

//====================================================
// Constructors
//====================================================

PowerSpectrum::PowerSpectrum(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec, 
    Perturbations *pert,
    double A_s,
    double n_s,
    double kpivot_mpc) : 
  cosmo(cosmo), 
  rec(rec), 
  pert(pert),
  A_s(A_s),
  n_s(n_s),
  kpivot_mpc(kpivot_mpc)
{}

//====================================================
// Do all the solving
//====================================================
void PowerSpectrum::solve(){

  //=========================================================================
  // TODO: Choose the range of k's and the resolution to compute Theta_ell(k)
  //=========================================================================
  Vector k_array = Utils::linspace(k_min, k_max, n_k);  
  Vector log_k_array = log(k_array);


  //=========================================================================
  // TODO: Make splines for j_ell. 
  // Implement generate_bessel_function_splines
  //=========================================================================
  generate_bessel_function_splines();
    std::cout << "You are at :  Checkpoint 00!" << "\n";

  // TODO: Line of sight integration to get Theta_ell(k)
  // Implement line_of_sight_integration
  //=========================================================================
  //=========================================================================
  line_of_sight_integration(k_array); // create thetaT_ell and thetaE_ell;
  std::cout << "You are at :  Checkpoint 01!" << "\n";

  //=========================================================================
  // TODO: Integration to get Cell by solving dCell^f/dlogk = Delta(k) * f_ell(k)^2
  // Implement solve_for_cell
  //=========================================================================
  auto cell_TT = solve_for_cell(log_k_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
  
  auto cell_TE = solve_for_cell(log_k_array, thetaT_ell_of_k_spline, thetaE_ell_of_k_spline);
  cell_TE_spline.create(ells, cell_TE, "Cell_TE_of_ell");
  
  auto cell_EE = solve_for_cell(log_k_array, thetaE_ell_of_k_spline, thetaE_ell_of_k_spline);
  cell_EE_spline.create(ells, cell_EE, "Cell_TT_of_ell");
  
  std::cout << "You are at :  Checkpoint 02!" << "\n";
  
}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================

void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("besselspline");
  
  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(ells.size());
    
  //=============================================================================
  // TODO: Compute splines for bessel functions j_ell(z)
  // Choose a suitable range for each ell
  // NB: you don't want to go larger than z ~ 40000, then the bessel routines
  // might break down. Use j_ell(z) = Utils::j_ell(ell, z)
  //=============================================================================
  double n_z = 10001;
  // Vector log_z_array = Utils::linspace(log(1e-8), log(4e4), n_z);
  Vector z_array = Utils::linspace(0.0, 5000, n_z);
  Vector j_ell_z_array(n_z);

  for(size_t i = 0; i < ells.size(); i++){
    const int ell = ells[i];
    
    for(int iz = 0; iz < n_z; iz++)
    {
      j_ell_z_array[iz] = Utils::j_ell(ell, z_array[iz]);
    }

    // Make the j_ell_splines[i] spline
    j_ell_splines[i].create(z_array, j_ell_z_array, "Spline of j_ell(xs)");
  }

  Utils::EndTiming("besselspline");
}

//====================================================
// Do the line of sight integration for a single
// source function
//====================================================

Vector2D PowerSpectrum::line_of_sight_integration_single(
    Vector & k_array, 
    std::function<double(double,double)> &source_function){
  Utils::StartTiming("lineofsight");

  // Make storage for the results
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));
  double eta0 = cosmo->eta_of_x(0.0);
    
  double Theta_ell_ini = 0.0;
  Vector Theta_ell_ic{Theta_ell_ini};
  
  #pragma omp parallel for schedule(dynamic, 4) num_threads(16)
  for(size_t ik = 0; ik < k_array.size(); ik++){

    //=============================================================================
    // TODO: Implement to solve for the general line of sight integral 
    // F_ell(k) = Int dx jell(k(eta-eta0)) * S(x,k) for all the ell values for the 
    // given value of k
    //=============================================================================
    double k = k_array[ik];
    int nx = 500;
    auto x_array = Utils::linspace(x_start, x_end, nx);
    
    for(int iell = 0; iell < ells.size(); iell++){
         
      ODEFunction dTheta_elldx = [&](double x, const double *Theta_ell_k, double *dTheta_elldx){
      double eta = cosmo->eta_of_x(x);
      dTheta_elldx[0] = source_function(x, k)*j_ell_splines[iell](k*(eta0 - eta));          
      return GSL_SUCCESS;
      };

      ODESolver ode;
          // Setting up accuracy parameters
      double hstart = 1e-3, abserr = 1e-10, relerr = 1e-10;
      ode.set_accuracy(hstart, abserr, relerr);
      ode.solve(dTheta_elldx, x_array, Theta_ell_ic, gsl_odeiv2_step_rkf45);

      // result[iell][ik] = ode.get_final_data()[0]; //Theta_ell at the last point, x = 0;
      auto Theta_ell_today = ode.get_data_by_xindex(nx - 1);
      result[iell][ik] = Theta_ell_today[0];
      
    }

    // Store the result for Source_ell(k) in results[ell][ik]
  }

  Utils::EndTiming("lineofsight");
  return result;
}

//====================================================
// Do the line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration(Vector & k_array){
  const int n_k        = k_array.size();
  const int n          = 100;
  const int nells      = ells.size();
  
  // Make storage for the splines we are to create
  thetaT_ell_of_k_spline = std::vector<Spline>(nells);

  //============================================================================
  // TODO: Solve for Theta_ell(k) and spline the result
  //============================================================================

  // Make a function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
  };

  // Do the line of sight integration
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);

  // Spline the result and store it in thetaT_ell_of_k_spline
  for(int iell = 0; iell < ells.size(); iell++){
    thetaT_ell_of_k_spline[iell].create(k_array, thetaT_ell_of_k[iell]);
  }

  //============================================================================
  // TODO: Solve for ThetaE_ell(k) and spline
  //============================================================================
  if(Constants.polarization){
    
  // Make storage for the splines we are to create
  thetaE_ell_of_k_spline = std::vector<Spline>(nells);
  std::function<double(double,double)> source_function_E = [&](double x, double k){
    return pert->get_Source_E(x,k);
  };

  // Do the line of sight integration
  Vector2D thetaE_ell_of_k = line_of_sight_integration_single(k_array, source_function_E);

  // Spline the result and store it in thetaT_ell_of_k_spline
  for(int iell = 0; iell < ells.size(); iell++){
    thetaE_ell_of_k_spline[iell].create(k_array, thetaE_ell_of_k[iell]);
  }

  }
}

//====================================================
// Compute Cell (could be TT or TE or EE) 
// Cell = Int_0^inf 4 * pi * P(k) f_ell g_ell dk/k
//====================================================
Vector PowerSpectrum::solve_for_cell(
    Vector & log_k_array,
    std::vector<Spline> & f_ell_spline,
    std::vector<Spline> & g_ell_spline){
  const int nells      = ells.size();
  const double pi = 3.14159265359;

  //============================================================================
  // TODO: Integrate Cell = Int 4 * pi * P(k) f_ell g_ell dk/k
  // or equivalently solve the ODE system dCell/dlogk = 4 * pi * P(k) * f_ell * g_ell
  //============================================================================

  Vector result(nells);
  for(int iell = 0; iell < nells; iell++){

    ODEFunction dC_elldlogk = [&](double k, const double *C_ell, double *dC_elldlogk){
    dC_elldlogk[0] = 4*pi*A_s*pow(Constants.Mpc*exp(k)/kpivot_mpc, n_s-1.0)*f_ell_spline[iell](exp(k))*g_ell_spline[iell](exp(k));          
    return GSL_SUCCESS;
    };

    double C_ell_ini = 0.0;
    Vector C_ell_ic{C_ell_ini};

    ODESolver ode;
    ode.solve(dC_elldlogk, log_k_array, C_ell_ic);

    // result[iell]= ode.get_final_data()[0]; //C_ell at the last point, k_max;
    
    auto Cell = ode.get_data_by_xindex(log_k_array.size() - 1);
    // Saving result
    result[iell] = Cell[0];

  }

  return result;
}

//====================================================
// The primordial power-spectrum
//====================================================

double PowerSpectrum::primordial_power_spectrum(const double k) const{
  return A_s * pow( Constants.Mpc * k / kpivot_mpc , n_s - 1.0);
}

//====================================================
// P(k) in units of (Mpc)^3
//====================================================

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k_mpc) const{

  //=============================================================================
  // TODO: Compute the matter power spectrum
  //=============================================================================
  double OmegaCDM0 = cosmo -> get_OmegaCDM(0);
  double OmegaB0   = cosmo -> get_OmegaB(0);
  double OmegaM0   = OmegaB0 + OmegaCDM0;
  double H0        = cosmo -> get_H0();
  
  double Phi = pert -> get_Phi(x, k_mpc); // Metric potential perturbation 
  double c = Constants.c;         
  
  double DeltaM = 2.0 * c * c * k_mpc * k_mpc * Phi / (3.0 * OmegaM0 * exp(-x) * H0 * H0);
  double P_primordial = primordial_power_spectrum(k_mpc);
  double pofk = std::fabs(DeltaM) * std::fabs(DeltaM) * P_primordial;

  return pofk;
}

//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_cell_TT(const double ell) const{
  return cell_TT_spline(ell);
}
double PowerSpectrum::get_cell_TE(const double ell) const{
  return cell_TE_spline(ell);
}
double PowerSpectrum::get_cell_EE(const double ell) const{
  return cell_EE_spline(ell);
}

//====================================================
// Output the cells to file
//====================================================

void PowerSpectrum::output(std::string filename) const{
  // Output in standard units of muK^2
  std::ofstream fp(filename.c_str());
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);
  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    double normfactorN = (ell * (ell+1)) / (2.0 * M_PI) 
      * pow(1e6 * cosmo->get_TCMB() *  pow(4.0/11.0, 1.0/3.0), 2);
    double normfactorL = (ell * (ell+1)) * (ell * (ell+1)) / (2.0 * M_PI);
    fp << ell                                 << " ";
    fp << cell_TT_spline( ell ) * normfactor  << " ";
    if(Constants.polarization){
      fp << cell_EE_spline( ell ) * normfactor  << " ";
      fp << cell_TE_spline( ell ) * normfactor  << " ";
    }
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
}

