#include"Perturbations.h"

//====================================================
// Constructors
//====================================================

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec) : 
  cosmo(cosmo), 
  rec(rec)
{}

//====================================================
// Do all the solving
//====================================================

void Perturbations::solve(){

  // Integrate all the perturbation equation and spline the result
  integrate_perturbations();

  // Compute source functions and spline the result
  compute_source_functions();
}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================

void Perturbations::integrate_perturbations(){
  Utils::StartTiming("integrateperturbation");

  //===================================================================
  // TODO: Set up the k-array for the k's we are going to integrate over
  // Start at k_min end at k_max with n_k points with either a
  // quadratic or a logarithmic spacing
  //===================================================================
  Vector k_array(n_k); //n_k defined in the header
  Vector log_k_array = Utils::linspace(log(k_min), log(k_max), n_k);
  for(int i=0; i<n_k; i++){
    k_array[i] = exp(log_k_array[i]);}
  
  Vector x_array = Utils::linspace(x_start, x_end, n_x);
  
  Vector delta_cdm_x_k(n_x * n_k);  // n_x rows and n_k columns
  Vector delta_b_x_k(n_x * n_k); 
  Vector v_cdm_x_k(n_x * n_k);
  Vector v_b_x_k(n_x * n_k);
  Vector Phi_x_k(n_x * n_k);
  Vector Theta0_x_k(n_x * n_k);
  Vector Theta1_x_k(n_x * n_k);  
  Vector Nu2_x_k(n_x * n_k);  
  
  Vector Psi_x_k(n_x * n_k);
  Vector Pi_x_k(n_x * n_k, 0);
  
  const double c       = Constants.c;
  const double H0           = cosmo->get_H0();
  const double OmegaR0  = cosmo->get_OmegaR(0.0);
  const double OmegaNu0  = cosmo->get_OmegaNu(0.0);
  
  // Loop over all wavenumbers
  #pragma omp parallel for schedule(dynamic, 2) num_threads(16)
  // #pragma omp parallel for schedule(static) num_threads(16)
  for(int ik = 0; ik < n_k; ik++){

    // Progress bar...
    if( (10*ik) / n_k != (10*ik+10) / n_k ) {
      std::cout << (100*ik+100)/n_k << "% " << std::flush;
      if(ik == n_k-1) std::cout << std::endl;
    }

    // Current value of k
    double k = k_array[ik];

    // Find value to integrate to, as a function k
    double x_end_tight = get_tight_coupling_time(k);
    // std::cout << "x_tight_end "<<x_end_tight <<std::endl;
    int idx_end_tight = (int) ((x_end_tight - x_start)/(x_end - x_start) * (n_x - 1.0));
    // std::cout << "x_array_tight_end "<<x_array[idx_end_tight] <<std::endl;

    //===================================================================
    // TODO: Tight coupling integration
    // Remember to implement the routines:
    // set_ic : The IC at the start
    // rhs_tight_coupling_ode : The dydx for our coupled ODE system
    //===================================================================

    // Set up initial conditions in the tight coupling regime
    auto y_tight_coupling_ini = set_ic(x_start, k);

    // The tight coupling ODE system
    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };
 
    // Integrate from x_start -> x_end_tight
    Vector x_tc = Utils::linspace(x_start, x_array[idx_end_tight], idx_end_tight+1);
    
    ODESolver tc_ode;
    tc_ode.solve(dydx_tight_coupling, x_tc, y_tight_coupling_ini);
    
    auto d_cdm_tc_arr = tc_ode.get_data_by_component(0);
    auto d_b_tc_arr = tc_ode.get_data_by_component(1);
    auto v_cdm_tc_arr = tc_ode.get_data_by_component(2);
    auto v_b_tc_arr = tc_ode.get_data_by_component(3);
    auto Phi_tc_arr = tc_ode.get_data_by_component(4);
    auto Theta0_tc_arr = tc_ode.get_data_by_component(5);
    auto Theta1_tc_arr = tc_ode.get_data_by_component(6);
    auto Nu2_tc_arr = tc_ode.get_data_by_component(9);

    // std::cout << "Theta1 at x_end_tight" << Theta1_tc_arr.back() <<"\n";
    //====i===============================================================
    // TODO: Full equation integration
    // Remember to implement the routines:
    // set_ic_after_tight_coupling : The IC after tight coupling ends
    // rhs_full_ode : The dydx for our coupled ODE system
    //===================================================================

    // Set up initial conditions (y_tight_coupling is the solution at the end of tight coupling)
    auto y_tight_coupling = tc_ode.get_final_data();
    auto y_full_ini = set_ic_after_tight_coupling(y_tight_coupling, x_end_tight, k);

    // The full ODE system
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };

    // Integrate from x_end_tight -> x_end
    Vector x_full = Utils::linspace(x_array[idx_end_tight], x_end, n_x - idx_end_tight); //in total nx+1 points with x_end_tight twice;

    ODESolver full_ode;
    full_ode.solve(dydx_full, x_full, y_full_ini);
    
    auto d_cdm_full_arr = full_ode.get_data_by_component(0);
    auto d_b_full_arr = full_ode.get_data_by_component(1);
    auto v_cdm_full_arr = full_ode.get_data_by_component(2);
    auto v_b_full_arr = full_ode.get_data_by_component(3);
    auto Phi_full_arr = full_ode.get_data_by_component(4);
    auto Theta0_full_arr = full_ode.get_data_by_component(5);
    auto Theta1_full_arr = full_ode.get_data_by_component(6);
    auto Theta2_full_arr = full_ode.get_data_by_component(7);
    auto Thetap0_full_arr = full_ode.get_data_by_component(13);
    auto Thetap2_full_arr = full_ode.get_data_by_component(14);
    auto Nu2_full_arr = full_ode.get_data_by_component(23);
    
    // x_full.erase(x_full.begin());
    d_cdm_full_arr.erase(d_cdm_full_arr.begin());
    d_b_full_arr.erase(d_b_full_arr.begin());
    v_cdm_full_arr.erase(v_cdm_full_arr.begin());
    v_b_full_arr.erase(v_b_full_arr.begin());
    Phi_full_arr.erase(Phi_full_arr.begin());
    Theta0_full_arr.erase(Theta0_full_arr.begin());
    Theta1_full_arr.erase(Theta1_full_arr.begin());
    Theta2_full_arr.erase(Theta2_full_arr.begin());
    Thetap0_full_arr.erase(Thetap0_full_arr.begin());   
    Thetap2_full_arr.erase(Thetap2_full_arr.begin());   
    // std::cout << "Theta1 at x_end_tight + 0.00001" << Theta1_full_arr.front() <<"\n";
    // combine x_tc and x_full, d_tc and d_full;
    Vector x_all = x_tc;
    Vector d_cdm_all = d_cdm_tc_arr;
    Vector d_b_all = d_b_tc_arr;
    Vector v_cdm_all = v_cdm_tc_arr;
    Vector v_b_all = v_b_tc_arr;
    Vector Phi_all = Phi_tc_arr;
    Vector Theta0_all = Theta0_tc_arr;
    Vector Theta1_all = Theta1_tc_arr;
    Vector Nu2_all = Nu2_tc_arr;
    
    x_all.insert( x_all.end(), x_full.begin(), x_full.end());
    d_cdm_all.insert( d_cdm_all.end(), d_cdm_full_arr.begin(), d_cdm_full_arr.end());
    d_b_all.insert( d_b_all.end(), d_b_full_arr.begin(), d_b_full_arr.end());
    v_cdm_all.insert( v_cdm_all.end(), v_cdm_full_arr.begin(), v_cdm_full_arr.end());
    v_b_all.insert( v_b_all.end(), v_b_full_arr.begin(), v_b_full_arr.end());
    Phi_all.insert( Phi_all.end(), Phi_full_arr.begin(), Phi_full_arr.end());
    Theta0_all.insert(Theta0_all.end(), Theta0_full_arr.begin(), Theta0_full_arr.end());
    Theta1_all.insert(Theta1_all.end(), Theta1_full_arr.begin(), Theta1_full_arr.end());
    Nu2_all.insert(Nu2_all.end(), Nu2_full_arr.begin(), Nu2_full_arr.end());

    //===================================================================
    // TODO: remember to store the data found from integrating so we can
    // spline it below
    //
    // To compute a 2D spline of a function f(x,k) the data must be given 
    // to the spline routine as a 1D array f_array with the points f(ix, ik) 
    // stored as f_array[ix + n_x * ik]
    // Example:

    for(int ix = 0; ix<n_x; ix++)
    {
      double a = exp(x_array[ix]);
      delta_cdm_x_k[ix + n_x*ik] = d_cdm_all[ix];
      delta_b_x_k[ix + n_x*ik] = d_b_all[ix];
      v_cdm_x_k[ix + n_x*ik] = v_cdm_all[ix];
      v_b_x_k[ix + n_x*ik] = v_b_all[ix];
      Phi_x_k[ix + n_x*ik] = Phi_all[ix];
      Theta0_x_k[ix + n_x*ik] = Theta0_all[ix];
      Theta1_x_k[ix + n_x*ik] = Theta1_all[ix];
      Nu2_x_k[ix + n_x*ik] = Nu2_all[ix];
      
      if(ix <= idx_end_tight)
      {
        Psi_x_k[ix + n_x*ik]  = -1.0*Phi_all[ix];
      }
      else
      {
        // std::cout<<"ix: " << ix <<"idx_end_tight: "<< idx_end_tight;
        Pi_x_k[ix + n_x*ik] = Theta2_full_arr[ix - idx_end_tight -1] + Thetap0_full_arr[ix - idx_end_tight -1] + Thetap2_full_arr[ix - idx_end_tight -1];
        Psi_x_k[ix + n_x*ik]  = -1.0*Phi_all[ix] - 12.0*H0*H0/pow(c*k*a, 2)*(OmegaR0*Theta2_full_arr[ix - idx_end_tight -1] + OmegaNu0*Nu2_all[ix]);
      }
    }
    
    // We can now use the spline as f_spline(x, k)
    //
    // NB: If you use Theta_spline then you have to allocate it first,
    // before using it e.g.
    // Theta_spline = std::vector<Spline2D>(n_ell_theta);

  }
  const int n_ell_theta         = Constants.n_ell_theta;
  // const int n_ell_thetap        = Constants.n_ell_thetap;
  // const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  Theta_spline = std::vector<Spline2D>(n_ell_theta);
  // Theta_p_spline = std::vector<Spline2D>(n_ell_thetap)
  // Nu_spline   = std::vector<Spline2D>(n_ell_neutrinos)
        
  delta_cdm_spline.create(x_array, k_array, delta_cdm_x_k, "delt_cdm__x_k");
  delta_b_spline.create(x_array, k_array, delta_b_x_k);
  v_cdm_spline.create(x_array, k_array, v_cdm_x_k);
  v_b_spline.create(x_array, k_array, v_b_x_k);
  Phi_spline.create(x_array, k_array, Phi_x_k);
  Theta_spline[0].create(x_array, k_array, Theta0_x_k);                    
  Theta_spline[1].create(x_array, k_array, Theta1_x_k); 
  Psi_spline.create(x_array, k_array, Psi_x_k);
  Pi_spline.create(x_array, k_array, Pi_x_k);
                      
  Utils::EndTiming("integrateperturbation");

  //=============================================================================
  // TODO: Make all splines needed: Theta0,Theta1,Theta2,Phi,Psi,...
  //=============================================================================
  // ...
  // ...
  // ...
}

//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{
  
  // std::cout << "You are at :  Checkpoint 00!" << "\n";

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);

  //=============================================================================
  // Compute where in the y_tc array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const int n_ell_tot_tc        = Constants.n_ell_tot_tc;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // References to the tight coupling quantities
  double &delta_cdm    =  y_tc[Constants.ind_deltacdm_tc];
  double &delta_b      =  y_tc[Constants.ind_deltab_tc];
  double &v_cdm        =  y_tc[Constants.ind_vcdm_tc];
  double &v_b          =  y_tc[Constants.ind_vb_tc];
  double &Phi          =  y_tc[Constants.ind_Phi_tc];
  double *Theta        = &y_tc[Constants.ind_start_theta_tc];
  double *Nu           = &y_tc[Constants.ind_start_nu_tc];
  
  const double c       = Constants.c;
  //=============================================================================
  // TODO: Set the initial conditions in the tight coupling regime
  //=============================================================================

  const double H0           = cosmo->get_H0();
  const double OmegaR0  = cosmo->get_OmegaR(0.0);
  const double OmegaNu0 = cosmo->get_OmegaNu(0.0);
  const double Hp       = cosmo->Hp_of_x(x);
  const double fNu      = OmegaNu0/(OmegaR0 + OmegaNu0);
  const double Psi      = -1.0/(3.0/2.0 + 2.0*fNu/5.0);

  // SET: Scalar quantities (Gravitational potential, baryons and CDM)

  Phi = -1*(1.0 + 2.0*fNu/5.0)*Psi;
  delta_cdm = -1*3.0/2.0*Psi;
  delta_b = -1*3.0/2.0*Psi;
  v_cdm = -1*c*k/2/Hp*Psi; 
  v_b = -1*c*k/2/Hp*Psi;

  // SET: Photon temperature perturbations (Theta_ell)

  Theta[0] = -1.0/2.0*Psi;
  Theta[1] = c*k/6.0/Hp*Psi;

  // SET: Neutrino perturbations (N_ell)
  if(neutrinos){
    // n_ell_neutrinos_tc = 8
    Nu[0] = -0.5*Psi;
    Nu[1] = c*k/6.0/Hp*Psi;
    Nu[2] = -1.0*pow(c*k*exp(x), 2)*(Phi + Psi)/12.0/pow(H0, 2)/OmegaNu0;
    
    for(int ell = 3; ell < 8; ell++)
    {
      Nu[ell] = c*k/(2.0*ell +1.0)/Hp*Nu[ell - 1];
    }
  }
  return y_tc;
}

//====================================================
// Set IC for the full ODE system after tight coupling 
// regime ends
//====================================================

Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x, 
    const double k) const{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);
  
  //=============================================================================
  // Compute where in the y array each component belongs and where corresponding
  // components are located in the y_tc array
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Number of multipoles we have in the full regime
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;
  
  const double c       = Constants.c;
  double Hp         = cosmo->Hp_of_x(x);
  double eta        = cosmo->eta_of_x(x);
  double dtau       = rec->dtaudx_of_x(x);
  

  // Number of multipoles we have in the tight coupling regime
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;

  // References to the tight coupling quantities
  const double &delta_cdm_tc    =  y_tc[Constants.ind_deltacdm_tc];
  const double &delta_b_tc      =  y_tc[Constants.ind_deltab_tc];
  const double &v_cdm_tc        =  y_tc[Constants.ind_vcdm_tc];
  const double &v_b_tc          =  y_tc[Constants.ind_vb_tc];
  const double &Phi_tc          =  y_tc[Constants.ind_Phi_tc];
  const double *Theta_tc        = &y_tc[Constants.ind_start_theta_tc];
  const double *Nu_tc           = &y_tc[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set
  double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  double &delta_b         =  y[Constants.ind_deltab_tc];
  double &v_cdm           =  y[Constants.ind_vcdm_tc];
  double &v_b             =  y[Constants.ind_vb_tc];
  double &Phi             =  y[Constants.ind_Phi_tc];
  double *Theta           = &y[Constants.ind_start_theta_tc];
  double *Theta_p         = &y[Constants.ind_start_thetap_tc];
  double *Nu              = &y[Constants.ind_start_nu_tc];

  //=============================================================================
  // TODO: fill in the initial conditions for the full equation system below
  // NB: remember that we have different number of multipoles in the two
  // regimes so be careful when assigning from the tc array
  //=============================================================================

  // SET: Scalar quantities (Gravitational potental, baryons and CDM)
  delta_cdm = delta_cdm_tc;
  delta_b   = delta_b_tc;
  v_cdm     = v_cdm_tc;
  v_b       = v_b_tc;
  Phi       = Phi_tc;

  // SET: Photon temperature perturbations (Theta_ell)
  Theta[0] = Theta_tc[0];
  Theta[1] = Theta_tc[1];
  if(polarization){Theta[2] = -8.0*c*k/15.0/Hp/dtau*Theta[1];}
  else{Theta[2] = -20.0*c*k/45.0/Hp/dtau*Theta[1];}
  
  for(int ell = 3; ell < 8; ell++)
  {
    Theta[ell] = -1.0*ell/(2.0*ell+1.0)*c*k/Hp/dtau*Theta[ell-1];
  }
    
  // SET: Photon polarization perturbations (Theta_p_ell)
  if(polarization){
    Theta_p[0] = 5.0/4.0*Theta[2];
    Theta_p[1] = -1.0*c*k/4.0/Hp/dtau*Theta[2];
    Theta_p[2] = 0.25*Theta[2];
    for(int ell = 3; ell < 8; ell++)
    {
      Theta_p[ell] = -1.0*ell/(2.0*ell+1.0)*c*k/Hp/dtau*Theta_p[ell-1];
    }    
  }

  // SET: Neutrino perturbations (N_ell)
  if(neutrinos){
     // n_ell_neutrinos_tc = 8    
    for(int ell = 0; ell < 8; ell++)
    {
      Nu[ell] = Nu_tc[ell];
    }    
  }

  return y;
}

//====================================================
// The time when tight coupling end
//====================================================

double Perturbations::get_tight_coupling_time(const double k) const{
  double x_tight_coupling_end = 0.0;

  //=============================================================================
  // TODO: compute and return x for when tight coupling ends
  // Remember all the three conditions in Callin
  // =============================================================================  
  int npts = 1000;
  Vector x_arr = Utils::linspace(x_start, 0, npts);
  for(int i=0; i<npts; i++){
    double x = x_arr[i];
    double cond1 = abs(rec->dtaudx_of_x(x)/10.0);
    double cond2 = abs(rec->dtaudx_of_x(x)*cosmo->Hp_of_x(x)/(k*Constants.c*10.0));
    double cond3 = rec->Xe_of_x(x)/0.99;

    if((cond1 < 1) || (cond2 < 1) || (cond3 < 1)){
      return x;
    }
  }
  printf("ERROR: Could not find condition for end of tight coupling regime.");
  return -9999999;
}

//====================================================
// After integrsating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(){
  Utils::StartTiming("source");

  //=============================================================================
  // TODO: Make the x and k arrays to evaluate over and use to make the splines
  //=============================================================================
  // ...
  // ...
  Vector k_array(n_k); //n_k defined in the header
  Vector log_k_array = Utils::linspace(log(k_min), log(k_max), n_k);
  for(int i=0; i<n_k; i++){
    k_array[i] = exp(log_k_array[i]);}
  
  Vector x_array = Utils::linspace(x_start, x_end, n_x);
  
  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());
  Vector SE_array(k_array.size() * x_array.size());

  // Compute source functions
  for(auto ix = 0; ix < x_array.size(); ix++){
    const double x = x_array[ix];
    for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];

      // NB: This is the format the data needs to be stored 
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      const int index = ix + n_x * ik;

      //=============================================================================
      // TODO: Compute the source functions
      //=============================================================================
      // Fetch all the things we need...
      const double c    = Constants.c;
      double a          = exp(x);
      double H0         = cosmo->get_H0();
      double Hp         = cosmo->Hp_of_x(x);
      double dHp        = cosmo->dHpdx_of_x(x);
      double ddHp       = cosmo->ddHpddx_of_x(x);
      double eta        = cosmo->eta_of_x(x);
      double eta0       = cosmo->eta_of_x(0.0);
      double tau        = rec->tau_of_x(x);
      double dtau       = rec->dtaudx_of_x(x);
      double ddtau      = rec->ddtauddx_of_x(x);
     
      double g_tilde    = rec->g_tilde_of_x(x);
      double dg_tilde   = rec->dgdx_tilde_of_x(x);
      double ddg_tilde  = rec->ddgddx_tilde_of_x(x);
      
      double v_b        = get_v_b(x, k);
      double dv_bdx     = get_dv_bdx(x, k);
      
      double Pi         = get_Pi(x,k);
      double dPidx      = get_dPidx(x, k);
      double ddPiddx    = get_ddPiddx(x, k);
      double Theta0     = get_Theta(x, k, 0);
      double Psi        = get_Psi(x,k);
      double dPsidx     = get_dPsidx(x, k);
      double dPhidx     = get_dPhidx(x, k);

      // Temperatur source
      double SW, ISW, Dop, Quad;
      SW = g_tilde*(Theta0 + Psi + 0.25*Pi);
      ISW = exp(-tau)*(dPsidx - dPhidx);
      Dop = 1.0/c/k*(dHp*g_tilde*v_b + Hp*dg_tilde*v_b + Hp*g_tilde*dv_bdx);
      Quad = 0.75/(c*k*c*k) * (dHp*dHp + Hp*ddHp)*g_tilde*Pi + 3.0*Hp*dHp*(dg_tilde*Pi + g_tilde*dPidx) + Hp*Hp*(ddg_tilde*Pi + 2.0*dg_tilde*dPidx + g_tilde*ddPiddx);

      ST_array[index] = SW + ISW - Dop + Quad;

      // Polarization source
      if(Constants.polarization){
        double delta_eta = eta0 - eta;
      if(delta_eta/Constants.Mpc < 1e-19){
        SE_array[index] = 0;
        }
      else{
        SE_array[index] = 0.75*g_tilde*Pi/pow(k*(delta_eta),2);
        }
      }
    }
  }

  // Spline the source functions
  std::cout << "S_E: " << SE_array[0] << SE_array[100]<<SE_array[1000]<<"\n";
  ST_spline.create (x_array, k_array, ST_array, "Source_Temp_x_k");
  if(Constants.polarization){
    SE_spline.create (x_array, k_array, SE_array, "Source_Pol_x_k");
  }

  Utils::EndTiming("source");
}

//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const bool neutrinos          = Constants.neutrinos;
  const bool polarization       = Constants.polarization;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  const double &delta_b         =  y[Constants.ind_deltab_tc];
  const double &v_cdm           =  y[Constants.ind_vcdm_tc];
  const double &v_b             =  y[Constants.ind_vb_tc];
  const double &Phi             =  y[Constants.ind_Phi_tc];
  const double *Theta           = &y[Constants.ind_start_theta_tc];
  const double *Nu              = &y[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm_tc];
  double &dv_bdx          =  dydx[Constants.ind_vb_tc];
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];
  double *dThetadx        = &dydx[Constants.ind_start_theta_tc];
  double *dNudx           = &dydx[Constants.ind_start_nu_tc];

  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================

  const double c    = Constants.c;
  double a          = exp(x);
  double H0         = cosmo->get_H0();
  double OmegaB0    = cosmo->get_OmegaB(0.0);
  double OmegaCDM0  = cosmo->get_OmegaCDM(0.0);
  double OmegaR0    = cosmo->get_OmegaR(0.0);
  double OmegaNu0   = cosmo->get_OmegaNu(0.0);
  double Hp         = cosmo->Hp_of_x(x);
  double eta        = cosmo->eta_of_x(x);
  double dtau       = rec->dtaudx_of_x(x);
  double ddtau      = rec->ddtauddx_of_x(x);
  double dHp        = cosmo->dHpdx_of_x(x);
  double Theta2;
    
  if(polarization)
  {Theta2 = -8.0*c*k/15.0/Hp/dtau*Theta[1];}
  else
  {Theta2 = -20.0*c*k/45.0/Hp/dtau*Theta[1];}
  
  double Psi      = -1.0*Phi - 12.0*H0*H0/pow(c*k*a, 2)*(OmegaR0*Theta2 + OmegaNu0*Nu[2]); // Theta[2] = 0??
  double R = 4.0*OmegaR0/3.0/OmegaB0/a;
  double q = (-1.0*((1.0-R)*dtau+(1.0+R)*ddtau)*(3.0*Theta[1]+v_b)-c*k/Hp*Psi+(1.0-dHp/Hp)*c*k/Hp*(-1.0*Theta[0] +2*Theta2)-c*k/Hp*dThetadx[0])/((1.0+R)*dtau+dHp/Hp-1);
  // SET: Scalar quantities (Phi, delta, v, ...)
  dPhidx = Psi - pow(c*k, 2)/(3.0*Hp*Hp)*Phi + H0*H0/(2.0*Hp*Hp)*(OmegaCDM0/a*delta_cdm + OmegaB0/a*delta_b + 4.0*OmegaR0/a/a*Theta[0] + 4.0*OmegaNu0/a/a*Nu[0]);
  ddelta_cdmdx = c*k/Hp*v_cdm - 3*dPhidx;
  dv_cdmdx = -v_cdm - c*k/Hp*Psi;
  ddelta_bdx = c*k/Hp*v_b - 3*dPhidx;
  dv_bdx = 1.0/(1.0+R)*(-1.0*v_b - c*k/Hp*Psi + R*(q + c*k/Hp*(-1*Theta[0] + 2*Theta2) - c*k/Hp*Psi));  // set Theta_2 = 0!!!??

  // SET: Photon multipoles (Theta_ell)
  dThetadx[0] = -1.0*c*k/Hp*Theta[1] - dPhidx;
  dThetadx[1] = 1.0/3.0*(q - dv_bdx);

  // SET: Neutrino mutlipoles (Nu_ell)
  if(neutrinos){
  
  dNudx[0] = -1.0*c*k/Hp*Nu[1] - dPhidx;
  dNudx[1] = c*k/3.0/Hp*Nu[0] - 2.0*c*k/3.0/Hp*Nu[2] + c*k/3.0/Hp*Psi;
  
    for(int ell = 2; ell < 7; ell++)
    {
      dNudx[ell] = 1.0*ell*c*k/(2.0*ell + 1.0)/Hp*Nu[ell-1] - (ell+1)*c*k/(2.0*ell+1.0)/Hp*Nu[ell+1];
    }
    
  dNudx[7] = c*k/Hp*Nu[6] - (7+1)*c/eta/Hp*Nu[7];
  }
  return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================

int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){
  
  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Index and number of the different quantities
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm];
  const double &delta_b         =  y[Constants.ind_deltab];
  const double &v_cdm           =  y[Constants.ind_vcdm];
  const double &v_b             =  y[Constants.ind_vb];
  const double &Phi             =  y[Constants.ind_Phi];
  const double *Theta           = &y[Constants.ind_start_theta];
  const double *Theta_p         = &y[Constants.ind_start_thetap];
  const double *Nu              = &y[Constants.ind_start_nu];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm];
  double &dv_bdx          =  dydx[Constants.ind_vb];
  double &dPhidx          =  dydx[Constants.ind_Phi];
  double *dThetadx        = &dydx[Constants.ind_start_theta];
  double *dTheta_pdx      = &dydx[Constants.ind_start_thetap];
  double *dNudx           = &dydx[Constants.ind_start_nu];

  double a          = exp(x);
    const double c    = Constants.c;
  // Cosmological parameters and variables
  double H0         = cosmo->get_H0();
  double OmegaB0    = cosmo->get_OmegaB(0.0);
  double OmegaCDM0  = cosmo->get_OmegaCDM(0.0);
  double OmegaR0    = cosmo->get_OmegaR(0.0);
  double OmegaNu0   = cosmo->get_OmegaNu(0.0);
  double Hp         = cosmo->Hp_of_x(x);
  double eta        = cosmo->eta_of_x(x);

  // Recombination variables
  
  double dtau       = rec->dtaudx_of_x(x);
  double ddtau      = rec->ddtauddx_of_x(x);

  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================
  double Psi      = -1*Phi - 12*H0*H0/pow(c*k*a, 2)*(OmegaR0*Theta[2] + OmegaNu0*Nu[2]); // Theta[2] = 0??
  double R = 4*OmegaR0/3.0/OmegaB0/a;
  
  // SET: Scalar quantities (Phi, delta, v, ...)
  dPhidx = Psi - pow(c*k, 2)/(3*Hp*Hp)*Phi + H0*H0/(2.0*Hp*Hp)*(OmegaCDM0/a*delta_cdm + OmegaB0/a*delta_b + 4*OmegaR0/a/a*Theta[0] + 4*OmegaNu0/a/a*Nu[0]);
  ddelta_cdmdx = c*k/Hp*v_cdm - 3*dPhidx;
  dv_cdmdx = -v_cdm - c*k/Hp*Psi;
  ddelta_bdx = c*k/Hp*v_b - 3*dPhidx;
  dv_bdx = -1*v_b - c*k/Hp*Psi + dtau*R*(3*Theta[1] + v_b);

  // SET: Photon multipoles (Theta_ell)
  dThetadx[0] = -1*c*k/Hp*Theta[1] - dPhidx;
  dThetadx[1] = c*k/3.0/Hp*Theta[0] -2.0*c*k/3.0/Hp*Theta[2] + c*k/3.0/Hp*Psi + dtau*(Theta[1] + 1.0/3.0*v_b);
  
  double II = Theta[2] + Theta_p[0] + Theta_p[2];
  for(int ell = 2; ell<7; ell++)
  {
    double kd = (ell == 2) ? 1 : 0;  // kronecker delta for l==2.
    dThetadx[ell] = 1.0*ell*c*k/(2.0*ell+1.0)/Hp*Theta[ell-1] - (ell+1)/(2.0*ell+1.0)*c*k/Hp*Theta[ell+1] + dtau*(Theta[ell] - 0.1*II*kd);
  }
  dThetadx[7] = c*k/Hp*Theta[6] - c*(7+1)/Hp/eta*Theta[7] + dtau*Theta[7];
  
  // SET: Photon polarization multipoles (Theta_p_ell)
  if(polarization){

    dTheta_pdx[0] = -1*c*k/Hp*Theta_p[1] + dtau*(Theta_p[0] - 0.5*II);
    for(int ell=1; ell<7;ell++)
    {
      double kd = (ell == 2) ? 1 : 0;  // kronecker delta for l==2.
      dTheta_pdx[ell] = 1.0*ell*c*k/(2.0*ell+1.0)/Hp*Theta_p[ell-1] - (ell+1)/(2.0*ell+1.0)*c*k/Hp*Theta_p[ell+1] +dtau*(Theta_p[ell] - 0.1*II*kd);
    }
    dTheta_pdx[7] = c*k/Hp*Theta_p[6] - c*(7+1)/Hp/eta*Theta_p[7] + dtau*Theta_p[7];
  }

  // SET: Neutrino mutlipoles (Nu_ell)
  if(neutrinos){
    
    dNudx[0] = -1*c*k/Hp*Nu[1] - dPhidx;
    dNudx[1] = 1.0*c*k/3.0/Hp*Nu[0] - 2.0*c*k/3.0/Hp*Nu[2] + c*k/3.0/Hp*Psi;
  
    for(int ell = 2; ell < 7; ell++)
    {
      dNudx[ell] = 1.0*ell*c*k/(2.0*ell + 1)/Hp*Nu[ell-1] - (ell+1.0)*c*k/(2.0*ell+1.0)/Hp*Nu[ell+1];
    }
    
    dNudx[7] = c*k/Hp*Nu[6] - (7+1.0)*c/eta/Hp*Nu[7];
  }

  return GSL_SUCCESS;
}

//====================================================
// Get methods
//====================================================

double Perturbations::get_delta_cdm(const double x, const double k) const{
  return delta_cdm_spline(x,k);
}
double Perturbations::get_delta_b(const double x, const double k) const{
  return delta_b_spline(x,k);
}
double Perturbations::get_v_cdm(const double x, const double k) const{
  return v_cdm_spline(x,k);
}
double Perturbations::get_v_b(const double x, const double k) const{
  return v_b_spline(x,k);
}
double Perturbations::get_dv_bdx(const double x, const double k) const{
  return v_b_spline.deriv_x(x,k);
}
double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}
double Perturbations::get_dPhidx(const double x, const double k) const{
  return Phi_spline.deriv_x(x,k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}
double Perturbations::get_dPsidx(const double x, const double k) const{
  return Psi_spline.deriv_x(x,k);
}
double Perturbations::get_Pi(const double x, const double k) const{
  return Pi_spline(x,k);
}
double Perturbations::get_dPidx(const double x, const double k) const{
  return Pi_spline.deriv_x(x,k);
}
double Perturbations::get_ddPiddx(const double x, const double k) const{
  return Pi_spline.deriv_xx(x,k);
}
double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}
double Perturbations::get_Source_E(const double x, const double k) const{
  return SE_spline(x,k);
}
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return Theta_spline[ell](x,k);
}
double Perturbations::get_Theta_p(const double x, const double k, const int ell) const{
  return Theta_p_spline[ell](x,k);
}
double Perturbations::get_Nu(const double x, const double k, const int ell) const{
  return Nu_spline[ell](x,k);
}

// Test

//====================================================
// Print some useful info about the class
//====================================================

void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:       " << x_start                << "\n";
  std::cout << "x_end:         " << x_end                  << "\n";
  std::cout << "n_x:     " << n_x              << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:     " << n_k              << "\n";
  if(Constants.polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(Constants.neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm         << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab           << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm             << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb               << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi              << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta      << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta          << "\n";
  if(Constants.polarization){
    std::cout << "ind_start_thetap:   " << Constants.ind_start_thetap   << "\n";
    std::cout << "n_ell_thetap:       " << Constants.n_ell_thetap       << "\n";
  }
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu       << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos    << "\n";
  }
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full       << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc      << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc        << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm_tc          << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb_tc            << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc           << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc   << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc       << "\n";
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc    << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc         << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some results to file for a given value of k
//====================================================

void Perturbations::output(const double k, const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts = 5000;
  auto x_array = Utils::linspace(x_start, x_end, npts);
  auto print_data = [&] (const double x) {
    double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                  << " ";  // 0
    fp << get_delta_cdm(x,k) << " ";  // 1
    fp << get_delta_b(x,k)   << " ";  // 2
    fp << get_v_cdm(x,k)     << " ";  // 3
    fp << get_v_b(x,k)       << " ";  // 4
    fp << get_Phi(x,k)       << " ";  // 5
    fp << get_Theta(x,k,0)   << " ";  // 6
    fp << get_Theta(x,k,1)   << " ";  // 7
    // fp << get_Theta(x,k,2)   << " ";

    fp << get_Psi(x,k)       << " ";  // 8
    fp << get_Pi(x,k)        << " ";  // 9
    fp << get_Source_T(x,k)  << " ";  // 10
    fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " "; //11
    fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " ";//12
    fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " ";//13
    fp << get_Source_E(x,k)  << " ";  // 14
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

