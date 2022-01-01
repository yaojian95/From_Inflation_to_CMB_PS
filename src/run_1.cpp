#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"
#include "PowerSpectrum.h"

int main(int argc, char **argv){
  Utils::StartTiming("Everything");

  //=========================================================================
  // Parameters
  //=========================================================================

  // Background parameters
  double h           = 0.674;
  double OmegaB      = 0.05;
  double OmegaCDM    = 0.267;
  double OmegaK      = 0.0;
  double Neff        = 3.046;
  double TCMB        = 2.7255;

  // Recombination parameters
  double Yp            = 0.245;
  bool Reionization    = true;
  double z_reion       = 8;
  double delta_z_reion = 0.5;
  double z_Hereion     = 3.5;  
  double delta_z_Hereion = 0.5;

  // Power-spectrum parameters
  double A_s         = 2.1e-9;
  double n_s         = 0.965;
  double kpivot_mpc  = 0.05;

  //=========================================================================
  // Module I
  //=========================================================================

  // Set up and solve the background
  BackgroundCosmology cosmo(h, OmegaB, OmegaCDM, OmegaK, Neff, TCMB);
  cosmo.solve();
  cosmo.info();
    
  const double x_min = -10.0;
  const double x_max =  0.0;
  const int    n_pts =  100;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

    
    double H = cosmo.Hp_of_x(0);
      std::cout << "H:        " << H      << "\n";
  
  // Output background evolution quantities
  // cosmo.output("cosmology.txt");
  //  return 0; 
   //=========================================================================
  // Module II
  //=========================================================================
  
  // Solve the recombination history
  RecombinationHistory rec(&cosmo, Reionization, z_reion, delta_z_reion, z_Hereion, delta_z_Hereion, Yp);
  rec.solve();
  rec.info();

  // Output recombination quantities
//   rec.output("recombination.txt");
  
//   // Remove when module is completed
//   return 0;

  
  //=========================================================================
  // Module III
  //=========================================================================
 
  // Solve the perturbations
  Perturbations pert(&cosmo, &rec);
  pert.solve();
          std::cout << "You are at here!!!" << std::endl;
  pert.info();
  
  // Output perturbation quantities
  double kvalue = 0.01 / Constants.Mpc;
  // pert.output(kvalue, "perturbations_k0.01.txt");
  
  // Remove when module is completed
  // return 0;
  
  //=========================================================================
  // Module IV
  //=========================================================================
  std::cout << "Begin PowerSpectrum!" << std::endl;
  PowerSpectrum power(&cosmo, &rec, &pert, A_s, n_s, kpivot_mpc);
  //------- adjust for other ----------//
    // PowerSpectrum power(&cosmo, &rec, &pert);
  std::cout << "Begin PowerSpectrum.solve!" << std::endl;
  power.solve();
  power.output("cells.txt");
  
  // Remove when module is completed
  return 0;
  Utils::EndTiming("Everything");
}