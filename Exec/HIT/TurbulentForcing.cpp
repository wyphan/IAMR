#include<TurbulentForcing_params.H>
#include <AMReX_ParmParse.H>

void
TurbulentForcing::read_turbulent_forcing_params()
{
  // read in parameters for turbulent forcing
  ParmParse pp("turb");

  // Make user set nmodes. This must match up with nmodes used in making MagicFile
  TurbulentForcing::nmodes = -1;
  pp.qet("nmodes", TurbulentForcing::nmodes);

  TurbulentForcing::div_free_force = true;
  pp.query("div_free_force", TurbulentForcing::div_free_force);

  TurbulentForcing::ff_factor = 4;
  pp.query("ff_factor", TurbulentForcing::ff_factor);

  TurbulentForcing::mode_start = 0;
  pp.query("mode_start", TurbulentForcing::mode_start);
  
  // Now read in values from the MagicFile.
  
}
