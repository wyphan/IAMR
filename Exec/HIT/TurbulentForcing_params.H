#ifndef _TurbulentForcing_params_H_
#define _TurbulentForcing_params_H_

namespace TurbulentForcing {
  // function to read in the parameters
  void read_turbulent_forcing_params();
  
  // factor by which to reduce sampling for faster performance
  extern AMREX_GPU_MANAGED int ff_factor;
  // make the forcing divergence free?
  extern AMREX_GPU_MANAGED bool div_free_force;
  // how many modes to use
  extern AMREX_GPU_MANAGED int nmodes;
  // don't use any modes below mode_start. We probably don't need this
  extern AMREX_GPU_MANAGED int mode_start;

  constexpr int array_size = 33;
  extern amrex::GpuArray<amrex::GpuArray<amrex::Real, array_size>, AMREX_SPACEDIM> FTX;
  extern amrex::GpuArray<amrex::GpuArray<amrex::Real, array_size>, AMREX_SPACEDIM> TAT;
  extern amrex::GpuArray<amrex::GpuArray<amrex::Real, array_size>, AMREX_SPACEDIM> FAX;
  extern amrex::GpuArray<amrex::GpuArray<amrex::Real, array_size>, AMREX_SPACEDIM> FAY;
  extern amrex::GpuArray<amrex::GpuArray<amrex::Real, array_size>, AMREX_SPACEDIM> FAZ;
  extern amrex::GpuArray<amrex::GpuArray<amrex::Real, array_size>, AMREX_SPACEDIM> FPX;
  extern amrex::GpuArray<amrex::GpuArray<amrex::Real, array_size>, AMREX_SPACEDIM> FPY;
  extern amrex::GpuArray<amrex::GpuArray<amrex::Real, array_size>, AMREX_SPACEDIM> FPZ;
  extern amrex::GpuArray<amrex::GpuArray<amrex::Real, array_size>, AMREX_SPACEDIM> FPXX;
  extern amrex::GpuArray<amrex::GpuArray<amrex::Real, array_size>, AMREX_SPACEDIM> FPXY;
  extern amrex::GpuArray<amrex::GpuArray<amrex::Real, array_size>, AMREX_SPACEDIM> FPXZ;
  extern amrex::GpuArray<amrex::GpuArray<amrex::Real, array_size>, AMREX_SPACEDIM> FPYX;
  extern amrex::GpuArray<amrex::GpuArray<amrex::Real, array_size>, AMREX_SPACEDIM> FPYY;
  extern amrex::GpuArray<amrex::GpuArray<amrex::Real, array_size>, AMREX_SPACEDIM> FPYZ;
  extern amrex::GpuArray<amrex::GpuArray<amrex::Real, array_size>, AMREX_SPACEDIM> FPZX;
  extern amrex::GpuArray<amrex::GpuArray<amrex::Real, array_size>, AMREX_SPACEDIM> FPZY;
  extern amrex::GpuArray<amrex::GpuArray<amrex::Real, array_size>, AMREX_SPACEDIM> FPZZ;
}
#endif
