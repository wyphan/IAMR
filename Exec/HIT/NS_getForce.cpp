
#include <NavierStokesBase.H>
#include <AMReX_BLFort.H>

// For now, define pi here, but maybe later make iamr_constants.H
namespace {
  constexpr Real Pi    = 3.141592653589793238462643383279502884197;
  constexpr Real TwoPi = 2.0 * 3.141592653589793238462643383279502884197;
}


using namespace amrex;

//
// Virtual access function for getting the forcing terms for the
// velocities and scalars.  The base version computes a buoyancy.
//
// NOTE: This function returns a rho weighted source term.
//
// For conservative (i.e. do_mom_diff=1, do_cons_trac=1), velocities
// are integrated according to the equation
//
//     ui_t + uj ui_j = S_ui        ===> tforces = rho S_ui
//
// and scalars psi where (psi = rho q = Scal) as
//
//     psi_t + (uj psi)_j = S_psi   ===> tforces = S_psi = rho S_q
//
// For non-conservative, this rho-weighted source term will get divided
// by rho in the predict_velocity, velocity_advection, scalar_advection,
// and advection_update routines.
//
// For temperature (which is always non-conservative), we evolve
//
//     dT/dt - U dot grad T = [del dot lambda grad T + S_T] / (rho*c_p)
//     ===> tforces =  S_T/c_p
//
//
// For user-defined forcing, this means
//   - For conservative variables, the force term computed here gets used
//     as-is
//   - For non-conservative variables, the force term computed here is
//     divided by rho before use
//

void
NavierStokesBase::getForce (FArrayBox&       force,
                            const Box&       bx,
                            int              scomp,
                            int              ncomp,
                            const Real       time,
                            const FArrayBox& Vel,
                            const FArrayBox& Scal,
                            int              scalScomp,
                            const MFIter&    mfi)
{

   const Real* VelDataPtr  = Vel.dataPtr();
   const Real* ScalDataPtr = Scal.dataPtr(scalScomp);

   const Real  grav     = gravity;
   const int*  f_lo     = force.loVect();
   const int*  f_hi     = force.hiVect();
   const int*  v_lo     = Vel.loVect();
   const int*  v_hi     = Vel.hiVect();
   const int*  s_lo     = Scal.loVect();
   const int*  s_hi     = Scal.hiVect();

   if (ParallelDescriptor::IOProcessor() && getForceVerbose) {
      amrex::Print() << "NavierStokesBase::getForce(): Entered..." << std::endl 
                     << "time      = " << time << std::endl
                     << "scomp     = " << scomp << std::endl
                     << "ncomp     = " << ncomp << std::endl
                     << "scalScomp = " << scalScomp << std::endl;

      if (scomp==0)
	if  (ncomp==3) amrex::Print() << "Doing velocities only" << std::endl;
	else           amrex::Print() << "Doing all components" << std::endl;
      else if (scomp==3)
	if  (ncomp==1) amrex::Print() << "Doing density only" << std::endl;
	else           amrex::Print() << "Doing all scalars" << std::endl;
      else if (scomp==4) amrex::Print() << "Doing tracer only" << std::endl;
      else               amrex::Print() << "Doing individual scalar" << std::endl;

      amrex::Print() << "NavierStokesBase::getForce(): Filling Force on box:"
		     << bx << std::endl;
#if (AMREX_SPACEDIM == 3)
      amrex::Print() << "NavierStokesBase::getForce(): Force Domain:" << std::endl;
      amrex::Print() << "(" << f_lo[0] << "," << f_lo[1] << "," << f_lo[2] << ") - "
                     << "(" << f_hi[0] << "," << f_hi[1] << "," << f_hi[2] << ")" << std::endl;
      amrex::Print() << "NavierStokesBase::getForce(): Vel Domain:" << std::endl;
      amrex::Print() << "(" << v_lo[0] << "," << v_lo[1] << "," << v_lo[2] << ") - "
                     << "(" << v_hi[0] << "," << v_hi[1] << "," << v_hi[2] << ")" << std::endl;
      amrex::Print() << "NavierStokesBase::getForce(): Scal Domain:" << std::endl;
      amrex::Print() << "(" << s_lo[0] << "," << s_lo[1] << "," << s_lo[2] << ") - "
                     << "(" << s_hi[0] << "," << s_hi[1] << "," << s_hi[2] << ")" << std::endl;
#else
      amrex::Print() << "NavierStokesBase::getForce(): Force Domain:" << std::endl;
      amrex::Print() << "(" << f_lo[0] << "," << f_lo[1] << ") - "
                     << "(" << f_hi[0] << "," << f_hi[1] << ")" << std::endl;
      amrex::Print() << "NavierStokesBase::getForce(): Vel Domain:" << std::endl;
      amrex::Print() << "(" << v_lo[0] << "," << v_lo[1] << ") - "
                     << "(" << v_hi[0] << "," << v_hi[1] << ")" << std::endl;
      amrex::Print() << "NavierStokesBase::getForce(): Scal Domain:" << std::endl;
      amrex::Print() << "(" << s_lo[0] << "," << s_lo[1] << ") - "
                     << "(" << s_hi[0] << "," << s_hi[1] << ")" << std::endl;
#endif

      Vector<Real> velmin(AMREX_SPACEDIM), velmax(AMREX_SPACEDIM);
      Vector<Real> scalmin(NUM_SCALARS), scalmax(NUM_SCALARS);
      for (int n=0; n<AMREX_SPACEDIM; n++) {
          velmin[n]= 1.e234;
          velmax[n]=-1.e234;
      }
      int ix = v_hi[0]-v_lo[0]+1;
      int jx = v_hi[1]-v_lo[1]+1;
#if (AMREX_SPACEDIM == 3)
      int kx = v_hi[2]-v_lo[2]+1;
      for (int k=0; k<kx; k++) {
#endif
         for (int j=0; j<jx; j++) {
            for (int i=0; i<ix; i++) {
               for (int n=0; n<AMREX_SPACEDIM; n++) {
#if (AMREX_SPACEDIM == 3)
                  int cell = ((n*kx+k)*jx+j)*ix+i;
#else
                  int cell = (n*jx+j)*ix+i;
#endif
                  Real v = VelDataPtr[cell];
                  if (v<velmin[n]) velmin[n] = v;
                  if (v>velmax[n]) velmax[n] = v;
               }
            }
         }
#if (AMREX_SPACEDIM == 3)
      }
#endif
      for (int n=0; n<AMREX_SPACEDIM; n++) 
         amrex::Print() << "Vel  " << n << " min/max " 
                        << velmin[n] << " / " << velmax[n] << std::endl;

      for (int n=0; n<NUM_SCALARS; n++) {
         scalmin[n]= 1.e234;
         scalmax[n]=-1.e234;
      }
      ix = s_hi[0]-s_lo[0]+1;
      jx = s_hi[1]-s_lo[1]+1;
#if (AMREX_SPACEDIM == 3)
      kx = s_hi[2]-s_lo[2]+1;
      for (int k=0; k<kx; k++) {
#endif
         for (int j=0; j<jx; j++) {
            for (int i=0; i<ix; i++) {
               for (int n=0; n<NUM_SCALARS; n++) {
#if (AMREX_SPACEDIM == 3)
                  int cell = ((n*kx+k)*jx+j)*ix+i;
#else
                  int cell = (n*jx+j)*ix+i;
#endif
                  Real s = ScalDataPtr[cell];
                  if (s<scalmin[n]) scalmin[n] = s;
                  if (s>scalmax[n]) scalmax[n] = s;
               }
            }
         }
#if (AMREX_SPACEDIM == 3)
      }
#endif
      for (int n=0; n<NUM_SCALARS; n++) 
         amrex::Print() << "Scal " << n << " min/max " << scalmin[n] 
                        << " / " << scalmax[n] << std::endl;
   } //end if(getForceVerbose)

   //
   // Here's the meat
   //
   // Velocity forcing
   //
   if ( scomp<AMREX_SPACEDIM ){
     AMREX_ALWAYS_ASSERT(scomp==Xvel);
     AMREX_ALWAYS_ASSERT(ncomp>=AMREX_SPACEDIM);
   }

   if ( scomp==Xvel ){
     //
     // TODO: add some switch for user-supplied/problem-dependent forcing
     //
     auto const& frc  = force.array(scomp);
     auto const& scal = Scal.array(scalScomp);

     if ( std::abs(grav) > 0.0001) {
       amrex::ParallelFor(bx, [frc, scal, grav]
       AMREX_GPU_DEVICE(int i, int j, int k) noexcept
       {
	 frc(i,j,k,0) = Real(0.0);
#if ( AMREX_SPACEDIM == 2 )
         frc(i,j,k,1) = grav*scal(i,j,k,0);
#elif ( AMREX_SPACEDIM == 3 )
         frc(i,j,k,1) = Real(0.0);
         frc(i,j,k,2) = grav*scal(i,j,k,0);
#endif
       });
     }
     else {
       force.setVal<RunOn::Gpu>(0.0, bx, Xvel, AMREX_SPACEDIM);
     }

#ifdef USE_TURBULENT_FORCING
     //
     // Homogeneous Isotropic Turbulence
     //

     // For now, only works in 3D and without tiling
     AMREX_ALWAYS_ASSERT(AMREX_SPACEDIM==3);
     AMREX_ALWAYS_ASSERT(bx==force.box())
     
     // force array bounds
     // FIXME -- think about how bx (the box we want to fill) may not be the same as the force box!
     int ilo = f_lo[0];
     int jlo = f_lo[1];
     int klo = f_lo[2];
     int ihi = f_hi[0];
     int jhi = f_hi[1];
     int khi = f_hi[2];

     // Physical coordinates of the lower left corner of the domain
     auto const& problo = geom.ProbLoArray();
     // Physical coordinates of the upper right corner of the domain
     auto const& probhi = geom.ProbHiArray();

     Real Lx = probhi[0]-problo[0];
     Real Ly = probhi[1]-problo[1];
     Real Lz = probhi[2]-problo[2];
     Real Lmin = min(Lx,Ly,Lz);
     
     int xstep = static_cast<int>(Lx/Lmin+0.5);
     int ystep = static_cast<int>(Ly/Lmin+0.5);
     int zstep = static_cast<int>(Lz/Lmin+0.5);

     Real kappaMax = nmodes/Lmin + 1.0e-8;

     auto const&  dx = geom.CellSizeArray();
     Real hx = dx[0];
     Real hy = dx[1];
     Real hz = dx[2]; 
     
     // coarse cell size
     Real ff_hx = hx*ff_factor;
     Real ff_hy = hy*ff_factor;
     Real ff_hz = hz*ff_factor;

     // coarse bounds (without accounting for ghost cells)
     // FIXME -- think about how bx (the box we want to fill) may not be the same as the force box!
     int ff_ilo = ilo/ff_factor;
     int ff_jlo = jlo/ff_factor;
     int ff_klo = klo/ff_factor;

     int ff_ihi = (ihi+1)/ff_factor;
     int ff_jhi = (jhi+1)/ff_factor;
     int ff_khi = (khi+1)/ff_factor;
           
     // adjust for ghost cells
     if (ilo < (ff_ilo*ff_factor)) {
       ff_ilo=ff_ilo-1;
     }
     if (jlo < (ff_jlo*ff_factor)) {
       ff_jlo=ff_jlo-1;
     }
     if (klo < (ff_klo*ff_factor)) {
       ff_klo=ff_klo-1;
     }
     if (ihi == (ff_ihi*ff_factor)) {
       ff_ihi=ff_ihi+1;
     }
     if (jhi == (ff_jhi*ff_factor)) {
       ff_jhi=ff_jhi+1;
     }
     if (khi == (ff_khi*ff_factor)) {
       ff_khi=ff_khi+1;
     }
     
     // allocate coarse force array     
     Box ffbx(IntVect(ff_ilo, ff_jlo, ff_klo), IntVect(ff_ihi, ff_jhi, ff_khi));
     // not sure if want elixir, gpu::sync, or async_arena here...
     FArrayBox ff_force(ffbx,AMREX_SPACEDIM);
     const auto& ffarr = ff_force.array();

     RealBox loc = RealBox(bx,geom.CellSize(),geom.ProbLo());
     const auto& loc_lo = loc.lo();     
     amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> xlo = {AMREX_D_DECL(loc_lo[0],loc_lo[1],loc_lo[2])};
     
     // Construct node-based coarse forcing
     amrex::ParallelFor(ffbx, [ = ]
     AMREX_GPU_DEVICE (int i, int j, int k ) noexcept
     {
       Real z = xlo(2) + ff_hz*(k-ff_klo);
       Real y = xlo(1) + ff_hy*(j-ff_jlo);
       Real x = xlo(0) + ff_hx*(i-ff_ilo);

       for (int n = 0; n < AMREX_SPACEDIM; n++)
	 ffarr(i,j,k,n) = 0.0;
       
       for ( int kx = mode_start*xstep; kx <= nmodes*xstep; kx += xstep) {
	 for ( int ky = mode_start*ystep; ky <= nmodes*ystep; ky += ystep) {
	   for ( int kz = mode_start*zstep; kz <= nmodes*zstep; kz += zstep)
	   {
	     Real kappa = sqrt( (kx*kx)/(Lx*Lx) + (ky*ky)/(Ly*Ly) + (kz*kz)/(Lz*Lz) );

	     if (kappa <= kappaMax)
	     {
	       Real xT = cos(FTX(kx,ky,kz)*time + TAT(kx,ky,kz));

	       if ( div_free_force )
	       {
		 ffarr(i,j,k,0) += xT * 
		      ( FAZ(kx,ky,kz)*TwoPi*(ky/Ly) 
			*   sin(TwoPi*kx*x/Lx+FPZX(kx,ky,kz)) 
			*   cos(TwoPi*ky*y/Ly+FPZY(kx,ky,kz)) 
			*   sin(TwoPi*kz*z/Lz+FPZZ(kx,ky,kz)) 
			- FAY(kx,ky,kz)*TwoPi*(kz/Lz) 
			*   sin(TwoPi*kx*x/Lx+FPYX(kx,ky,kz)) 
			*   sin(TwoPi*ky*y/Ly+FPYY(kx,ky,kz)) 
			*   cos(TwoPi*kz*z/Lz+FPYZ(kx,ky,kz)) );
		   
		 ffarr(i,j,k,1) += xT * 
		      ( FAX(kx,ky,kz)*TwoPi*(kz/Lz) 
			*   sin(TwoPi*kx*x/Lx+FPXX(kx,ky,kz)) 
			*   sin(TwoPi*ky*y/Ly+FPXY(kx,ky,kz)) 
			*   cos(TwoPi*kz*z/Lz+FPXZ(kx,ky,kz)) 
			- FAZ(kx,ky,kz)*TwoPi*(kx/Lx) 
			*   cos(TwoPi*kx*x/Lx+FPZX(kx,ky,kz)) 
			*   sin(TwoPi*ky*y/Ly+FPZY(kx,ky,kz)) 
			*   sin(TwoPi*kz*z/Lz+FPZZ(kx,ky,kz)) );
		 
		 ffarr(i,j,k,2) += xT * 
		      ( FAY(kx,ky,kz)*TwoPi*(kx/Lx) 
			*   cos(TwoPi*kx*x/Lx+FPYX(kx,ky,kz)) 
			*   sin(TwoPi*ky*y/Ly+FPYY(kx,ky,kz)) 
			*   sin(TwoPi*kz*z/Lz+FPYZ(kx,ky,kz)) 
			- FAX(kx,ky,kz)*TwoPi*(ky/Ly) 
			*   sin(TwoPi*kx*x/Lx+FPXX(kx,ky,kz)) 
			*   cos(TwoPi*ky*y/Ly+FPXY(kx,ky,kz)) 
			*   sin(TwoPi*kz*z/Lz+FPXZ(kx,ky,kz)) );
	       }
	       else
	       {

		 ffarr(i,j,k,0) += xT*FAX(kx,ky,kz)*cos(TwoPi*kx*x/Lx+FPX(kx,ky,kz)) 
		      *                             sin(TwoPi*ky*y/Ly+FPY(kx,ky,kz)) 
                      *                             sin(TwoPi*kz*z/Lz+FPZ(kx,ky,kz));

		 ffarr(i,j,k,1) += xT*FAY(kx,ky,kz)*sin(TwoPi*kx*x/Lx+FPX(kx,ky,kz)) 
		      *                             cos(TwoPi*ky*y/Ly+FPY(kx,ky,kz)) 
		      *                             sin(TwoPi*kz*z/Lz+FPZ(kx,ky,kz));

		 ffarr(i,j,k,2) += xT*FAZ(kx,ky,kz)*sin(TwoPi*kx*x/Lx+FPX(kx,ky,kz)) 
		      *                             sin(TwoPi*ky*y/Ly+FPY(kx,ky,kz)) 
		      *                             cos(TwoPi*kz*z/Lz+FPZ(kx,ky,kz));
	       }
	     } 
	   }
	 }
       }
       
       for ( int kx = mode_start; kx <= nmodes*xstep; kx++) {
	 for ( int ky = mode_start; ky <= nmodes*ystep; ky++) {
	   for ( int kz = mode_start; kz <= nmodes*zstep; kz++)
	   {
	     Real kappa = sqrt( (kx*kx)/(Lx*Lx) + (ky*ky)/(Ly*Ly) + (kz*kz)/(Lz*Lz) );

	     if (kappa.le.kappaMax)
	     {
	       Real xT = cos(FTX(kx,ky,kz)*time + TAT(kx,ky,kz));

	       if ( div_free_force )
	       {
		 ffarr(i,j,k,0) += xT * 
		      ( FAZ(kx,ky,kz)*TwoPi*(ky/Ly) 
			*   sin(TwoPi*kx*x/Lx+FPZX(kx,ky,kz)) 
			*   cos(TwoPi*ky*y/Ly+FPZY(kx,ky,kz)) 
			*   sin(TwoPi*kz*z/Lz+FPZZ(kx,ky,kz)) 
			- FAY(kx,ky,kz)*TwoPi*(kz/Lz) 
			*   sin(TwoPi*kx*x/Lx+FPYX(kx,ky,kz)) 
			*   sin(TwoPi*ky*y/Ly+FPYY(kx,ky,kz)) 
			*   cos(TwoPi*kz*z/Lz+FPYZ(kx,ky,kz)) );
		   
		 ffarr(i,j,k,1) += xT * 
		      ( FAX(kx,ky,kz)*TwoPi*(kz/Lz) 
			*   sin(TwoPi*kx*x/Lx+FPXX(kx,ky,kz)) 
			*   sin(TwoPi*ky*y/Ly+FPXY(kx,ky,kz)) 
			*   cos(TwoPi*kz*z/Lz+FPXZ(kx,ky,kz)) 
			- FAZ(kx,ky,kz)*TwoPi*(kx/Lx) 
			*   cos(TwoPi*kx*x/Lx+FPZX(kx,ky,kz)) 
			*   sin(TwoPi*ky*y/Ly+FPZY(kx,ky,kz)) 
			*   sin(TwoPi*kz*z/Lz+FPZZ(kx,ky,kz)) );
		 
		 ffarr(i,j,k,2) += xT * 
		      ( FAY(kx,ky,kz)*TwoPi*(kx/Lx) 
			*   cos(TwoPi*kx*x/Lx+FPYX(kx,ky,kz)) 
			*   sin(TwoPi*ky*y/Ly+FPYY(kx,ky,kz)) 
			*   sin(TwoPi*kz*z/Lz+FPYZ(kx,ky,kz)) 
			- FAX(kx,ky,kz)*TwoPi*(ky/Ly) 
			*   sin(TwoPi*kx*x/Lx+FPXX(kx,ky,kz)) 
			*   cos(TwoPi*ky*y/Ly+FPXY(kx,ky,kz)) 
			*   sin(TwoPi*kz*z/Lz+FPXZ(kx,ky,kz)) );
	       }
	       else
	       {

		 ffarr(i,j,k,0) += xT*FAX(kx,ky,kz)*cos(TwoPi*kx*x/Lx+FPX(kx,ky,kz)) 
		      *                             sin(TwoPi*ky*y/Ly+FPY(kx,ky,kz)) 
                      *                             sin(TwoPi*kz*z/Lz+FPZ(kx,ky,kz));

		 ffarr(i,j,k,1) += xT*FAY(kx,ky,kz)*sin(TwoPi*kx*x/Lx+FPX(kx,ky,kz)) 
		      *                             cos(TwoPi*ky*y/Ly+FPY(kx,ky,kz)) 
		      *                             sin(TwoPi*kz*z/Lz+FPZ(kx,ky,kz));

		 ffarr(i,j,k,2) += xT*FAZ(kx,ky,kz)*sin(TwoPi*kx*x/Lx+FPX(kx,ky,kz)) 
		      *                             sin(TwoPi*ky*y/Ly+FPY(kx,ky,kz)) 
		      *                             cos(TwoPi*kz*z/Lz+FPZ(kx,ky,kz));
	       }
	     } 
	   }
	 }
       }
       
     });

     // Now interpolate onto fine grid
     // bx is the box we want to fill and may be smaller than force
     auto const& frc  = force.array(scomp);
     auto const& dens = Scal.array(scalScomp);

     amrex::ParallelFor(bx, AMREX_SPACEDIM, [ = ]
     AMREX_GPU_DEVICE (int i, int j, int k, int n ) noexcept
     { 
       Real zd = ( hz*(k-klo + 0.5) - ff_hz*(ff_k-ff_klo) )/ff_hz; 
       Real yd = ( hy*(j-jlo + 0.5) - ff_hy*(ff_j-ff_jlo) )/ff_hy;
       Real xd = ( hx*(i-ilo + 0.5) - ff_hx*(ff_i-ff_ilo) )/ff_hx;

       int ff_k = k/ff_factor;
       int ff_j = j/ff_factor;
       int ff_i = i/ff_factor;


       Real ff00 =  ffarr(ff_i  ,ff_j  ,ff_k  ,n) * (1. - xd)
	          + ffarr(ff_i+1,ff_j  ,ff_k  ,n) * xd;
       Real ff01 =  ffarr(ff_i  ,ff_j  ,ff_k+1,n) * (1. - xd)
	          + ffarr(ff_i+1,ff_j  ,ff_k+1,n) * xd;       
       Real ff10 =  ffarr(ff_i  ,ff_j+1,ff_k  ,n) * (1. - xd)
	          + ffarr(ff_i+1,ff_j+1,ff_k  ,n) * xd;
       Real ff11 =  ffarr(ff_i  ,ff_j+1,ff_k+1,n) * (1. - xd)
	          + ffarr(ff_i+1,ff_j+1,ff_k+1,n) * xd;

       Real ff =  ( ff00*(one-yd)+ff10*yd ) * (1. - zd)
	        + ( ff01*(one-yd)+ff11*yd ) * zd
       
       frc(i,j,k,n) += dens(i,j,k,0) * ff;
     });

     amrex::Gpu::synchronize();

#endif // Turbulent forcing

   }
   
   
   //
   // Scalar forcing
   //
   if ( scomp >= AMREX_SPACEDIM ) {
     // Doing only scalars
     force.setVal<RunOn::Gpu>(0.0, bx, 0, ncomp);

     //
     // Or create user-defined forcing.
     // Recall we compute a density-weighted forcing term.
     //
     // auto const& frc  = force.array();
     // amrex::ParallelFor(bx, ncomp, [frc]
     // AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept
     // {
     //          frc(i,j,k,n) = ;
     //          frc(i,j,k,n) *= rho;
     // });
   }
   else if ( scomp+ncomp > AMREX_SPACEDIM) {
     // Doing scalars with vel
     force.setVal<RunOn::Gpu>(0.0, bx, Density, ncomp-Density);

     //
     // Or create user-defined forcing.
     // Recall we compute a density-weighted forcing term.
     //
     // auto const& frc  = force.array(Density);
     // amrex::ParallelFor(bx, ncomp-Density, [frc]
     // AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept
     // {
     //          frc(i,j,k,n) = ;
     //          frc(i,j,k,n) *= rho;
     // });
   }
     
   if (ParallelDescriptor::IOProcessor() && getForceVerbose) {
      Vector<Real> forcemin(ncomp);
      Vector<Real> forcemax(ncomp);
      for (int n=0; n<ncomp; n++) {
         forcemin[n]= 1.e234;
         forcemax[n]=-1.e234;
      }

      int ix = f_hi[0]-f_lo[0]+1;
      int jx = f_hi[1]-f_lo[1]+1;
#if (AMREX_SPACEDIM == 3)
      int kx = f_hi[2]-f_lo[2]+1;
      for (int k=0; k<kx; k++) {
#endif
         for (int j=0; j<jx; j++) {
            for (int i=0; i<ix; i++) {
               for (int n=0; n<ncomp; n++) {
#if (AMREX_SPACEDIM == 3)
                  int cell = ((n*kx+k)*jx+j)*ix+i;
#else
                  int cell = (n*jx+j)*ix+i;
#endif
                  Real f = force.dataPtr()[cell];
                  if (f<forcemin[n]) forcemin[n] = f;
                  if (f>forcemax[n]) forcemax[n] = f;
               }
            }
         }
#if (AMREX_SPACEDIM == 3)
      }
#endif
      for (int n=0; n<ncomp; n++) 
         amrex::Print() << "Force " << n+scomp << " min/max " << forcemin[n] 
                        << " / " << forcemax[n] << std::endl;

      amrex::Print() << "NavierStokesBase::getForce(): Leaving..." 
                     << std::endl << "---" << std::endl;
   }
}
