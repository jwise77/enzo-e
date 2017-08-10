// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodRayTracer.cpp
/// @author   John H. Wise (jwise2gatech.edu)
/// @date     2017-08-07
/// @brief    Implements the EnzoMethodRayTracer class

#include "cello.hpp"
#include "enzo.hpp"

#define DEBUG_RAYS

#ifdef DEBUG_RAYS
#  define TRACE_RT(MESSAGE)						\
  CkPrintf ("%s:%d %s\n",						\
	    __FILE__,__LINE__,MESSAGE);				
#else
#  define TRACE_RT(MESSAGE) /* ... */
#endif

//----------------------------------------------------------------------

EnzoMethodRayTracer::EnzoMethodRayTracer 
(
 const FieldDescr * field_descr,
 const ParticleDescr * particle_descr, 
 EnzoConfig * enzo_config
) 
  : Method()
{
  // Initialize default Refresh object

  const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier);
  refresh(ir)->add_all_particles(particle_descr->num_types());
  refresh(ir)->add_all_fields(field_descr->field_count());

  // Parameters initialized in EnzoBlock::initialize()
  
}

//----------------------------------------------------------------------

void EnzoMethodRayTracer::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

}

//----------------------------------------------------------------------

void EnzoMethodRayTracer::compute ( Block * block) throw()
{

  TRACE_RT("compute()");

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

  const EnzoConfig * enzo_config = static_cast<const EnzoConfig*>
    (enzo_block->simulation()->config());

  Particle particle (block->data()->particle());
  Field    field    (block->data()->field());
  
  //enzo_float *density = (enzo_float *) field.values("density");
  enzo_float *ray_count = (enzo_float *) field.values("ray_count");

  TRACE_RT("compute1()");
  double dx, dy, dz;
  double lx, ly, lz;
  double ux, uy, uz;
  block->cell_width(&dx, &dy, &dz);
  block->lower(&lx, &ly, &lz);
  block->upper(&ux, &uy, &uz);

  int mx, my, mz;
  int Nx, Ny, Nz;
  int gx, gy, gz;
  field.size(&Nx, &Ny, &Nz);
  field.dimensions(0, &mx, &my, &mz);
  field.ghost_depth(0, &gx, &gy, &gz);

  // Just use a huge number for now so that all rays track through the
  // entire block. Should be dtPhoton instead of 1e20.
  const enzo_float EndTime = enzo_block->time() + 1e20;
  const double c = 1;  // worry about units later
  const double c_inv = 1.0 / c;
  
  const int rank = enzo_block->rank();

  TRACE_RT("compute2()");

  // Obtain ray attributes
  const int it = particle.type_index("rays");
  const int ia_x = (rank >= 1) ? particle.attribute_index(it, "x") : -1;
  const int ia_y = (rank >= 2) ? particle.attribute_index(it, "y") : -1;
  const int ia_z = (rank >= 3) ? particle.attribute_index(it, "z") : -1;
  const int ia_nx = (rank >= 1) ? particle.attribute_index(it, "normal_x") : -1;
  const int ia_ny = (rank >= 2) ? particle.attribute_index(it, "normal_y") : -1;
  const int ia_nz = (rank >= 3) ? particle.attribute_index(it, "normal_z") : -1;
  const int ia_sx = (rank >= 1) ? particle.attribute_index(it, "source_x") : -1;
  const int ia_sy = (rank >= 2) ? particle.attribute_index(it, "source_y") : -1;
  const int ia_sz = (rank >= 3) ? particle.attribute_index(it, "source_z") : -1;
  const int ia_f = particle.attribute_index(it, "flux");
  const int ia_r = particle.attribute_index(it, "radius");
  const int ia_time = particle.attribute_index(it, "time");

  // Obtain particle strides
  const int ps = particle.stride(it, ia_x);

  double *x = 0, *y = 0, *z = 0;
  double *nx = 0, *ny = 0, *nz = 0;
  double *sx = 0, *sy = 0, *sz = 0;
  double *flux = 0, *radius = 0, *time = 0;
  
  // Simple test (TO BE REMOVED): Create one ray every timestep
  int ib0, ipp0;
  double start[3] = {0.01, 0.01, 0.01};
  double unit[3] = {0.5, 0.866, 0.0};

  if (start[0] >= lx && start[0] <= ux &&
      start[0] >= ly && start[0] <= uy &&
      start[0] >= lz && start[0] <= uz) {

    int part0 = particle.insert_particles(it,1);
    particle.index(part0, &ib0, &ipp0);
    x = (double *) particle.attribute_array(it, ia_x, ib0);
    y = (double *) particle.attribute_array(it, ia_y, ib0);
    z = (double *) particle.attribute_array(it, ia_z, ib0);
    nx = (double *) particle.attribute_array(it, ia_nx, ib0);
    ny = (double *) particle.attribute_array(it, ia_ny, ib0);
    nz = (double *) particle.attribute_array(it, ia_nz, ib0);
    sx = (double *) particle.attribute_array(it, ia_sx, ib0);
    sy = (double *) particle.attribute_array(it, ia_sy, ib0);
    sz = (double *) particle.attribute_array(it, ia_sz, ib0);
    flux = (double *) particle.attribute_array(it, ia_f, ib0);
    time = (double *) particle.attribute_array(it, ia_time, ib0);
    radius = (double *) particle.attribute_array(it, ia_r, ib0);

    x[ipp0*ps] = start[0];
    y[ipp0*ps] = start[1];
    z[ipp0*ps] = start[2];
    nx[ipp0*ps] = unit[0];
    ny[ipp0*ps] = unit[1];
    nz[ipp0*ps] = unit[2];
    sx[ipp0*ps] = start[0];
    sy[ipp0*ps] = start[1];
    sz[ipp0*ps] = start[2];
    flux[ipp0*ps] = 1.0;
    time[ipp0*ps] = enzo_block->time();
    radius[ipp0*ps] = 0.0;

    printf("[0] Created a ray. ib0=%d, ipp0=%d, r0=%g\n", ib0, ipp0, radius[ipp0*ps]);
    
  }
  
  const int nb = particle.num_batches(it);
  printf("[0] le %g %g %g // re = %g %g %g\n", lx, ly, lz, ux, uy, uz);
  printf("[0] np=%d, nb=%d, ia_f=%d\n", particle.num_particles(it,0), nb, ia_f);
  
  TRACE_RT("compute3()");
  
  // Loop over ray particle batches and then individual rays
  for (int ib = 0; ib < nb; ib++) {

    const int np = particle.num_particles(it,ib);
    TRACE_RT("compute4()");
    
    if (rank >= 1) {
      x  = (double *) particle.attribute_array(it, ia_x, ib);
      nx = (double *) particle.attribute_array(it, ia_nx, ib);
      sx = (double *) particle.attribute_array(it, ia_sx, ib);
    }
    if (rank >= 2) {
      y  = (double *) particle.attribute_array(it, ia_y, ib);
      ny = (double *) particle.attribute_array(it, ia_ny, ib);
      sy = (double *) particle.attribute_array(it, ia_sy, ib);
    }
    if (rank >= 3) {
      z  = (double *) particle.attribute_array(it, ia_z, ib);
      nz = (double *) particle.attribute_array(it, ia_nz, ib);
      sz = (double *) particle.attribute_array(it, ia_sz, ib);
    }

    flux = (double *) particle.attribute_array(it, ia_f, ib);
    radius = (double *) particle.attribute_array(it, ia_r, ib);
    time = (double *) particle.attribute_array(it, ia_time, ib);

    if (rank == 1) {
      for (int ip = 0; ip < np; ip++) {

	/* Pre-computing quantities before the ray trace */
	
	// Ray can only go left or right, so we take the sign of the
	// supplied normal
	int ipps = ip*ps;
	int n_sign = SIGN(nx[ipps]);
	int n_dir = (n_sign+1)/2;

	// Current cell in integer units
	int ix = floor((x[ipps] - lx) / dx);
	double fx = lx + ix * dx;  // cell boundary
	
	// On cell boundaries, the index will change in negative directions
	if (x[ipps] == fx)
	  ix += (n_sign-1)/2;

	int keep_walking = 1;
	double oldr, thisr, min_r, ce, nce, dr, ddr, dt;
	while (keep_walking) {
	  
	  // This and the next cell edge
	  ce = lx + ix*dx;
	  nce = lx + (ix+n_dir)*dx;

	  // Radius at the next tranversal
	  oldr = radius[ipps];
	  thisr = nce - sx[ipps];
	  dr = thisr - oldr;
	  // avoid being on the edge
	  dr += OMEGA_TOLERANCE*dx;

	  // Do not transport longer than the timestep
	  ddr = MIN(dr, c*(EndTime - time[ipps]));
	  dt = ddr * c_inv;

	  // Update quantities
	  time[ipps] += dt;
	  flux[ipps] += 0.0;  // attenuate
	  radius[ipps] += ddr;
	  
	  ray_count[gx+ix] += 1.0;
	  x[ipps] += n_dir * ddr;
	  ix += n_dir;

	  // Finished (exit the grid, attenuated, or cdt)?
	  if (ix < gx || ix >= Nx)
	    keep_walking = 0;

	  // tiny number for now. Change to physically motivated one later.
	  else if (flux[ipps] < 1e-20)
	    keep_walking = 0;

	  else if (time[ipps] >= EndTime)
	    keep_walking = 0;
	  
	} // ENDWHILE keep_walking
	
      }  // ENDFOR rays
    } // ENDIF rank == 1

    else if (rank == 2) {
      for (int ip = 0; ip < np; ip++) {

	/* Pre-computing quantities before the ray trace */
	
	// Ray can only go left or right, so we take the sign of the
	// supplied normal
	int ipps = ip*ps;
	int nx_sign = SIGN(nx[ipps]);
	int ny_sign = SIGN(ny[ipps]);

	// Avoid zeros in the normal directions and dividing by zero
	if (nx_sign == 0) { nx[ipps] = ETA_TOLERANCE; nx_sign = 1; }
	if (ny_sign == 0) { ny[ipps] = ETA_TOLERANCE; ny_sign = 1; }

	int nx_dir = (nx_sign+1)/2;
	int ny_dir = (ny_sign+1)/2;

	// Current cell in integer units
	int ix = floor((x[ipps] - lx) / dx);
	int iy = floor((y[ipps] - ly) / dy);
	double fx = lx + ix * dx;  // cell boundary
	double fy = ly + iy * dy;  // cell boundary
	
	// On cell boundaries, the index will change in negative directions
	if (x[ipps] == fx) ix += (nx_sign-1)/2;
	if (y[ipps] == fy) iy += (ny_sign-1)/2;

	// Inverse normals (for determining directions)
	double nx_inv = 1.0 / nx[ipps];
	double ny_inv = 1.0 / ny[ipps];
	
	int i, direction, keep_walking = 1;
	double cex, cey, ncex, ncey, drx, dry;
	double oldr, min_r, dr, ddr, dt;
	while (keep_walking) {
	  
	  // This and the next cell edge
	  cex = lx + ix*dx;
	  cey = ly + iy*dy;
	  ncex = lx + (ix+nx_dir)*dx;
	  ncey = ly + (iy+ny_dir)*dy;

	  // Radius of the next cell edge crossing in each dimension
	  drx = nx_inv * (ncex - sx[ipps]);
	  dry = ny_inv * (ncey - sy[ipps]);

	  // The closest one is the next edge crossing
	  if (drx < dry) {
	    min_r = drx;
	    direction = 0;
	  } else {
	    min_r = dry;
	    direction = 1;
	  }
	  
	  // Radius at the next tranversal
	  oldr = radius[ipps];
	  dr = min_r - oldr;
	  // avoid being on the edge
	  dr += OMEGA_TOLERANCE*dx;

	  // Do not transport longer than the timestep
	  ddr = MIN(dr, c*(EndTime - time[ipps]));
	  dt = ddr * c_inv;

	  // Update quantities
	  time[ipps] += dt;
	  flux[ipps] += 0.0;  // attenuate
	  radius[ipps] += ddr;

	  i = gx+ix + mx*(gy+iy);
	  ray_count[i] += 1.0;
	  x[ipps] += nx[ipps] * ddr;
	  y[ipps] += ny[ipps] * ddr;

	  if (direction == 0)
	    ix += nx_dir;
	  else
	    iy += ny_dir;

	  // Finished (exit the grid, attenuated, or cdt)?
	  if (ix < gx || ix >= Nx || iy < gy || iy >= Ny)
	    keep_walking = 0;

	  // tiny number for now. Change to physically motivated one later.
	  else if (flux[ipps] < 1e-20)
	    keep_walking = 0;

	  else if (time[ipps] >= EndTime)
	    keep_walking = 0;
	  
	} // ENDWHILE keep_walking
	
      }  // ENDFOR rays
    } // ENDIF rank == 2

    else if (rank == 3) {
      for (int ip = 0; ip < np; ip++) {

	/* Pre-computing quantities before the ray trace */
	
	// Ray can only go left or right, so we take the sign of the
	// supplied normal
	int ipps = ip*ps;
	int nx_sign = SIGN(nx[ipps]);
	int ny_sign = SIGN(ny[ipps]);
	int nz_sign = SIGN(nz[ipps]);

	// Avoid zeros in the normal directions and dividing by zero
	if (nx_sign == 0) { nx[ipps] = ETA_TOLERANCE; nx_sign = 1; }
	if (ny_sign == 0) { ny[ipps] = ETA_TOLERANCE; ny_sign = 1; }
	if (nz_sign == 0) { nz[ipps] = ETA_TOLERANCE; nz_sign = 1; }

	int nx_dir = (nx_sign+1)/2;
	int ny_dir = (ny_sign+1)/2;
	int nz_dir = (nz_sign+1)/2;

	// Current cell in integer units
	int ix = floor((x[ipps] - lx) / dx);
	int iy = floor((y[ipps] - ly) / dy);
	int iz = floor((z[ipps] - lz) / dz);
	double fx = lx + ix * dx;  // cell boundary
	double fy = ly + iy * dy;  // cell boundary
	double fz = lz + iz * dz;  // cell boundary
	
	// On cell boundaries, the index will change in negative directions
	if (x[ipps] == fx) ix += (nx_sign-1)/2;
	if (y[ipps] == fy) iy += (ny_sign-1)/2;
	if (z[ipps] == fz) iz += (nz_sign-1)/2;

	// Inverse normals (for determining directions)
	double nx_inv = 1.0 / nx[ipps];
	double ny_inv = 1.0 / ny[ipps];
	double nz_inv = 1.0 / nz[ipps];

	printf("[a] %d %d %d %d %d %d\n", nx_sign, nx_dir, ny_sign, ny_dir,
	       nz_sign, nz_dir);
	printf("[b] %d %d %d %g %g %g\n", ix, iy, iz, fx, fy, fz);
	printf("[c] %d %d // %g %g %g\n", ib, ip, nx_inv, ny_inv, nz_inv);

	int i, direction, keep_walking = 1;
	double cex, cey, cez, ncex, ncey, ncez, drx, dry, drz;
	double oldr, min_r, dr, ddr, dt;
	while (keep_walking) {
	  
	  // This and the next cell edge
	  cex = lx + ix*dx;
	  cey = ly + iy*dy;
	  cez = lz + iz*dz;
	  ncex = lx + (ix+nx_dir)*dx;
	  ncey = ly + (iy+ny_dir)*dy;
	  ncez = lz + (iz+nz_dir)*dz;

	  // Radius of the next cell edge crossing in each dimension
	  drx = nx_inv * (ncex - sx[ipps]);
	  dry = ny_inv * (ncey - sy[ipps]);
	  drz = nz_inv * (ncez - sz[ipps]);

	  // The closest one is the next edge crossing
	  if (drx < dry && drx < drz) {
	    min_r = drx;
	    direction = 0;
	  } else if (dry < drx && dry < drz) {
	    min_r = dry;
	    direction = 1;
	  } else {
	    min_r = drz;
	    direction = 2;
	  }
	  
	  // Radius at the next tranversal
	  oldr = radius[ipps];
	  //thisr = min_r;
	  dr = min_r - oldr;
	  // avoid being on the edge
	  dr += OMEGA_TOLERANCE*dx;

	  // Do not transport longer than the timestep
	  ddr = MIN(dr, c*(EndTime - time[ipps]));
	  dt = ddr * c_inv;

	  // Update quantities
	  time[ipps] += dt;
	  flux[ipps] += 0.0;  // attenuate
	  radius[ipps] += ddr;
	  
	  i = gx+ix + mx*(gy+iy + my*(gz+iz));
	  ray_count[i] += 1.0;
	  x[ipps] += nx[ipps] * ddr;
	  y[ipps] += ny[ipps] * ddr;
	  z[ipps] += nz[ipps] * ddr;

	  printf("[A] %g %g %g // %g %g %g\n", cex, cey, cez, ncex, ncey, ncez);
	  printf("[B] %g %g %g // %g %d\n", drx, dry, drz, min_r, direction);
	  printf("[C] %d %g :: %g %g %g %g\n", i, ray_count[i], oldr,
		 dr, ddr, dt);
	  printf("[D] %g %g %g\n", x[ipps], y[ipps], z[ipps]);
	  

	  if (direction == 0)
	    ix += nx_dir;
	  else if (direction == 1)
	    iy += ny_dir;
	  else
	    iz += nz_dir;

	  // Finished (exit the grid, attenuated, or cdt)?
	  if (ix < 0 || ix >= Nx ||
	      iy < 0 || iy >= Ny ||
	      iz < 0 || iz >= Nz) {
	    keep_walking = 0;
	    printf("Exited grid\n");
	  }
	  // tiny number for now. Change to physically motivated one later.
	  else if (flux[ipps] < 1e-20) {
	    keep_walking = 0;
	    printf("Low flux\n");
	  }
	  else if (time[ipps] >= EndTime) {
	    keep_walking = 0;
	    printf("Time's up!\n");
	  }
	} // ENDWHILE keep_walking
	
      }  // ENDFOR rays
    } // ENDIF rank == 3

  } // ENDFOR blocks
    
  block->compute_done(); 
}

//----------------------------------------------------------------------

double EnzoMethodRayTracer::timestep ( Block * block ) const throw()
{

  TRACE_RT("timestep()");

  const double dt = 0.1;
  return dt;
  
}
