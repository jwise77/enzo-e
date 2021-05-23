// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodRayTracer.cpp
/// @author   John H. Wise (jwise2gatech.edu)
/// @date     2017-08-07
/// @brief    Implements the EnzoMethodRayTracer class

#include "cello.hpp"
#include "enzo.hpp"

#define DEBUG_RAYS

#ifdef DEBUG_RAYS
#  define TRACE_RT(MESSAGE, BLOCK)					\
  CkPrintf ("[%s] %s:%d %s\n",						\
	    (BLOCK != NULL) ? BLOCK->name().c_str() : "root", __FILE__,__LINE__,MESSAGE);				
#else
#  define TRACE_RT(MESSAGE, BLOCK) /* ... */
#endif

//----------------------------------------------------------------------

EnzoMethodRayTracer::EnzoMethodRayTracer 
(
 const EnzoConfig * enzo_config
) 
  : Method()
{
  // Initialize default Refresh object

  cello::simulation()->new_refresh_set_name(ir_post_, name());
  Refresh * refresh = cello::refresh(ir_post_);

  // For now, refresh all particles and fields
  // In the future, trim to the essential data
  refresh->add_all_particles();
  refresh->add_all_fields();
  
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

  TRACE_RT("compute()", block);

  EnzoBlock * enzo_block = enzo::block(block);

  // Setup an quiescence (all work/msgs finished) exit call
  // CkStartQD(CkCallback(CkIndex_EnzoBlock::p_method_rt_end(NULL),
  // 		       block->thisProxy));
  
  if (block->is_leaf()) {
    setup_attributes(enzo_block);
    generate_rays(enzo_block);
    process_local(enzo_block);
  }

  block->compute_done();
  
}

//----------------------------------------------------------------------

void EnzoMethodRayTracer::process_local( EnzoBlock * enzo_block) throw()
{

  TRACE_RT("process_local()", enzo_block);

  // Ray trace
  int n_ray_exit;
  n_ray_exit = trace_rays(enzo_block);

  // Pack and send any rays that have exited the block. Model after
  // the particle refreshes in Cello/control_refresh.cpp.
  
  /*
    for (index_neighbor = ...) {
      packed_rays = ...;
      proxy_block[index_neighbor].p_receive_rays(packed_rays);
    }
   */
  
  return;
  
}

//----------------------------------------------------------------------

void EnzoBlock::p_method_rt_receive_rays(CkReductionMsg * msg)
{

  TRACE_RT("p_method_receive_rays()", this);

  delete msg;

  // Unpack rays


  // Process the rays
  //EnzoMethodRayTracer * method =
  //  static_cast<EnzoMethodRayTracer *> (this->method());
  //int n_ray_exit = method->trace_rays(this);
  
  return;
  
}

//----------------------------------------------------------------------

void EnzoBlock::p_method_rt_end(CkReductionMsg * msg)
{

  TRACE_RT("p_method_raytracer_end()", this);

  delete msg;

  // TODO: call grackle from here, updating all of the cells with
  // radiation
  
  compute_done();
  
  return;
  
}

//----------------------------------------------------------------------

void EnzoMethodRayTracer::setup_attributes ( EnzoBlock * enzo_block) throw()
{

  const int rank = cello::rank();

  // Obtain ray attributes
  Particle particle (enzo_block->data()->particle());
  it = particle.type_index("rays");
  ia_x = (rank >= 1) ? particle.attribute_index(it, "x") : -1;
  ia_y = (rank >= 2) ? particle.attribute_index(it, "y") : -1;
  ia_z = (rank >= 3) ? particle.attribute_index(it, "z") : -1;
  ia_nx = (rank >= 1) ? particle.attribute_index(it, "normal_x") : -1;
  ia_ny = (rank >= 2) ? particle.attribute_index(it, "normal_y") : -1;
  ia_nz = (rank >= 3) ? particle.attribute_index(it, "normal_z") : -1;
  ia_sx = (rank >= 1) ? particle.attribute_index(it, "source_x") : -1;
  ia_sy = (rank >= 2) ? particle.attribute_index(it, "source_y") : -1;
  ia_sz = (rank >= 3) ? particle.attribute_index(it, "source_z") : -1;
  ia_f = particle.attribute_index(it, "flux");
  ia_r = particle.attribute_index(it, "radius");
  ia_time = particle.attribute_index(it, "time");

  // Obtain particle strides
  ps = particle.stride(it, ia_x);

  return;
  
}

int EnzoMethodRayTracer::trace_rays ( EnzoBlock * enzo_block) throw()
{

  Particle particle (enzo_block->data()->particle());
  Field    field    (enzo_block->data()->field());

  // If no rays to trace, exit immediately
  if (particle.num_particles(it) == 0)
    return 0;

  //enzo_float *density = (enzo_float *) field.values("density");
  enzo_float *ray_count = (enzo_float *) field.values("ray_count");

  double dx, dy, dz;
  double lx, ly, lz;
  double ux, uy, uz;
  enzo_block->cell_width(&dx, &dy, &dz);
  enzo_block->lower(&lx, &ly, &lz);
  enzo_block->upper(&ux, &uy, &uz);

  int mx, my, mz;
  int Nx, Ny, Nz;
  int gx, gy, gz;
  field.size(&Nx, &Ny, &Nz);
  field.dimensions(0, &mx, &my, &mz);
  field.ghost_depth(0, &gx, &gy, &gz);

  // Just use sqrt(3) for now so that all rays track through the
  // entire block (corner to corner). Should be dtPhoton instead.
  const enzo_float EndTime = enzo_block->time() + sqrt(3.0);
  const double c = 1;  // worry about units later
  const double c_inv = 1.0 / c;
  
  const int rank = cello::rank();
  double *x = 0, *y = 0, *z = 0;
  double *nx = 0, *ny = 0, *nz = 0;
  double *sx = 0, *sy = 0, *sz = 0;
  double *flux = 0, *radius = 0, *time = 0;

  // Number of rays exiting the block
  int n_ray_exit = 0;
  
  const int nb = particle.num_batches(it);
  printf("[0] le %g %g %g // re = %g %g %g\n", lx, ly, lz, ux, uy, uz);
  printf("[0] np=%d, nb=%d, ia_f=%d\n", particle.num_particles(it,0), nb, ia_f);
  
  TRACE_RT("trace_rays()", enzo_block);
  
  // Loop over ray particle batches and then individual rays
  for (int ib = 0; ib < nb; ib++) {

    const int np = particle.num_particles(it,ib);
    TRACE_RT("trace_rays_loop()", enzo_block);
    
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
	  if (ix < gx || ix >= Nx) {
	    n_ray_exit++;
	    keep_walking = 0;
	  }
	  // tiny number for now. Change to physically motivated one later.
	  else if (flux[ipps] < 1e-20) {
	    keep_walking = 0;
	  }
	  else if (time[ipps] >= EndTime) {
	    keep_walking = 0;
	  }
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
	  if (ix < gx || ix >= Nx || iy < gy || iy >= Ny) {
	    keep_walking = 0;
	    n_ray_exit++;
	  }
	  // tiny number for now. Change to physically motivated one later.
	  else if (flux[ipps] < 1e-20) {
	    keep_walking = 0;
	  }
	  else if (time[ipps] >= EndTime) {
	    keep_walking = 0;
	  }
	  
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
	  printf("[D] %g %g %g, t=%g/%g\n", x[ipps], y[ipps], z[ipps], time[ipps], EndTime);
	  

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
	    n_ray_exit++;
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

    // Check whether any rays have exited the domain.  Need to reset
    // their source positions so that they appear "outside" the domain.

    //     ... domain extents
    double dxm,dym,dzm;
    double dxp,dyp,dzp;
    cello::hierarchy()->lower(&dxm,&dym,&dzm);
    cello::hierarchy()->upper(&dxp,&dyp,&dzp);
    
    for (int ip = 0; ip < np; ip++) {
      int ipps = ip*ps;
      if (x[ipps] < dxm) sx[ipps] += +(dxp - dxm);
      if (x[ipps] > dxp) sx[ipps] += -(dxp - dxm);
      if (rank >= 2) {
	if (y[ipps] < dym) sy[ipps] += +(dyp - dym);
	if (y[ipps] > dyp) sy[ipps] += -(dyp - dym);
      }
      if (rank >= 3) {
	if (z[ipps] < dzm) sz[ipps] += +(dzp - dzm);
	if (z[ipps] > dzp) sz[ipps] += -(dzp - dzm);
      }
    } // ENDFOR particles
    
  } // ENDFOR blocks

  return n_ray_exit;

}

//----------------------------------------------------------------------

void EnzoMethodRayTracer::generate_rays ( EnzoBlock * enzo_block) throw()
{

  TRACE_RT("generate_rays()", enzo_block);
  Particle particle (enzo_block->data()->particle());
  
  const int rank = cello::rank();
  double lx, ly, lz;
  double ux, uy, uz;
  enzo_block->lower(&lx, &ly, &lz);
  enzo_block->upper(&ux, &uy, &uz);

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
  
  return;
  
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------

double EnzoMethodRayTracer::timestep ( Block * block ) const throw()
{

  TRACE_RT("timestep()", block);

  const double dt = 0.1;
  return dt;
  
}
