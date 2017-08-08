// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPpm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoMethodPpm class

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodRayTracer::EnzoMethodRayTracer 
(
 const FieldDescr * field_descr,
 EnzoConfig * enzo_config
) 
  : Method()
{
  // Initialize default Refresh object

  const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier);
  refresh(ir)->add_all_fields(field_descr->field_count());

  // PPM parameters initialized in EnzoBlock::initialize()
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

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

  const EnzoConfig * enzo_config = static_cast<const EnzoConfig*>
    (enzo_block->simulation()->config());

  Particle particle (block->data()->particle());
  Field    field    (block->data()->field());
  
  //enzo_float *density = (enzo_float *) field.values("density");
  enzo_float *ray_count = (enzo_float *) field.values("ray_count");

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
  const int nb = particle.num_batches(it);
  
  // Obtain particle strides
  const int ps = particle.stride(it, ia_x);

  // Loop over ray particle batches and then individual rays
  for (int ib = 0; ib < nb; ib++) {
    double *x = 0, *y = 0, *z = 0;
    double *nx = 0, *ny = 0, *nz = 0;
    double *sx = 0, *sy = 0, *sz = 0;
    double *flux = 0, *radius = 0, *time = 0;

    const int np = particle.num_particles(it,ib);
    
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
    radius = (double *) particle.attribute_array(it, ia_time, ib);
    time = (double *) particle.attribute_array(it, ia_r, ib);

    if (rank == 1) {
      for (int ip = 0; ip < np; ip++) {

	/* Pre-computing quantities before the ray trace */
	
	// Ray can only go left or right, so we take the sign of the
	// supplied normal
	int ipps = ip*ps;
	int n_sign = sign(nx[ipps]);
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
	  // avoid being on the edge
	  thisr = nce - sx[ipps] + ETA_TOLERANCE;
	  dr = thisr - oldr;

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
	  if (ix < gx || ix > Nx)
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
	int nx_sign = sign(nx[ipps]);
	int ny_sign = sign(ny[ipps]);
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
	if (fabs(nx[ipps]) < ETA_TOLERANCE) nx[ipps] = nx_sign * ETA_TOLERANCE;
	if (fabs(ny[ipps]) < ETA_TOLERANCE) ny[ipps] = ny_sign * ETA_TOLERANCE;
	double nx_inv = 1.0 / nx[ipps];
	double ny_inv = 1.0 / ny[ipps];
	
	int i, direction, keep_walking = 1;
	double cex, cey, ncex, ncey, drx, dry;
	double oldr, thisr, min_r, dr, ddr, dt;
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
	  // avoid being on the edge
	  thisr = min_r + ETA_TOLERANCE;
	  dr = thisr - oldr;

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
	  if (ix < gx || ix > Nx || iy < gy || iy > Ny)
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
	int nx_sign = sign(nx[ipps]);
	int ny_sign = sign(ny[ipps]);
	int nz_sign = sign(nz[ipps]);
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
	if (fabs(nx[ipps]) < ETA_TOLERANCE) nx[ipps] = nx_sign * ETA_TOLERANCE;
	if (fabs(ny[ipps]) < ETA_TOLERANCE) ny[ipps] = ny_sign * ETA_TOLERANCE;
	if (fabs(nz[ipps]) < ETA_TOLERANCE) nz[ipps] = nz_sign * ETA_TOLERANCE;
	double nx_inv = 1.0 / nx[ipps];
	double ny_inv = 1.0 / ny[ipps];
	double nz_inv = 1.0 / nz[ipps];
	
	int i, direction, keep_walking = 1;
	double cex, cey, cez, ncex, ncey, ncez, drx, dry, drz;
	double oldr, thisr, min_r, dr, ddr, dt;
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
	  // avoid being on the edge
	  thisr = min_r + ETA_TOLERANCE;
	  dr = thisr - oldr;

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

	  if (direction == 0)
	    ix += nx_dir;
	  else if (direction == 1)
	    iy += ny_dir;
	  else
	    iz += nz_dir;

	  // Finished (exit the grid, attenuated, or cdt)?
	  if (ix < gx || ix > Nx ||
	      iy < gy || iy > Ny ||
	      iz < gz || iz > Nz)
	    keep_walking = 0;

	  // tiny number for now. Change to physically motivated one later.
	  else if (flux[ipps] < 1e-20)
	    keep_walking = 0;

	  else if (time[ipps] >= EndTime)
	    keep_walking = 0;
	  
	} // ENDWHILE keep_walking
	
      }  // ENDFOR rays
    } // ENDIF rank == 3

  } // ENDFOR blocks
    
  block->compute_done(); 
}

//----------------------------------------------------------------------
