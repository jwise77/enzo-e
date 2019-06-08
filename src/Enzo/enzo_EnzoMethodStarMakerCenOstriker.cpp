// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPpm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoMethodPpm class

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodStarMakerCenOstriker::EnzoMethodStarMakerCenOstriker 
(
 const EnzoConfig * enzo_config
)
  : Method()
{
  // Initialize default Refresh object

  const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier,
                             enzo_sync_id_method_cen_ostriker);
  refresh(ir)->add_all_fields();

  // PPM parameters initialized in EnzoBlock::initialize()
}

//----------------------------------------------------------------------

void EnzoMethodStarMakerCenOstriker::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

}

//----------------------------------------------------------------------

void EnzoMethodStarMakerCenOstriker::compute ( Block * block) throw()
{
  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

  const EnzoConfig * enzo_config = enzo::config();

  int star_counter = 0;
  // Are we at the highest level?
  if (block->is_leaf()) {

    Particle particle (block->data()->particle());
    Field field = enzo_block->data()->field();

    double dx, dy, dz;

    block->cell_width(&dx, &dy, &dz);

    double lx, ly, lz;
    block->lower(&lx, &ly, &lz);

    // declare particle position arrays

    const int it = particle.type_index ("stars");

    const int ia_m = particle.attribute_index (it,"mass");
    const int ia_x = particle.attribute_index (it,"x");
    const int ia_y = particle.attribute_index (it,"y");
    const int ia_z = particle.attribute_index (it,"z");
    const int ia_vx = particle.attribute_index (it,"vx");
    const int ia_vy = particle.attribute_index (it,"vy");
    const int ia_vz = particle.attribute_index (it,"vz");


    int ib=0;  // batch counter
    int ipp=0;

    double * mass = 0;
    double * px = 0;
    double * py = 0;
    double * pz = 0;
    double * vx = 0;
    double * vy = 0;
    double * vz = 0;

    const int ps  = particle.stride(it,ia_m);

    int rank = cello::rank();

    enzo_float * density    = (enzo_float *)field.values("density");
    enzo_float * velocity_x = (rank >= 1) ? 
      (enzo_float *)field.values("velocity_x") : NULL;
    enzo_float * velocity_y = (rank >= 2) ? 
      (enzo_float *)field.values("velocity_y") : NULL;
    enzo_float * velocity_z = (rank >= 3) ? 
      (enzo_float *)field.values("velocity_z") : NULL;

    int gx,gy,gz;
    field.ghost_depth (0,&gx,&gy,&gz);

    int mx,my,mz;
    field.dimensions (0,&mx,&my,&mz);

    int nx,ny,nz;
    field.size (&nx,&ny,&nz);

    for (int iz=gz; iz<nz+gz; iz++) {
      for (int iy=gy; iy<ny+gy; iy++) {
        for (int ix=gx; ix<nx+gx; ix++) {

          int i = ix + mx*(iy + my*iz);

          if (density[i] >= enzo_config->star_maker_co_density_threshold) {

            // insert a new particle and grab the particle index
            int my_particle = particle.insert_particles (it,1);

            // find the block index (ib) and position in block (ipp)
            particle.index(my_particle, &ib, &ipp);

            // pointer to mass array in block
            mass = (double *) particle.attribute_array(it,ia_m,ib);

            // update particle mass
            mass[ipp*ps] = density[i] * dx*dy*dz * enzo_config->star_maker_co_efficiency;
            density[i] *= (1 - enzo_config->star_maker_co_efficiency);

            // repeat for particle positions and velocities
            px = (double *) particle.attribute_array(it, ia_x, ib);
            py = (double *) particle.attribute_array(it, ia_y, ib);
            pz = (double *) particle.attribute_array(it, ia_z, ib);
            px[ipp*ps] = lx + (ix + 0.5) * dx;
            py[ipp*ps] = ly + (iy + 0.5) * dy;
            pz[ipp*ps] = lz + (iz + 0.5) * dz;

            vx = (double *) particle.attribute_array(it, ia_vx, ib);
            vy = (double *) particle.attribute_array(it, ia_vy, ib);
            vz = (double *) particle.attribute_array(it, ia_vz, ib);
            vx[ipp*ps] =  100.0 ; //velocity_x[i] ;
            vy[ipp*ps] =  100.0 ; //velocity_y[i] ;
            vz[ipp*ps] = velocity_z[i] ;

            star_counter++;
          }

        }
      }
    }



  }

  block->compute_done(); 
}

//----------------------------------------------------------------------
