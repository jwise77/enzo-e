// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoInitialDisk.cpp
/// @author	John Wise (jwise@physics.gatech.edu)
/// @date	Fri 24 Feb 2023
/// @brief	Implementation of a rotating disk around a central point mass
///
/// Description TBW

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------
#define NO_DEBUG_PERFORMANCE
//----------------------------------------------------------------------

EnzoInitialDisk::EnzoInitialDisk
(const EnzoConfig * config) throw()
: Initial(config->initial_cycle, config->initial_time)
{
  // read in parameter settings from config and
  // set corresponding member variables
  for(int i = 0; i < 3; i ++){
    center_position_[i] = config->initial_IG_center_position[i];
    bfield_[i] = config->initial_IG_bfield[i];
  }

  EnzoUnits * enzo_units = enzo::units();
  FieldDescr * field_descr = cello::field_descr();

  // Store variables locally with code units for convenience
  //     this is not strictly necessary, but makes routine a little cleaner

  this->scale_height_            = config->initial_IG_scale_height * enzo_constants::kpc_cm /
                                   enzo_units->length();
  this->scale_length_	           = config->initial_IG_scale_length * enzo_constants::kpc_cm /
                                   enzo_units->length();
  this->disk_mass_               = config->initial_IG_disk_mass * enzo_constants::mass_solar /
                                   enzo_units->mass();
  this->gas_fraction_            = config->initial_IG_gas_fraction;
  this->disk_temperature_        = config->initial_IG_disk_temperature /
                                   enzo_units->kelvin_per_energy_units();
  this->disk_metal_fraction_     = config->initial_IG_disk_metal_fraction;
  this->gas_halo_mass_           = config->initial_IG_gas_halo_mass * enzo_constants::mass_solar /
                                   enzo_units->mass();
  this->gas_halo_temperature_    = config->initial_IG_gas_halo_temperature /
                                   enzo_units->kelvin_per_energy_units();
  this->gas_halo_metal_fraction_ = config->initial_IG_gas_halo_metal_fraction;
  this->gas_halo_density_        = config->initial_IG_gas_halo_density /
                                   enzo_units->density();
  this->gas_halo_radius_         = config->initial_IG_gas_halo_radius * enzo_constants::kpc_cm /
                                   enzo_units->length();

  // on / off settings for IC
  this->use_gas_particles_       = config->initial_IG_use_gas_particles;
  this->live_dm_halo_            = config->initial_IG_live_dm_halo;
  this->stellar_disk_            = config->initial_IG_stellar_disk;
  this->stellar_bulge_           = config->initial_IG_stellar_bulge;
  this->analytic_velocity_       = config->initial_IG_analytic_velocity;

  // AE: NOTE: This is a bit of a hack at the moment -
  //           this grouping should be registered elsewhere (I think??)
  //           but is apparently not - hard coding this for now
  //
  //           This likely means this will need to be done for
  //           color fields as well, but maybe that is taken care of
  //           properly
  ParticleDescr * particle_descr = cello::particle_descr();
  if (this->stellar_disk_ || this->stellar_bulge_) {
    particle_descr->check_particle_attribute("star","mass");
    particle_descr->groups()->add("star","is_gravitating"); // hack
  }
  if (this->live_dm_halo_) {
    particle_descr->check_particle_attribute("dark","mass");
    particle_descr->groups()->add("dark","is_gravitating");
  }

  // Compute halo density / mass
  if ((this->gas_halo_density_ == 0.0) && this->gas_halo_mass_ > 0)
  {
    if (this->gas_halo_radius_ > 0.0)
    {
      this->gas_halo_density_ = this->gas_halo_mass_ /
                    (4.0 * cello::pi / 3.0 * pow( this->gas_halo_radius_, 3));
    } else { // assume full box of LUxLUxLU
      this->gas_halo_density_ = this->gas_halo_mass_; // In code density  mass / 1**3
    }
  } // AE: Need an error check to make sure these are set correctly

  // gather other parameters not associated with this IC
  this->uniform_density_    = config->field_uniform_density / enzo_units->density();
  this->dual_energy_        = field_descr->is_field("internal_energy") &&
                              field_descr->is_field("total_energy");
  this->gamma_              = enzo::fluid_props()->gamma();
  this->mu_                 = enzo::fluid_props()->mol_weight();

  // read in data for initialization
  this->ntypes_            = 0;              // num of IC particle types
  this->ndim_              = cello::rank();  //
  this->nparticles_        = 0;              // max num of IC particles per type

  //
  // Parameters for model of recent star formation
  //    goes back to recent_SF_age and selects stars in bins up to present
  //    where N_select_i = recent_SF_SFR * recent_SF_bin_size
  //
  //    grabs them randomly over galaxy with exponential decay probability
  //    and assigns lifetimes. This allows for IMMEDIATE initial stellar feedback
  //    in the galaxy
  //
  this->include_recent_SF   = config->initial_IG_include_recent_SF;
  this->recent_SF_start     = config->initial_IG_recent_SF_start; // age of oldest newly formed star
  this->recent_SF_end       = config->initial_IG_recent_SF_end;
  this->recent_SF_bin_size  = config->initial_IG_recent_SF_bin_size; // time resolution
  this->recent_SF_SFR       = config->initial_IG_recent_SF_SFR;
  this->recent_SF_seed      = config->initial_IG_recent_SF_seed;

  if (this->include_recent_SF) srand(this->recent_SF_seed);

  this->tiny_number_       = 1.0E-10;

  this->ReadInVcircData();        // circular velocity curve
  this->ReadParticles();          // IC particles - only if available / used

  return;
}

//----------------------------------------------------------------------

void EnzoInitialDisk::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP; //

  Initial::pup(p); //

  PUParray(p,center_position_,3);
  PUParray(p, bfield_,3);
  p | scale_length_;
  p | disk_mass_;
  p | gas_fraction_;
  p | disk_temperature_;
  p | disk_metal_fraction_;
  p | gas_halo_mass_;
  p | gas_halo_temperature_;
  p | gas_halo_metal_fraction_;
  p | gas_halo_density_;
  p | gas_halo_radius_;

  p | uniform_density_;
  p | dual_energy_;
  p | gamma_;
  p | mu_;
  p | tiny_number_;

  p | use_gas_particles_;
  p | live_dm_halo_;
  p | stellar_bulge_;
  p | stellar_disk_;
  p | analytic_velocity_;

  p | ntypes_;
  p | ndim_;
  p | nparticles_;

  p | include_recent_SF;
  p | recent_SF_start;
  p | recent_SF_end;
  p | recent_SF_bin_size;
  p | recent_SF_SFR;
  p | recent_SF_seed;

  if (p.isUnpacking()){
    // we know aboutntype, ndim, and nparticles..
    // allocate particle IC arrays
    allocateParticles();
  }

  // Pup particle IC arrays element by element
  for (int k = 0; k < ntypes_; k++){
    for (int j = 0; j < ndim_; j ++){
      PUParray(p, particleIcPosition[k][j], nparticles_);
      PUParray(p, particleIcVelocity[k][j], nparticles_);
    }
    PUParray(p, particleIcMass[k], nparticles_);
    PUParray(p, particleIcCreationTime[k], nparticles_);
    PUParray(p, particleIcLifetime[k], nparticles_);
  }
  PUParray(p, particleIcTypes, ntypes_);

  // vector can just be used without anything special
  p | particleIcFileNames;

  return;
}

void EnzoInitialDisk::enforce_block
(
 Block * block,
 const Hierarchy * hierarchy
 ) throw()
{

  //
  // Make sure we can operate on this block
  //
  // remove once parent-to-child particle ICs working
  if (!block->is_leaf()) {
    block->initial_done();
    return;
  }

  Timer timer;
  timer.start();

  ASSERT("EnzoInitialDisk",
         "Block does not exist",
         block != NULL);

  Particle particle = block->data()->particle();
  InitializeParticles(block, &particle);


#ifdef CONFIG_USE_GRACKLE
   grackle_field_data grackle_fields_;
   const EnzoMethodGrackle * grackle_method = enzo::grackle_method();
   grackle_method->setup_grackle_fields(block, &grackle_fields_);
#endif

  if (this->use_gas_particles_){
    this->InitializeGasFromParticles(block);
  } else {
    this->InitializeExponentialGasDistribution(block);
  }


  // Update temperature field if it exists
  Field field = block->data()->field();
  const EnzoConfig * enzo_config = enzo::config();

#ifdef CONFIG_USE_GRACKLE
  /* Assign grackle chemistry fields to default fractions based on density */
  const size_t num_method = enzo_config->method_list.size();

  for (size_t index_method=0; index_method < num_method ; index_method++) {
    std::string name = enzo_config->method_list[index_method];

    if (name == "grackle"){

      grackle_method->update_grackle_density_fields(block, &grackle_fields_);
    }
  }
#endif

  enzo_float * temperature = field.is_field("temperature") ?
                   (enzo_float*) field.values("temperature") : NULL;

  if (temperature) {
    EnzoComputeTemperature compute_temperature(enzo::fluid_props(),
                                               enzo_config->physics_cosmology);

    compute_temperature.compute(block);
  }

  // now initialize particles
  // remove once parent-to-child particle ICs working
//  Particle particle = block->data()->particle();
//  InitializeParticles(block, &particle);

#ifdef CONFIG_USE_GRACKLE
  grackle_method->delete_grackle_fields(&grackle_fields_);
#endif

  return;
}

void EnzoInitialDisk::InitializeGasDistribution(Block * block){

  // Initialize gas distribution using a double exponential (like AGORA)

  EnzoUnits * enzo_units = enzo::units();
  const EnzoConfig * enzo_config = enzo::config();
  Field field = block->data()->field();

  //
  // Get Grid and Field parameters
  //

  // Block sizes (excluding ghost zones)
  int nx,ny,nz;
  field.size(&nx,&ny,&nz);

  // Ghost zone depth (gx + gx total ghost zones in x dimension)
  int gx, gy, gz;
  field.ghost_depth(0,&gx,&gy,&gz);

  // total number of cells in each dimension
  int mx = nx + 2*gx;
  int my = ny + 2*gy;
  int mz = nz + 2*gz;

  // Grab information about grid properties
  //     Cell min coords (xm,ym,zm)  - of min active zone cell (not ghost)
  //     Cell max coords (xp,yp,zp)
  //     Cell widths (hx,hy,hz)
  double xm, ym, zm, xp, yp, zp, hx, hy, hz;
  block->data()->lower(&xm,&ym,&zm);
  block->data()->upper(&xp,&yp,&zp);
  block->cell_width(&hx,&hy,&hz);

  // Get Fields
  enzo_float * d           = (enzo_float *) field.values ("density");
  enzo_float * p           = field.is_field("pressure") ?
                             (enzo_float *) field.values ("pressure") : NULL;
  enzo_float * a3[3]       = { (enzo_float *) field.values("acceleration_x"),
                               (enzo_float *) field.values("acceleration_y"),
                               (enzo_float *) field.values("acceleration_z")};
  enzo_float * v3[3]       = { (enzo_float *) field.values("velocity_x"),
                               (enzo_float *) field.values("velocity_y"),
                               (enzo_float *) field.values("velocity_z")};

  enzo_float * te          = (enzo_float *) field.values("total_energy");
  enzo_float * ge          = (enzo_float *) field.values("internal_energy");
  enzo_float * pot         = (enzo_float *) field.values("potential");

  enzo_float * metal       = field.is_field("metal_density") ?
                             (enzo_float *) field.values("metal_density") : NULL;

  //
  // Now lets calculate some physical properties of the galaxy and halo
  //
  double rho_zero = this->disk_mass_ * this->gas_fraction_ / (4.0 * cello::pi)/
      (pow((this->scale_length_),2)*(this->scale_height_));

  double halo_gas_energy = this->gas_halo_temperature_ / this->mu_ / (this->gamma_ -1);

  double disk_gas_energy = this->disk_temperature_ / this->mu_ / (this->gamma_ -1) ;
  
  // initialize fields to something
  for (int iz=0; iz<mz; iz++){
    for (int iy=0; iy<my; iy++){
      for (int ix=0; ix<mx; ix++){
        int i = INDEX(ix,iy,iz,mx,my);
        d[i]   = this->uniform_density_;
        te[i]  = halo_gas_energy;
        p[i]   = (this->gamma_ - 1.0) * te[i] * d[i];

        if(pot) pot[i] = 0.0;

        if (this->dual_energy_)
        {
          ge[i] = halo_gas_energy;
        }

        for (int dim = 0; dim < 3; dim++){
          a3[dim][i] = 0.0;
          v3[dim][i] = 0.0;
        }

        if(metal) metal[i] = this->tiny_number_ * d[i];

      }
    }
  } // end loop over all cells for background values


  for (int iz=0; iz<mz; iz++){
    // compute z coordinate of cell center
    double z = zm + (iz - gz + 0.5)*hz - this->center_position_[2];
    z *= enzo_units->length(); // convert to cgs units

    for (int iy=0; iy<my; iy++){
      double y = ym + (iy - gy + 0.5)*hy - this->center_position_[1];
      y *= enzo_units->length();

      for (int ix=0; ix<mx; ix++){
        double x = xm + (ix - gx + 0.5)*hx - this->center_position_[0];
        x *= enzo_units->length();

        // 1D index of current cell
        int i = INDEX(ix,iy,iz,mx,my);

        // compute spherical and cylindrical radii (in cgs)
        double radius = sqrt(x*x + y*y + z*z);
        double r_cyl  = sqrt(x*x + y*y);

        // compute the disk density (in code units)
        double disk_density = this->gauss_mass(rho_zero, x/enzo_units->length(), 
          y/enzo_units->length(), z/enzo_units->length(), hx) / (hx*hy*hz);

        if ((this->gas_halo_density_ * this->gas_halo_temperature_ > disk_density*this->disk_temperature_) &&
            (radius < this->gas_halo_radius_*enzo_units->length())){
          // in the halo, set the halo properties
          d[i]  = this->gas_halo_density_;
          te[i] = halo_gas_energy;
          p[i]  = (this->gamma_ - 1.0) * te[i] * d[i];

          if (this->dual_energy_) {
            ge[i] = halo_gas_energy;
          }

          // set metal fraction here
          if(metal) metal[i] = this->gas_halo_metal_fraction_ * d[i];

        }
        else if ( radius < this->gas_halo_radius_*enzo_units->length() ) // we are in the disk
        {
          // in the disk, set the disk properties
          d[i]   = disk_density;

          double vcirc = 0.0;
          if (this->analytic_velocity_){
            double rcore_cgs = enzo_config->method_background_acceleration_core_radius * enzo_constants::kpc_cm;
            double rvir_cgs = enzo_config->method_background_acceleration_DM_mass_radius * enzo_constants::kpc_cm;
            double Mvir_cgs = enzo_config->method_background_acceleration_DM_mass * enzo_constants::mass_solar;
            ASSERT1("Enzo::InitialIsolatedGalaxy", "DM halo mass (=%e g) must be positive and in units of solar masses", Mvir_cgs, (Mvir_cgs > 0));

            double conc = rvir_cgs  / rcore_cgs;
            double   rx = r_cyl / rvir_cgs;

            vcirc = (std::log(1.0 + conc*rx) - (conc*rx)/(1.0+conc*rx))/
                      (std::log(1.0 + conc) - (conc / (1.0 + conc))) / rx;
            vcirc = std::sqrt(vcirc * enzo_constants::grav_constant * Mvir_cgs / rvir_cgs);

          } else {
            vcirc = this->InterpolateVcircTable(r_cyl);
          }

          /* Assume counter-clockwise rotation */
          v3[0][i] = -(vcirc*(y/r_cyl))/enzo_units->velocity();
          v3[1][i] =  (vcirc*(x/r_cyl))/enzo_units->velocity();
          v3[2][i] = 0.0;

          te[i]  = disk_gas_energy;
          p[i]   = (this->gamma_ - 1.0) * te[i] * d[i];

          for (int dim = 0; dim < 3; dim++) // AE: do I need a check for not ppm?
          {
            te[i] += 0.5*v3[dim][i]*v3[dim][i];
          }

          if (this->dual_energy_)
          {
            ge[i] = disk_gas_energy;
          }

          if(metal) metal[i] = this->disk_metal_fraction_ * d[i];

        } // end disk / halo check

      }
    }
  } // end loop over all cells

  block->initial_done();

  return;
}
