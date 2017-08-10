// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodRayTracer.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Thu Apr  1 16:14:38 PDT 2010
/// @brief    [\ref Enzo] Implementation of Enzo PPM hydro method

#ifndef ENZO_ENZO_METHOD_RAYTRACER_HPP
#define ENZO_ENZO_METHOD_RAYTRACER_HPP

class EnzoMethodRayTracer : public Method {

  /// @class    EnzoMethodRayTracer
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate Enzo's ray tracer

public: // interface

  /// Create a new EnzoMethodRayTracer object
  EnzoMethodRayTracer(const FieldDescr * field_descr,
		      const ParticleDescr * particle_descr, 
		      EnzoConfig * enzo_config);

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodRayTracer);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodRayTracer (CkMigrateMessage *m)
    : rays_per_cell_(5.1)
  {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  
  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual double timestep( Block * block) const throw();

  void setup_attributes( EnzoBlock * enzo_block) throw();
  void trace_rays( EnzoBlock * enzo_block) throw();
  void generate_rays( EnzoBlock * enzo_block) throw();

  virtual std::string name () throw () 
  { return "ray_tracer"; }

protected: // interface

  int rays_per_cell_;

private:
  int it, ia_x, ia_y, ia_z, ia_nx, ia_ny, ia_nz, ia_sx,
    ia_sy, ia_sz, ia_f, ia_r, ia_time;
  int ps;
  
};

#endif /* ENZO_ENZO_METHOD_RAYTRACER_HPP */
