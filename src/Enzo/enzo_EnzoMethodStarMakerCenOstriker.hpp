// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodStarMakerCenOstriker.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Thu Apr  1 16:14:38 PDT 2010
/// @brief    [\ref Enzo] Implementation of Enzo PPM hydro method

#ifndef ENZO_ENZO_METHOD_STARMAKER_CO_HPP
#define ENZO_ENZO_METHOD_STARMAKER_CO_HPP

class EnzoMethodStarMakerCenOstriker : public Method {

  /// @class    EnzoMethodStarMakerCenOstriker
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate Enzo's PPM hydro method

public: // interface

  /// Create a new EnzoMethodStarMakerCenOstriker object
  EnzoMethodStarMakerCenOstriker(const EnzoConfig * enzo_config);

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodStarMakerCenOstriker);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodStarMakerCenOstriker (CkMigrateMessage *m)
    : comoving_coordinates_(false)
  {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  
  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "star_maker_co"; }

protected: // interface

  int comoving_coordinates_;
};

#endif /* ENZO_ENZO_METHOD_STARMAKER_CO_HPP */
