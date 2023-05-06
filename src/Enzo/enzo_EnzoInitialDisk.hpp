#ifndef ENZO_ENZO_INITIAL_DISK_HPP
#define ENZO_ENZO_INITIAL_DISK_HPP

class EnzoInitialDisk : public Initial {

private:

public: // interface

  /// CHARM++ constructor
  EnzoInitialDisk(const EnzoConfig * enzo_config) throw();

  PUPable_decl(EnzoInitialDisk);

  /// Charm++ PUP::able migration constructor
  EnzoInitialDisk (CkMigrateMessage *m)
    : Initial (m)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Initialize block
  virtual void enforce_block
  ( Block * block,
    const Hierarchy * hierarchy ) throw();

  void InitialGasDistribution(Block * block);

  /// Destructor
  virtual ~EnzoInitialDisk(void) throw();

private: // attributes

  double center_position_[3];
  double scale_length_;
  double scale_height_;
  double central_mass_;
  double disk_mass_;
  double toomre_parameter_;
  double disk_temperature_;
  double disk_metal_fraction_;
  double disk_midplane_density_;
  double disk_outer_radius_;
  double disk_inner_radius_;

  int dual_energy_;
  double uniform_density_;
  double gamma_;
  double mu_;
  double tiny_number_;

};

#endif /* ENZO_ENZO_INITIAL_ISOLATED_GALAXY_HPP */
