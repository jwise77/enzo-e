// See LICENSE_CELLO file for license and copyright information
/// @file     enzo_EnzoInitialInclinedWave.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Fri May 24 2019
/// @brief    [\ref Enzo] Initialization routine for spherical cloud in a wind

#ifndef ENZO_ENZO_INITIAL_CLOUD_HPP
#define ENZO_ENZO_INITIAL_CLOUD_HPP

// tolerance to be used when checking the equivalence of 2 floating point
// numbers in EnzoInitialCloud. This is the same as ETA_TOLERANCE for doubles.
#define INIT_CLOUD_TOLERANCE 1.0e-10

// represents a spherical region
class SphereRegion;

class EnzoInitialCloud : public Initial {
  /// @class    EnzoInitialCloud
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Initializer a cloud crushing problem
  /// The wind is assumed to blow along the x-directions

public: // interface

  /// Constructor
  EnzoInitialCloud(int cycle, double time, ParameterGroup p) noexcept
    : Initial(cycle,time),
      subsample_n_(p.value_integer("subsample_n",0)),
      cloud_radius_(p.value_float("cloud_radius",0.0)),
      cloud_center_x_(p.value_float("cloud_center_x",0.0)),
      cloud_center_y_(p.value_float("cloud_center_y",0.0)),
      cloud_center_z_(p.value_float("cloud_center_z",0.0)),
      density_cloud_(p.value_float("cloud_density",0.0)),
      density_wind_(p.value_float("wind_density",0.0)),
      etot_wind_(p.value_float("wind_total_energy",0.0)),
      eint_wind_(p.value_float("wind_internal_energy",0.0)),
      velocity_wind_(p.value_float("wind_velocity",0.0)),
      metal_mass_frac_(p.value_float("metal_mass_fraction",0.0)),
      perturb_Nwaves_(p.value_integer("perturb_Nwaves", 0)),
      perturb_amplitude_(p.value_float("perturb_amplitude", 0.)),
      perturb_min_wavelength_(p.value_float
                              ("perturb_min_lambda",
                               std::numeric_limits<double>::min())),
      perturb_max_wavelength_(p.value_float
                              ("perturb_max_lambda",
                               std::numeric_limits<double>::min())),
      perturb_seed_(0)
  {
    if ( (perturb_Nwaves_ > 0) &&
         (perturb_max_wavelength_ == std::numeric_limits<double>::min()) ){
      perturb_max_wavelength_ = cloud_radius_;
    }

    int perturb_seed_tmp = p.value_integer("perturb_seed",0);
    ASSERT("EnzoInitialCloud",
           "Initial:cloud:perturb_seed must be a 32-bit unsigned integer",
           (perturb_seed_tmp >= 0) && (perturb_seed_tmp <= 4294967295L));
    perturb_seed_ = (unsigned int) perturb_seed_tmp;

    int uniform_bfield_length = p.list_length("uniform_bfield");
    if (uniform_bfield_length == 0) {
      initialize_uniform_bfield_ = false;
    } else if (uniform_bfield_length == 3) {
      initialize_uniform_bfield_ = true;
      for (int i = 0; i <3; i++) {
        uniform_bfield_[i] = p.list_value_float(i,"uniform_bfield");
      }
    } else {
      ERROR("EnzoInitialCloud",
            "Initial:cloud:uniform_bfield must contain 0 or 3 entries.");
    }

    ASSERT("EnzoInitialCloud", "subsample_n must be >=0", subsample_n_>=0);
    ASSERT("EnzoInitialCloud", "cloud_radius must be positive",
	   cloud_radius_>0);
    ASSERT("EnzoInitialCloud", "density_wind must be positive",
           density_wind_>0);
    ASSERT("EnzoInitialCloud", "density_cloud must be positive",
           density_cloud_>0);
    double ke_w = 0.5*velocity_wind_*velocity_wind_;
    ASSERT("EnzoInitialCloud", "etot_wind must exceed wind kinetic energy",
	   etot_wind_>ke_w);
    if (initialize_uniform_bfield_) {
      double me_w = 0.5*(uniform_bfield_[0]*uniform_bfield_[0]+
			 uniform_bfield_[1]*uniform_bfield_[1]+
			 uniform_bfield_[2]*uniform_bfield_[2])/density_wind_;
      ASSERT("EnzoInitialCloud",
	     "etot_wind must exceed wind specific kinetic and magnetic energy",
	     etot_wind_>(ke_w+me_w));
    }
    ASSERT("EnzoInitialCloud", "eint_wind must be zero or positive.",
	   eint_wind_>=0.);
    ASSERT("EnzoInitialCloud", "metal_mass_frac must be in [0,1]",
	   metal_mass_frac_ >=0 && metal_mass_frac_ <=1);
    if (perturb_Nwaves_ > 0){
      ASSERT("EnzoInitialCloud",
             "perturb_amplitude_ must be positive if perturb_Nwaves_>0",
             perturb_amplitude_ > 0);
      ASSERT("EnzoInitialCloud",
             "perturb_min_wavelength_ must be positive if perturb_Nwaves_>0",
             perturb_min_wavelength_ > 0);
      ASSERT("EnzoInitialCloud",
             "perturb_min_wavelength_ must be positive if perturb_Nwaves_>0",
             perturb_max_wavelength_ > perturb_min_wavelength_);
    }
  }

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialCloud);

  /// CHARM++ migration constructor
  EnzoInitialCloud(CkMigrateMessage *m)
    : Initial(m),
      subsample_n_(0),
      cloud_radius_(0.),
      cloud_center_x_(0.),
      cloud_center_y_(0.),
      cloud_center_z_(0.),
      density_cloud_(0.),
      density_wind_(0.),
      etot_wind_(0.),
      velocity_wind_(0.),
      metal_mass_frac_(0.),
      initialize_uniform_bfield_(false),
      perturb_Nwaves_(0),
      perturb_amplitude_(),
      perturb_min_wavelength_(0.),
      perturb_max_wavelength_(0.),
      perturb_seed_(0)
  {
    for (int i = 0; i < 3; i++){ uniform_bfield_[i] = 0.; }
  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    // NOTE: update whenever attributes change

    TRACEPUP;

    Initial::pup(p);
    p | subsample_n_;
    p | cloud_radius_;
    p | cloud_center_x_;
    p | cloud_center_y_;
    p | cloud_center_z_;
    p | density_cloud_;
    p | density_wind_;
    p | etot_wind_;
    p | velocity_wind_;
    p | metal_mass_frac_;
    p | initialize_uniform_bfield_;
    PUParray(p, uniform_bfield_, 3);
    p | perturb_Nwaves_;
    p | perturb_amplitude_;
    p | perturb_min_wavelength_;
    p | perturb_min_wavelength_;
    p | perturb_seed_;
  }

public: // virtual methods

  /// Initialize the block
  virtual void enforce_block
  ( Block * block, const Hierarchy * hierarchy ) throw();

private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// number of subsampled cells per cell side is (int)pow(2,subsample_n_);
  int subsample_n_;
  /// Cloud radius in code units
  double cloud_radius_;
  /// Cloud center
  double cloud_center_x_;
  double cloud_center_y_;
  double cloud_center_z_;

  /// cloud and wind densities
  double density_cloud_;
  double density_wind_;

  /// specific total energy of the wind
  double etot_wind_;

  /// specific internal energy of the wind
  double eint_wind_;

  /// velocity of the wind
  double velocity_wind_;

  double metal_mass_frac_;

  /// indicates whether to initialize a uniform bfield
  bool initialize_uniform_bfield_;

  /// values of the uniform bfield
  double uniform_bfield_[3];


  /// The number of sinusoidal density waves to use for the perturbation.
  int perturb_Nwaves_;

  /// The amplitude of the perturbation wavelength. This is a fraction that is
  /// multiplied by the cloud density
  double perturb_amplitude_;

  /// The minimum wavelength that the perturbations can have
  double perturb_min_wavelength_;

  /// The maximum wavelength that the perturbations can have
  double perturb_max_wavelength_;

  /// The random seed used to seed the density perturbations
  unsigned int perturb_seed_;
};

#endif //ENZO_ENZO_INITIAL_CLOUD_HPP
