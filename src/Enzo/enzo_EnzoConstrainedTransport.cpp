// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoConstrainedTransport.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Mon May 6 2019
/// @brief    [\ref Enzo] Implementation of EnzoConstrainedTransport

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

void EnzoConstrainedTransport::compute_center_efield
(Block *block, int dim, std::string center_efield_name, Grouping &prim_group,
 int stale_depth)
{
  // Load the E-field
  EnzoFieldArrayFactory array_factory(block,stale_depth);
  EFlt3DArray efield = array_factory.from_name(center_efield_name);

  EnzoPermutedCoordinates coord(dim);
  int j = coord.j_axis();
  int k = coord.k_axis();

  // Load the jth and kth components of the velocity and cell-centered bfield
  EFlt3DArray velocity_j, velocity_k, bfield_j, bfield_k;
  velocity_j = array_factory.from_grouping(prim_group, "velocity", j);
  velocity_k = array_factory.from_grouping(prim_group, "velocity", k);
  bfield_j = array_factory.from_grouping(prim_group, "bfield", j);
  bfield_k = array_factory.from_grouping(prim_group, "bfield", k);

  for (int iz=0; iz<efield.shape(0); iz++) {
    for (int iy=0; iy<efield.shape(1); iy++) {
      for (int ix=0; ix<efield.shape(2); ix++) {
	efield(iz,iy,ix) = (-velocity_j(iz,iy,ix) * bfield_k(iz,iy,ix) +
			    velocity_k(iz,iy,ix) * bfield_j(iz,iy,ix));
      }
    }
  }
}

//----------------------------------------------------------------------

// The following is a helper function that actually computes the component of
// the edge-centered E-field along dimension i.
//
// First, a note on weight arrays:
//   If the i-direction poins along z, we would need to compute
//   dEzdx(ix-1/4,iy-1/2)
//
//   if upwind in positive y direction
//     dEdx(ix-1/4,iy-1/2) = 2.(E(ix,iy-1) - E(ix-1/2,iy-1))/dx
//   if upwind in negative y direction
//     dEdx(ix-1/4,iy-1/2) = 2.(E(ix,iy) - E(ix-1/2,iy))/dx
//   otherwise:
//     dEdx(ix-1/4,iy-1/2) = [(E(ix,iy-1) - E(ix-1/2,iy-1))/dx
//                            + (E(ix,iy) - E(ix-1/2,iy))/dx ]
//   We pulled out a factor of 2 and dx from the derivatives (they cancel)
//  
//   If upwind is in the positive y direction W_y = 1, if downwind W_y = 0,
//     and otherwise W_y = 0.5 (W_y is centered on faces along y-direction)
//
// The form of equations from Stone & Gardiner 09 are as follows. This is
// all necessary to compute Ez(ix-1/2, iy-1/2, iz): 
//  dEdj_r = dEdx(ix-1/4,iy-1/2) =
//       W_y(    ix,iy-1/2)  * (E(    ix,  iy-1) - E(ix-1/2,  iy-1)) +
//    (1-W_y(    ix,iy-1/2)) * (E(    ix,    iy) - E(ix-1/2,    iy))
//  dEdj_l = dEdx(ix-3/4,iy-1/2) =
//       W_y(  ix-1,iy-1/2)  * (E(ix-1/2,  iy-1) - E(  ix-1,  iy-1)) +
//    (1-W_y(  ix-1,iy-1/2)) * (E(ix-1/2,    iy) - E(  ix-1,    iy))
//  dEdk_r = dEdy(ix-1/2,iy-1/4) =
//       W_x(ix-1/2,   iy)  * (E(   ix-1,    iy) - E(  ix-1,iy-1/2)) +
//    (1-W_x(ix-1/2,   iy)) * (E(     ix,    iy) - E(    ix,iy-1/2))
//  dEdk_l = dEdy(ix-1/2,iy-3/4) =
//       W_x(ix-1/2, iy-1)  * (E(   ix-1,iy-1/2) - E(  ix-1,  iy-1)) +
//    (1-W_x(ix-1/2, iy-1)) * (E(     ix,iy-1/2) - E(    ix,  iy-1))
//
// This solution to the problem is easier to program if we reframe it as the
// calculation of Ez(ix+1/2, iy+1/2, iz). The above equations are then
// rewritten as:
//  dEdj_r = dEdx(ix+3/4,iy+1/2) =
//       W_y(  ix+1,iy+1/2)  * (E(  ix+1,    iy) - E(ix+1/2,    iy)) +
//    (1-W_y(  ix+1,iy+1/2)) * (E(  ix+1,  iy+1) - E(ix+1/2,  iy+1))
//  dEdj_l = dEdx(ix+1/4,iy+1/2) =
//       W_y(    ix,iy+1/2)  * (E(ix+1/2,    iy) - E(    ix,    iy)) +
//    (1-W_y(    ix,iy+1/2)) * (E(ix+1/2,  iy+1) - E(    ix,  iy+1))
//  dEdk_r = dEdy(ix+1/2,iy+3/4) =
//       W_x(ix+1/2,  iy+1)  * (E(     ix,  iy+1) - E(    ix,iy+1/2)) +
//    (1-W_x(ix+1/2,  iy+1)) * (E(   ix+1,  iy+1) - E(  ix+1,iy+1/2))
//  dEdk_l = dEdy(ix+1/2,iy+1/4) =
//       W_x(ix+1/2,   iy)  * (E(     ix,iy+1/2) - E(    ix,    iy)) +
//    (1-W_x(ix+1/2,   iy)) * (E(   ix+1,iy+1/2) - E(  ix+1,    iy))
//
// Now let's rewrite in totally Generalized notation:  [z->i, x->j, y->k] 
//  dEdj_r = dEdj(k+1/2,j+3/4,i) =
//       W_k(k+1/2,  j+1)  * (E(    k,  j+1) - E(    k,j+1/2)) +
//    (1-W_k(k+1/2,  j+1)) * (E(  k+1,  j+1) - E(  k+1,j+1/2))
//  dEdj_l = dEdj(k+1/2,j+1/4,i) =
//       W_k(k+1/2,    j)  * (E(    k,j+1/2) - E(    k,    j)) +
//    (1-W_k(k+1/2,    j)) * (E(  k+1,j+1/2) - E(  k+1,    j))
//  dEdk_r = dEdk(k+3/4,j+1/2,i) =
//       W_j(  k+1,j+1/2)  * (E(  k+1,    j) - E(k+1/2,    j)) +
//    (1-W_j(  k+1,j+1/2)) * (E(  k+1,  j+1) - E(k+1/2,  j+1))
//  dEdk_l = dEdk(k+1/4,j+1/2,i) =
//       W_j(    k,j+1/2)  * (E(k+1/2,     j) - E(    k,    j)) +
//    (1-W_j(    k,j+1/2)) * (E(k+1/2,   j+1) - E(    k,  j+1))
//
//
// Define subarrays of W_j, W_k, E_cen, E_j, and E_k
//   Weights:
//     Wj(k,j,i)       =   W_j(    k,j+1/2,i)
//     Wj_kp1(k,j,i)   =   W_j(  k+1,j+1/2,i)
//     Wk(k,j,i)       =   W_k(k+1/2,    j,i)
//     Wk_jp1(k,j,i)   =   W_k(k+1/2,  j+1,i)
//
//   Cell-Centered Efield (i-component):
//     Ec(k,j,i)       =     E(    k,    j,i)
//     Ec_jp1(k,j,i)   =     E(    k,  j+1,i)
//     Ec_kp1(k,j,i)   =     E(  k+1,    j,i)
//     Ec_jkp1(k,j,i)  =     E(  k+1,  j+1,i)
//
//   face-centered E-field (i-component):
//     Ej(k,j,i)       =     E(    k,j+1/2,i)
//     Ej_kp1(k,j,i)   =     E(  k+1,j+1/2,i)
//     Ek(k,j,i)       =     E(k+1/2,    j,i)
//     Ek_jp1(k,j,i)   =     E(k+1/2,  j+1,i)
//
//   edge-centered E-field (i-component):
//     Eedge(k,j,i)    =     E(k+1/2,j+1/2,i)
// To be independent of reconstruction method, compute E-field at all edges
// that lie within the mesh.

void compute_edge_(int xstart, int ystart, int zstart,
		   int xstop, int ystop, int zstop,
		   EFlt3DArray &Eedge, EFlt3DArray &Wj, EFlt3DArray &Wj_kp1,
		   EFlt3DArray &Wk, EFlt3DArray &Wk_jp1, EFlt3DArray &Ec,
		   EFlt3DArray &Ec_jkp1, EFlt3DArray &Ec_jp1,
		   EFlt3DArray &Ec_kp1, EFlt3DArray &Ej, EFlt3DArray &Ej_kp1,
		   EFlt3DArray &Ek, EFlt3DArray &Ek_jp1)
{
  for (int iz = zstart; iz < zstop; iz++){
    for (int iy = ystart; iy < ystop; iy++){
      for (int ix = xstart; ix < xstop; ix++){

	enzo_float dEdj_l, dEdj_r, dEdk_l, dEdk_r;

	//  dEdj(k+1/2,j+3/4,i) =
	//       W_k(k+1/2,  j+1)  * (E(    k,  j+1) - E(    k,j+1/2)) +
	//    (1-W_k(k+1/2,  j+1)) * (E(  k+1,  j+1) - E(  k+1,j+1/2))
	dEdj_r =
	       Wk_jp1(iz,iy,ix)  * ( Ec_jp1(iz,iy,ix) -     Ej(iz,iy,ix)) +
	  (1 - Wk_jp1(iz,iy,ix)) * (Ec_jkp1(iz,iy,ix) - Ej_kp1(iz,iy,ix));

	//  dEdj(k+1/2,j+1/4,i) =
	//       W_k(k+1/2,    j)  * (E(    k,j+1/2) - E(    k,    j)) +
	//    (1-W_k(k+1/2,    j)) * (E(  k+1,j+1/2) - E(  k+1,    j))
	dEdj_l =
	           Wk(iz,iy,ix)  * (     Ej(iz,iy,ix) -     Ec(iz,iy,ix)) +
	  (1 -     Wk(iz,iy,ix)) * ( Ej_kp1(iz,iy,ix) - Ec_kp1(iz,iy,ix));

	//  dEdk(k+3/4,j+1/2,i) =
	//       W_j(  k+1,j+1/2)  * (E(  k+1,    j) - E(k+1/2,    j)) +
	//    (1-W_j(  k+1,j+1/2)) * (E(  k+1,  j+1) - E(k+1/2,  j+1))
	dEdk_r =
	       Wj_kp1(iz,iy,ix)  * ( Ec_kp1(iz,iy,ix) -     Ek(iz,iy,ix)) +
	  (1 - Wj_kp1(iz,iy,ix)) * (Ec_jkp1(iz,iy,ix) - Ek_jp1(iz,iy,ix));

	//  dEdk(k+1/4,j+1/2,i) =
	//       W_j(    k,j+1/2)  * (E(k+1/2,     j) - E(    k,    j)) +
	//    (1-W_j(    k,j+1/2)) * (E(k+1/2,   j+1) - E(    k,  j+1))
	dEdk_l =
	           Wj(iz,iy,ix)  * (     Ek(iz,iy,ix) -     Ec(iz,iy,ix)) +
	  (1 -     Wj(iz,iy,ix)) * ( Ek_jp1(iz,iy,ix) - Ec_jp1(iz,iy,ix));

	Eedge(iz,iy,ix) = 0.25*(Ej(iz,iy,ix) + Ej_kp1(iz,iy,ix) +
				Ek(iz,iy,ix) + Ek_jp1(iz,iy,ix) +
				(dEdj_l-dEdj_r) + (dEdk_l - dEdk_r));

      }
    }
  }
}

//----------------------------------------------------------------------

void inplace_entry_multiply_(EFlt3DArray &array, enzo_float val){
  for (int iz = 0; iz < array.shape(0); iz++){
    for (int iy = 0; iy < array.shape(1); iy++){
      for (int ix = 0; ix < array.shape(2); ix++){
	array(iz,iy,ix)*=val;
      }
    }
  }
}

//----------------------------------------------------------------------

// Computes the edge-centered E-fields pointing in the ith direction
// It uses the component of the cell-centered E-field pointing in that
// direction, and the face-centered E-field pointed in that direction
// the face-centered E-fields are given by elements of jflux_ids and
// kflux_ids. dim points along i.
// i, j, and k are any cyclic permutation of x, y, z
//
// This Method is applicable for:
//    - 2D array (z-axis only has 1 entry), with dim = 2
// 
//
// The Athena++ code calculates a quantity they refer to as v_over_c at all
// cell-faces.
//   - Basically, this tells them how much to weight derivatives while
//     computing the edge_efield. If deemed necessary, these values can be
//     passed as weight_group
//   - Weight_group includes 3 temporary fields centered on the faces looking
//     down the x, y, and z direction. Currently, it expects values of 1 and 0
//     to indicate that the upwind direction is in positive and negative
//     direction, or 0 to indicate no upwind direction.
void EnzoConstrainedTransport::compute_edge_efield
(Block *block, int dim, std::string center_efield_name, Grouping &efield_group,
 Grouping &jflux_group, Grouping &kflux_group, Grouping &weight_group,
 int stale_depth)
{

  EnzoPermutedCoordinates coord(dim);
  // determine components of j and k unit vectors:
  int j_x, j_y, j_z, k_x, k_y, k_z;
  coord.j_unit_vector(j_z, j_y, j_x);
  coord.k_unit_vector(k_z, k_y, k_x);

  // Initialize Cell-Centered E-fields
  EFlt3DArray Ec, Ec_jp1, Ec_kp1, Ec_jkp1;
  EnzoFieldArrayFactory array_factory(block,stale_depth);
  Ec = array_factory.from_name(center_efield_name);
  Ec_jp1  = coord.left_edge_offset(Ec, 0, 1, 0);
  Ec_kp1  = coord.left_edge_offset(Ec, 1, 0, 0);
  Ec_jkp1 = coord.left_edge_offset(Ec, 1, 1, 0);

  // Initialize edge-centered Efield [it maps (k,j,i) -> (k+1/2,j+1/2,i)]
  EFlt3DArray Eedge = array_factory.from_grouping(efield_group, "efield", dim);

  // Initialize face-centered E-fields
  EFlt3DArray Ej, Ej_kp1, Ek, Ek_jp1;

  // Ex(k,j+1/2,i) = -1.*yflux(Bz)
  Ej = array_factory.from_grouping(jflux_group, "bfield", coord.k_axis());
  inplace_entry_multiply_(Ej,-1.);
  Ej_kp1 = coord.left_edge_offset(Ej, 1, 0, 0);

  // Ex(k+1/2,j,i) = zflux(By)
  Ek = array_factory.from_grouping(kflux_group, "bfield", coord.j_axis());
  // No need to multiply entries by minus 1
  Ek_jp1 = coord.left_edge_offset(Ek, 0, 1, 0);

  // Initialize the weight arrays
  EFlt3DArray Wj, Wj_kp1, Wk, Wk_jp1;
  Wj = array_factory.from_grouping(weight_group, "weight", coord.j_axis());
  Wj_kp1 = coord.left_edge_offset(Wj, 1, 0, 0);
  Wk = array_factory.from_grouping(weight_group, "weight", coord.k_axis());
  Wk_jp1 = coord.left_edge_offset(Wk, 0, 1, 0);

  // Integration limits
  //
  // If computing the edge E-field along z-direction:
  //    - If grid is 3D (the grid has more than 1 cell along the z-component),
  //      there is no need to compute the e-field at iz=0 or iz= mz-1 since
  //      we have E-fields on the exterior of the mesh. (This logic applies to
  //      components of other dimensions). Generalizing to computing
  //      E-field along dimension i, then we need to compute it at i=1 up to
  //      (but not including imax-1)
  // For all cases, if we are computing the E-field along dimension i, then we
  // need to compute it at:  j = 1/2 up to (but not including) j = jmax-1/2
  //                         k = 1/2 up to (but not including) k = kmax-1/2
  //
  // Note if an quantitiy is face-centered along dimension dim:
  //    idim maps to idim+1/2
  //
  // To summarize:
  //    istart = 1     istop = imax - 1
  //    jstart = 0     jstop = jmax - 1
  //    kstart = 0     kstop = kmax - 1

  int xstart = 1 - j_x - k_x; // if dim==0: 1, otherwise: 0
  int ystart = 1 - j_y - k_y; // if dim==1: 1, otherwise: 0
  int zstart = 1 - j_z - k_z; // if dim==2: 1, otherwise: 0

  compute_edge_(xstart, ystart, zstart,
		Ec.shape(2) - 1, Ec.shape(1) - 1, Ec.shape(0) - 1,
		Eedge, Wj, Wj_kp1, Wk, Wk_jp1, Ec, Ec_jkp1, Ec_jp1, Ec_kp1,
		Ej, Ej_kp1, Ek, Ek_jp1);
}

//----------------------------------------------------------------------

void EnzoConstrainedTransport::compute_all_edge_efields
  (Block *block, Grouping &prim_group, Grouping &xflux_group,
   Grouping &yflux_group, Grouping &zflux_group, std::string center_efield_name,
   Grouping &efield_group, Grouping &weight_group, int stale_depth)
{

  for (int i = 0; i < 3; i++){
    compute_center_efield(block, i, center_efield_name, prim_group,
			  stale_depth);
    Grouping *jflux_group;
    Grouping *kflux_group;
    if (i == 0){
      jflux_group = &yflux_group;
      kflux_group = &zflux_group;
    } else if (i==1){
      jflux_group = &zflux_group;
      kflux_group = &xflux_group;
    } else {
      jflux_group = &xflux_group;
      kflux_group = &yflux_group;
    }

    compute_edge_efield(block, i, center_efield_name, efield_group,
			*jflux_group, *kflux_group, weight_group, stale_depth);
  }
}

//----------------------------------------------------------------------

// Compute the face-centered B-field component along the ith dimension
//
// Bnew_i(k, j, i-1/2) = Bold_i(k, j, i-1/2) -
//     dt/dj*(E_k(    k,j+1/2,i-1/2) - E_k(    k,j-1/2,i-1/2) +
//     dt/dk*(E_j(k+1/2,    j,i-1/2) - E_j(k-1/2,    j,i-1/2)
// [The positioning of dt/dj with respect to E_k is correct]
//
// Bnew_i(k, j, i+1/2) =
//   Bold_i(k, j, i+1/2) - E_k_term(k,j,i+1/2) + E_j_term(k,j,i+1/2)
//
// Assuming 3D:
//   E_k_term(k,j,i+1/2) = dt/dj*(E_k(k,j+1/2,i+1/2) - E_k(k,j-1/2,i+1/2))
//   E_j_term(k,j,i+1/2) = dt/dk*(E_j(k+1/2,j,i+1/2) - E_j(k-1/2,j,i+1/2))
// 
// We define: (notation DIFFERENT from compute_edge_efield & compute_edge_)
//   ej_Lk(k,j,i) = E_j(k-1/2,     j, i+1/2) 
//   ej_Rk(k,j,i) = E_j(k+1/2,     j, i+1/2)
//   ek_Lj(k,j,i) = E_k(    k, j-1/2, i+1/2) 
//   ek_Rj(k,j,i) = E_k(    k, j+1/2, i+1/2)
//
// Then:
//   E_k_term(k,j,i+1/2) = dt/dj*(ek_Rj(k,j,i) - ek_Lj(k,j,i))
//   E_j_term(k,j,i+1/2) = dt/dk*(ej_Rk(k,j,i) - ej_Lk(k,j,i))
void EnzoConstrainedTransport::update_bfield(Block *block, int dim,
					     Grouping &efield_group,
					     Grouping &cur_bfieldi_group,
					     Grouping &out_bfieldi_group,
					     enzo_float dt, int stale_depth)
{
  EnzoPermutedCoordinates coord(dim);

  // compute the ratios of dt to the widths of cells along j and k directions
  EnzoBlock *enzo_block = enzo::block(block);
  enzo_float dtdj = dt/enzo_block->CellWidth[coord.j_axis()];
  enzo_float dtdk = dt/enzo_block->CellWidth[coord.k_axis()];

  // The following comments all assume that we are talking about unstaled
  // region (and that we have dropped all staled cells)
  EnzoFieldArrayFactory array_factory(block,stale_depth);
  
  // Load edge centered efields 
  EFlt3DArray E_j, ej_Lk, ej_Rk, E_k, ek_Lj, ek_Rj;
  E_j = array_factory.from_grouping(efield_group, "efield", coord.j_axis());
  E_k = array_factory.from_grouping(efield_group, "efield", coord.k_axis());

  // Load interface bfields field (includes exterior faces)
  EFlt3DArray cur_bfield, bcur, out_bfield, bout;
  cur_bfield = array_factory.from_grouping(cur_bfieldi_group, "bfield", dim);
  out_bfield = array_factory.from_grouping(out_bfieldi_group, "bfield", dim);

  // Now to take slices. If the unstaled region of the grid has shape has shape
  // (mk, mj, mi) then:
  //   - E_j has shape (mk-1, mj, mi-1)
  //       ej_Lk includes k=1/2 up to (but not including) k=mk-3/2
  //       ej_Rk includes k=3/2 up to (but not including) k=mk-1/2
  //   - E_k has shape (mk, mj-1, mi-1)
  //       ek_Lj includes j=1/2 up to (but not including) j=mj-3/2
  //       ek_Rj includes j=3/2 up to (but not including) j=mj-1/2
  //   - cur_bfield and out_bfield each have shape (mk, mj, mi+1)
  //       bnew and bout only include interior faces
  //
  // Also need to omit outermost layer of cell-centered vals
  CSlice full_ax(nullptr, nullptr); // includes full axis
  CSlice inner_cent(1,-1);          // excludes outermost cell-centered values

  // the following arrays should all have the same shape
  ej_Lk = coord.get_subarray(E_j, CSlice(0,      -1), inner_cent, full_ax);
  ej_Rk = coord.get_subarray(E_j, CSlice(1, nullptr), inner_cent, full_ax);
  ek_Lj = coord.get_subarray(E_k, inner_cent, CSlice(0,      -1), full_ax);
  ek_Rj = coord.get_subarray(E_k, inner_cent, CSlice(1, nullptr), full_ax);
  bcur = coord.get_subarray(cur_bfield, inner_cent, inner_cent, CSlice(1,-1));
  bout = coord.get_subarray(out_bfield, inner_cent, inner_cent, CSlice(1,-1));

  // We could simplify this iteration by using subarrays - However, it would be
  // more complicated
  for (int iz=0; iz<bout.shape(0); iz++) {
    for (int iy=0; iy<bout.shape(1); iy++) {
      for (int ix=0; ix<bout.shape(2); ix++) {

	// E_k_term(k,j,i+1/2) = dt/dj*(E_k(k,j+1/2,i+1/2) - E_k(k,j-1/2,i+1/2))
	enzo_float E_k_term = dtdj*(ek_Rj(iz,iy,ix) - ek_Lj(iz,iy,ix));

	// E_j_term(k,j,i+1/2) = dt/dk*(E_j(k+1/2,j,i+1/2) - E_j(k-1/2,j,i+1/2))
	enzo_float E_j_term = dtdk*(ej_Rk(iz,iy,ix) - ej_Lk(iz,iy,ix)); 

	// Bnew_i(k, j, i+1/2) =
	//   Bold_i(k, j, i+1/2) - E_k_term(k,j,i+1/2) + E_j_term(k,j,i+1/2)
	bout(iz,iy,ix) = bcur(iz,iy,ix) - E_k_term + E_j_term;
      }
    }
  }
}

//----------------------------------------------------------------------

// This method also intentionally includes calculation of bfields in the
// outermost cells so that it can be used to initially setup the bfield.
//
// Compute cell-centered bfield along dimension i
//   B_i(k,j,i) = 0.5*(B_i(k,j,i+1/2) + B_i(k,j,i-1/2))
// For a simpler implementation, we will rewrite this as:
//   B_i(k,j,i+1) = 0.5*(B_i(k,j,i+3/2) + B_i(k,j,i+1/2))
// We define:
//   B_center(k,j,i)   ->  B_i(k,j,i+1)
//   Bi_left(k,j,i)    ->  B_i(k,j,i+1/2)
//   Bi_right(k,j,i)   ->  B_i(k,j,i+3/2)
void EnzoConstrainedTransport::compute_center_bfield(Block *block, int dim,
						     Grouping &bfieldc_group,
						     Grouping &bfieldi_group,
						     int stale_depth)
{
  EnzoPermutedCoordinates coord(dim);
  EnzoFieldArrayFactory array_factory(block,stale_depth);

  // Load cell-centerd field
  EFlt3DArray b_center = array_factory.from_grouping(bfieldc_group, "bfield",
						     coord.i_axis());
  // Load Face-centered fields
  EFlt3DArray bi_left = array_factory.from_grouping(bfieldi_group, "bfield",
						    coord.i_axis());
  // Get the view of the Face-center field that starting from i=1
  EFlt3DArray bi_right = coord.left_edge_offset(bi_left,0,0,1);

  // iteration limits are compatible with a 2D grid and 3D grid
  for (int iz=0; iz<b_center.shape(0); iz++) {
    for (int iy=0; iy<b_center.shape(1); iy++) {
      for (int ix=0; ix<b_center.shape(2); ix++) {
	b_center(iz,iy,ix) = 0.5*(bi_left(iz,iy,ix) + bi_right(iz,iy,ix));
      }
    }
  }
}
