// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Refresh.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2017-08-11
/// @brief    

#include "problem.hpp"

//----------------------------------------------------------------------

void Refresh::add_field(std::string field_name)
{
  const int id_field = cello::field_descr()->field_id(field_name);
  add_field(id_field);
}
//----------------------------------------------------------------------

void Refresh::add_field_src_dst(std::string field_src, std::string field_dst)
{
  const int id_field_src = cello::field_descr()->field_id(field_src);
  const int id_field_dst = cello::field_descr()->field_id(field_dst);
  add_field_src_dst(id_field_src,id_field_dst);
}

//----------------------------------------------------------------------

void Refresh::add_all_fields(std::string field_group)
{
  if (field_group == "") {
    all_fields_ = true;
  } else {
    Grouping * groups = cello::field_groups();
    int n = groups->size(field_group);
    for (int i=0; i<n; i++) {
      std::string field = groups->item(field_group,i);
      add_field(field);
    }
  }
}

//----------------------------------------------------------------------

int Refresh::data_size () const
{
  int count = 0;

  // WARNING: Skipping many fields since data methods are only called
  // when the Refresh object is a member of FieldFace, which in turn
  // only accesses field and particle lists and accumulate_

  SIZE_INT_ARRAY(&count,field_list_src_);
  SIZE_INT_ARRAY(&count,field_list_dst_);
  SIZE_INT_ARRAY(&count,particle_list_);
  
  SIZE_INT(&count,all_fields_);
  SIZE_INT(&count,all_particles_);
  SIZE_INT(&count,all_fluxes_);
  SIZE_INT(&count,accumulate_);

  return count;

}

//----------------------------------------------------------------------

char * Refresh::save_data (char * buffer) const
{
  char * p = buffer;

  SAVE_INT_ARRAY(&p,field_list_src_);
  SAVE_INT_ARRAY(&p,field_list_dst_);
  SAVE_INT_ARRAY(&p,particle_list_);
  
  SAVE_INT(&p,all_fields_);
  SAVE_INT(&p,all_particles_);
  SAVE_INT(&p,all_fluxes_);
  SAVE_INT(&p,accumulate_);

  ASSERT2 ("Refresh::save_data\n",
 	   "Actual size %ld does not equal computed size %d",
	   p-buffer,data_size(),
	   ((p-buffer)==data_size()));

  return p;
  
  
}

//----------------------------------------------------------------------

char * Refresh::load_data (char * buffer)
{
  char * p = buffer;

  LOAD_INT_ARRAY(&p,field_list_src_);
  LOAD_INT_ARRAY(&p,field_list_dst_);
  LOAD_INT_ARRAY(&p,particle_list_);

  LOAD_INT(&p,all_fields_);
  LOAD_INT(&p,all_particles_);
  LOAD_INT(&p,all_fluxes_);
  LOAD_INT(&p,accumulate_);

  ASSERT2 ("Refresh::load_data\n",
	   "Actual size %ld does not equal computed size %d",
	   p-buffer,data_size(),
	   ((p-buffer)==data_size()));

  return p;
}

