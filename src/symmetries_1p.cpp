
/************************************************************************************************//**
 *  		
 * 	file: 		symmetries_1p.cpp
 * 	contents:   	see symmetries_1p.h
 * 
 ****************************************************************************************************/

#include <symmetries_1p.h>

using namespace std;

void init_symm( shared_ptr<se_tensor> se_ptr, vector<index_1p_t>& ind_cpl_list )
{
   // Initialize vector containing all symmetry functions 
   vector<symm_func_1p_t> symm_func_list; // = { exch_in, exch_out, compl_conj, time_rev, particle_hole, rot_k, mirror_vert, mirror_diag };

   symm_func_list.push_back(compl_conj);
   symm_func_list.push_back(time_rev);
   symm_func_list.push_back(particle_hole);

#ifndef NO_MOMENTA
   symm_func_list.push_back(rot_k);
   symm_func_list.push_back(mirror_vert);
#endif

   for(int w = 0; w < FREQ_COUNT_SE; ++w)
      for(int k = 0; k < PATCH_COUNT; ++k)
	 for(int s_in = 0; s_in < QN_COUNT; ++s_in)
	    for(int s_out = 0; s_out < QN_COUNT; ++s_out)
	    {
	       index_1p_t ind(w, k, s_in, s_out );
	       if ( !(*se_ptr)(ind).checked )  							// if tensor object not yet related to any other
	       {
		  ind_cpl_list.push_back(ind);							// push standard representative into ind_cpl_list
		  int ind_cpl_list_pos =  ind_cpl_list.size() - 1 ; 					// position in ind_cpl_list is size - 1
		  (*se_ptr)(ind) = ind_cpl_t( ind_cpl_list_pos ); 				// save position with trivial operations
		  operation track_op(false,false);
		  iterate( ind, track_op, *se_ptr, symm_func_list, ind_cpl_list_pos ); 		// start iterating on index 
	       }
	    }


}

void iterate( const index_1p_t& ind, const operation& track_op, se_tensor& vertex, vector<symm_func_1p_t> symm_func_list , int ind_cpl_list_pos )
{
   for(auto symm_func: symm_func_list) 		// iterate over list of all symmetries specified
   {
      index_1p_t ind_it = ind;			// copy ind
      operation curr_op = symm_func(ind_it) * track_op;	// apply symmetry operation and track operations applied
      if( !vertex(ind_it).checked )		// if resulting tensor index not yet related to any other
      { 
	 vertex(ind_it) = ind_cpl_t(ind_cpl_list_pos, curr_op); 		// relate to position
	 iterate( ind_it, curr_op, vertex, symm_func_list, ind_cpl_list_pos );	// iterate further	
      }
   }
}

//----- Antisymmetry

operation compl_conj(index_1p_t& ind)
{
   freq_sign_change(ind.w);

   // Changing momenta sign is unnecessary since equal to rotating twice by 90 degrees

   swap(ind.s_in, ind.s_out);
   return operation(false,true);
}

operation time_rev(index_1p_t& ind)
{
   swap(ind.s_in, ind.s_out);
   return operation(false,false);
}

operation particle_hole(index_1p_t& ind)
{
#ifndef NO_MOMENTA
   mirror_mom_pipi(ind.k);
#endif
   return operation(false,true);
}

operation rot_k(index_1p_t& ind)
{
   ind.k = rot_k_ind_arr[ind.k];
   return operation(false,false);
}

operation mirror_vert(index_1p_t& ind)
{
   mirror_mom_vert(ind.k);
   return operation(false,false);
}
