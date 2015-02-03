
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
   symm_func_list.push_back(particle_hole);
   
   se_checked_tensor checked(FREQ_COUNT_SE); 

   for(int w = 0; w < FREQ_COUNT_SE; ++w)
   {
      index_1p_t ind(w);
      if ( !checked(ind) )  							// if tensor object not yet related to any other
      {
	 checked(ind) = true; 
	 ind_cpl_list.push_back(ind);							// push standard representative into ind_cpl_list
	 int ind_cpl_list_pos =  ind_cpl_list.size() - 1 ; 					// position in ind_cpl_list is size - 1
	 (*se_ptr)(ind) = ind_cpl_t( ind_cpl_list_pos ); 				// save position with trivial operations
	 operation track_op(false,false);
	 iterate( ind, track_op, *se_ptr, checked, symm_func_list, ind_cpl_list_pos ); 		// start iterating on index 
      }
   }


}

void iterate( const index_1p_t& ind, const operation& track_op, se_tensor& selfEn, se_checked_tensor& checked, vector<symm_func_1p_t> symm_func_list , int ind_cpl_list_pos )
{
   for(auto symm_func: symm_func_list) 		// iterate over list of all symmetries specified
   {
      index_1p_t ind_it = ind;			// copy ind
      operation curr_op = symm_func(ind_it) * track_op;	// apply symmetry operation and track operations applied
      if( !checked(ind_it) )		// if resulting tensor index not yet related to any other
      { 
	 checked(ind_it) = true; 
	 selfEn(ind_it) = ind_cpl_t(ind_cpl_list_pos, curr_op); 		// relate to position
	 iterate( ind_it, curr_op, selfEn, checked, symm_func_list, ind_cpl_list_pos );	// iterate further	
      }
   }
}

//----- Antisymmetry

operation compl_conj(index_1p_t& ind)
{
   freq_sign_change(ind.w, FREQ_COUNT_SE);

   return operation(false,true);
}

operation particle_hole(index_1p_t& ind)
{
   return operation(false,true);
}

