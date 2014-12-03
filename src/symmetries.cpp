
/************************************************************************************************//**
 *  		
 * 	file: 		symmetries.cpp
 * 	contents:   	see symmetries.h
 * 
 ****************************************************************************************************/

#include <symmetries.h>

using namespace std;

void init_symm( shared_ptr<vertex_tensor> vertex_ptr, vector<index_2p_t>& ind_cpl_list )
{
   // Initialize vector containing all symmetry functions 
   vector<symm_func_t> symm_func_list; // = { exch_in, exch_out, compl_conj, time_rev, particle_hole, rot_k, mirror_vert, mirror_diag };

   symm_func_list.push_back(exch_in);
   symm_func_list.push_back(exch_out);
   symm_func_list.push_back(compl_conj);
   symm_func_list.push_back(time_rev);
   symm_func_list.push_back(particle_hole);
   symm_func_list.push_back(rot_k);
   symm_func_list.push_back(mirror_vert);
   //symm_func_list.push_back(mirror_diag); // Unnecessary, generated by rot_k and mirror_vert


   for(int w1_in = 0; w1_in < FREQ_COUNT_VERT; ++w1_in)
      for(int w2_in = 0; w2_in < FREQ_COUNT_VERT; ++w2_in)
	 for(int w1_out = 0; w1_out < FREQ_COUNT_VERT; ++w1_out)
	    for(int k1_in = 0; k1_in < PATCH_COUNT; ++k1_in)
	       for(int k2_in = 0; k2_in < PATCH_COUNT; ++k2_in)
		  for(int k1_out = 0; k1_out < PATCH_COUNT; ++k1_out)
		     for(int s1_in = 0; s1_in < QN_COUNT; ++s1_in)
			for(int s2_in = 0; s2_in < QN_COUNT; ++s2_in)
			   for(int s1_out = 0; s1_out < QN_COUNT; ++s1_out)
			      for(int s2_out = 0; s2_out < QN_COUNT; ++s2_out)
			      {
				 index_2p_t ind(w1_in, w2_in, w1_out, k1_in, k2_in, k1_out, s1_in, s2_in, s1_out, s2_out);
				 if ( (*vertex_ptr)(ind).ind == -1 )  							// if tensor object not yet related to any other
				 {
				    ind_cpl_list.push_back(ind);							// push standard representative into ind_cpl_list
				    int ind_cpl_list_pos =  ind_cpl_list.size() - 1 ; 					// position in ind_cpl_list is size - 1
				    (*vertex_ptr)(ind) = ind_cpl_t( ind_cpl_list_pos ); 				// save position with trivial operations
				    operation track_op(false,false);
				    iterate( ind, track_op, *vertex_ptr, symm_func_list, ind_cpl_list_pos ); 		// start iterating on index 
				 }
			      }


}

void iterate( const index_2p_t& ind, const operation& track_op, vertex_tensor& vertex, vector<symm_func_t> symm_func_list , int ind_cpl_list_pos )
{
   for(auto symm_func: symm_func_list) 		// iterate over list of all symmetries specified
   {
      index_2p_t ind_it = ind;			// copy ind
      operation curr_op = symm_func(ind_it) * track_op;	// apply symmetry operation and track operations applied
      if( !vertex(ind_it).checked )		// if resulting tensor index not yet related to any other
      { 
	 vertex(ind_it) = ind_cpl_t(ind_cpl_list_pos, curr_op); 		// relate to position
	 iterate( ind_it, curr_op, vertex, symm_func_list, ind_cpl_list_pos );	// iterate further	
      }
   }
}

//----- Antisymmetry

operation exch_in(index_2p_t& ind)
{
   swap(ind.w1_in, ind.w2_in);
   swap(ind.k1_in, ind.k2_in);
   swap(ind.s1_in, ind.s2_in);
   return operation(true,false);
}

operation exch_out(index_2p_t& ind) 
{
   int w2_out = ind.w1_in + ind.w2_in - ind.w1_out; // calculate w2_out by means of frequency conservation
   if ( 0 <= w2_out && w2_out < FREQ_COUNT_VERT )	// check if w2_out inside the grid, otherwise skip symmetry
   {
      ind.w1_out = w2_out;
      ind.k1_out = sum_mom[ sum_mom[ind.k1_in][ind.k2_in] ] [ sign_change_k_ind_arr[ind.k1_out] ]; // calculate k2_out and assign to k1_out
      swap(ind.s1_out,ind.s2_out);
      return operation(true,false);
   }
   return operation(false,false); 
}

operation compl_conj(index_2p_t& ind)
{
   int w2_out = ind.w1_in + ind.w2_in - ind.w1_out; // calculate w2_out by means of frequency conservation
   if ( 0 <= w2_out && w2_out < FREQ_COUNT_VERT )	// check if w2_out inside the grid, otherwise skip symmetry
   {
      swap(ind.w1_in, ind.w1_out);
      ind.w2_in = w2_out;

      freq_sign_change(ind.w1_in);
      freq_sign_change(ind.w2_in);
      freq_sign_change(ind.w1_out);

      // Changing all momenta signs is unnecessary since equal to rotating twice by 90 degrees

      swap(ind.s1_in, ind.s2_in);
      swap(ind.s1_out, ind.s2_out);

      return operation(false,true);
   }
   return operation(false, false);
}

operation time_rev(index_2p_t& ind)
{
   int w2_out = ind.w1_in + ind.w2_in - ind.w1_out; // calculate w2_out by means of frequency conservation
   if ( 0 <= w2_out && w2_out < FREQ_COUNT_VERT )	// check if w2_out inside the grid, otherwise skip symmetry
   {
      swap(ind.w1_in, ind.w1_out);
      ind.w2_in = w2_out;
      
      // Changing all momenta signs is unnecessary since equal to rotating twice by 90 degrees

      int k2_out = sum_mom[ sum_mom[ind.k1_in][ind.k2_in] ] [ sign_change_k_ind_arr[ind.k1_out] ]; // calculate k2_out 
      swap(ind.k1_in, ind.k1_out); 	// swap k1_in and k1_out
      ind.k2_in = k2_out; 		// swap k2_in and k2_out

      swap(ind.s1_in, ind.s2_in);
      swap(ind.s1_out, ind.s2_out);

      // Note: no operation needed
   }
   return operation(false,false);
}

operation particle_hole(index_2p_t& ind)
{
   mirror_mom_pipi(ind.k1_in);
   mirror_mom_pipi(ind.k2_in);
   mirror_mom_pipi(ind.k1_out);
   return operation(false,true);
}

operation rot_k(index_2p_t& ind)
{
   ind.k1_in = rot_k_ind_arr[ind.k1_in];
   ind.k2_in = rot_k_ind_arr[ind.k2_in];
   ind.k1_out = rot_k_ind_arr[ind.k1_out];
   return operation(false,false);
}

operation mirror_vert(index_2p_t& ind)
{
   mirror_mom_vert(ind.k1_in);
   mirror_mom_vert(ind.k2_in);
   mirror_mom_vert(ind.k1_out);
   return operation(false,false);
}

operation mirror_diag(index_2p_t& ind)
{
   mirror_mom_diag(ind.k1_in);
   mirror_mom_diag(ind.k2_in);
   mirror_mom_diag(ind.k1_out);
   return operation(false,false);
}

void freq_sign_change(int& ind)
{
   ind = FREQ_COUNT_VERT - ind - 1;
}

void mom_sign_change(int& ind)
{
   ind = sign_change_k_ind_arr[ind];
}

void mirror_mom_vert(int& ind)
{
   ind = mirror_mom_vert_arr[ind];
}

void mirror_mom_diag(int& ind)
{
   ind = mirror_mom_diag_arr[ind];
}

void mirror_mom_pipi(int& ind)
{   
   ind = mirror_mom_pipi_arr[ind]; 
}

void swap(int& a, int& b)
{
   int temp = a;
   a = b;
   b = temp;
}


