
/************************************************************************************************//**
 *  		
 * 	file: 		symmetries.cpp
 * 	contents:   	see symmetries.h
 * 
 ****************************************************************************************************/

#include <symmetries.h>


void iterate( const index_t& ind, const operation& track_op, vertex_tensor& vertex, std::vector<symm_func_t> symm_func_list , int ind_cpl_list_pos )
{
   for(auto symm_func: symm_func_list) 		// iterate over list of all symmetries specified
   {
      index_t ind_it = ind;			// copy ind
      operation curr_op = symm_func(ind_it) * track_op;	// apply symmetry operation and track operations applied
      if( vertex(ind_it).ind == -1)		// if resulting tensor index not yet related to any other
      { 
	 vertex(ind_it) = ind_cpl_t(ind_cpl_list_pos, curr_op); 		// relate to position
	 iterate( ind_it, curr_op, vertex, symm_func_list, ind_cpl_list_pos );	// iterate further	
      }
   }
}

//----- Antisymmetry

operation exch_in(index_t& ind)
{
   swap(ind.w1_in, ind.w2_in);
   swap(ind.k1_in, ind.k2_in);
   swap(ind.s1_in, ind.s2_in);
   return operation(true,false);
}

operation exch_out(index_t& ind) 
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

operation compl_conj(index_t& ind)
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

      return operation(false,true);
   }
   return operation(false, false);
}

operation time_rev(index_t& ind)
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

operation particle_hole(index_t& ind)
{
   mirror_mom_pipi(ind.k1_in);
   mirror_mom_pipi(ind.k2_in);
   mirror_mom_pipi(ind.k1_out);
   return operation(false,true);
}

operation rot_k(index_t& ind)
{
   ind.k1_in = rot_k_ind_arr[ind.k1_in];
   ind.k2_in = rot_k_ind_arr[ind.k2_in];
   ind.k1_out = rot_k_ind_arr[ind.k1_out];
   return operation(false,false);
}

operation mirror_vert(index_t& ind)
{
   mirror_mom_vert(ind.k1_in);
   mirror_mom_vert(ind.k2_in);
   mirror_mom_vert(ind.k1_out);
   return operation(false,false);
}

operation mirror_diag(index_t& ind)
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


