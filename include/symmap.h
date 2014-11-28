
/************************************************************************************************//**
 *  		
 * 	file: 		symmap.h
 * 	contents:  	For a given set of symmetries, establish the linking between the 
 * 			vertex tensor and the independent couplings	
 * 
 ****************************************************************************************************/


#ifndef SYMMAP_H
#define SYMMAP_H

#include <symmetries.h>
#include <const.h>

typedef operation (*symm_func_t)(index_t&); 		///< Symmetry function that acts on an index_t object and alters it 

///< Iterate a symmetry on vertex object
void iterate( const index_t& ind, const operation& track_op, vertex_tensor& vertex, std::vector<symm_func_t> symm_func_list , int ind_cpl_list_pos );

// Symmetries
operation exch_in(index_t& ind);	///< Exchange ingoing lines
operation rot_k(index_t& ind);		///< Rotate all momenta by 90 degrees - IMPLEMENT USING STATIC VECTOR
operation compl_conj(index_t& ind);	///< Complex conjugation
operation exch_out(index_t& ind);	///< Exchange outgoing lines
operation time_rev(index_t& ind);	///< Time reversal symmetry


// Helper functions
void freq_sign_change(int& ind);	///< Change sign of signle frequency
void mom_sign_change(int& ind);		///< Change sign of single momentum
void swap(int& a, int& b);		///< Swap two numbers

const int rot_k_ind_arr[8] = {0, 2,3,4,1, 7, 6, 5}; ///< Array that specifies how to rotate single momentum index
const int sign_change_k_ind_arr[8] = {0, 3, 4, 1, 2, 5, 6, 7}; ///< Array that specifies how to change sign of single momentum index ( corresponds to 2 rotations)
/**
 *	Calculate the sum of two momenta
 */
const int sum_mom[8][8] = { 	{ 0, 1, 2, 3, 4, 5, 6, 7 }, { 1, 6, 7, 0, 5, 2, 3, 4 }, { 2, 7, 4, 5, 0, 1, 4, 3 }, { 3, 0, 5, 6, 7, 4, 1, 2 }, 
   					{ 4, 5, 0, 7, 6, 3, 2, 1 }, { 5, 2, 1, 4, 3, 0, 7, 6 }, { 6, 3, 4, 1, 2, 7, 0, 5 }, { 7, 4, 3, 2, 1, 6, 5, 0 } } ;


#endif 
