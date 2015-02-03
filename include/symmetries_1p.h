
/*******************************************************************************************//** @file
 *  		
 * 	file: 		symmetries_1p.h
 * 	contents:  	For a given set of symmetries, provides the functions to establish the linking 
 * 			between the self-energy tensor and the independent couplings.
 * 			Also defines all the available symmetry operations and activates / deactivates
 * 			them
 * 
 ****************************************************************************************************/


#ifndef SYMMETRIES_1P_H
#define SYMMETRIES_1P_H

#include <symmap_1p.h>
#include <const.h>

typedef operation (*symm_func_1p_t)(index_1p_t&); ///< Symmetry function that acts on an index_1p_t object and alters it 

/**
 *	Apply symmetries and establish mapping	
 */
void init_symm( std::shared_ptr<se_tensor> se_ptr, std::vector<index_1p_t>& ind_cpl_list ); 
/**
 *	Iterate a symmetry on vertex object	
 */
void iterate( const index_1p_t& ind, const operation& track_op, se_tensor& vertex, se_checked_tensor& checked, std::vector<symm_func_1p_t> symm_func_list , int ind_cpl_list_pos );

// Symmetries
operation compl_conj(index_1p_t& ind);		///< Complex conjugation 
operation particle_hole(index_1p_t& ind);	///< Particle hole symmetry

#endif 
