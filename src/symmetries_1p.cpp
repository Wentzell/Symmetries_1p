
/************************************************************************************************//**
 *  		
 * 	file: 		symmetries_1p.cpp
 * 	contents:   	see symmetries_1p.h
 * 
 ****************************************************************************************************/


#include <symmetries_1p.h>
#include <def.h>

using namespace std;

void init_symm( shared_ptr<se_tensor> se_ptr, vector<idx_1p_t>& ind_cpl_list )
{
   // Initialize vector containing all symmetry functions 
   vector<symm_func_1p_t> symm_func_list; // = { exch_in, exch_out, compl_conj, time_rev, particle_hole, rot_k, mirror_vert, mirror_diag };

   symm_func_list.push_back( compl_conj );
   symm_func_list.push_back( time_rev );
   symm_func_list.push_back( particle_hole );
   symm_func_list.push_back( spin_symm );

#ifndef NO_MOMENTA
   symm_func_list.push_back( rot_k );
   symm_func_list.push_back( mirror_vert );
#endif

   se_checked_tensor checked( POS_FREQ_COUNT_SE, PATCH_COUNT, QN_COUNT ); 

   for( int w = -POS_FREQ_COUNT_SE; w < POS_FREQ_COUNT_SE; ++w )
      for( int k = 0; k < PATCH_COUNT; ++k )
	 for( int s_in = 0; s_in < QN_COUNT; ++s_in )
	    for( int s_out = 0; s_out < QN_COUNT; ++s_out )
	    {
	       idx_1p_t idx( w, k, s_in, s_out );
	       if ( !checked( idx ) )  						// if tensor object not yet related to any other
	       {
		  checked( idx ) = true; 
		  ind_cpl_list.push_back( idx );				// push standard representative into ind_cpl_list
		  int ind_cpl_list_pos =  ind_cpl_list.size() - 1 ; 		// position in ind_cpl_list is size - 1
		  (*se_ptr )( idx ) = ind_cpl_t( ind_cpl_list_pos ); 		// save position with trivial operations
		  operation track_op( false, false );
		  iterate( idx, track_op, *se_ptr, checked, symm_func_list, ind_cpl_list_pos ); 		// start iterating on index 
	       }
	    }


}

void iterate( const idx_1p_t& idx, const operation& track_op, se_tensor& selfEn, se_checked_tensor& checked, vector<symm_func_1p_t> symm_func_list, int ind_cpl_list_pos )
{
   for( auto symm_func: symm_func_list ) 				// iterate over list of all symmetries specified
   {
      idx_1p_t idx_it = idx;						// copy ind
      operation curr_op = symm_func( idx_it ) * track_op;		// apply symmetry operation and track operations applied
      if( !checked( idx_it ) )						// if resulting tensor index not yet related to any other
      { 
	 checked( idx_it ) = true; 
	 selfEn( idx_it ) = ind_cpl_t( ind_cpl_list_pos, curr_op ); 	// relate to position
	 iterate( idx_it, curr_op, selfEn, checked, symm_func_list, ind_cpl_list_pos );	// iterate further	
      }
   }
}

//----- Antisymmetry

operation compl_conj( idx_1p_t& idx )
{
   freq_sign_change( idx( w ) );

   // Changing momenta sign is unnecessary since equal to rotating twice by 90 degrees

   swap( idx( s_in ), idx( s_out ) );
   return operation( false, true );
}

operation time_rev( idx_1p_t& idx )
{
   swap( idx( s_in ), idx( s_out ) );
   return operation( false, false );
}

operation particle_hole( idx_1p_t& idx )
{
#ifndef NO_MOMENTA
   mirror_mom_pipi( idx( k ) );
#endif
   return operation( false, true );
}

operation spin_symm( idx_1p_t& idx )
{
   freq_sign_change( idx( w ) );
   
   swap( idx( s_in ), idx( s_out ) );

   flip_spin( idx( s_in ) );
   flip_spin( idx( s_out ) );

   if ( ( idx( s_in ) + !( idx( s_out ) ) ) % 2  != 0 ) // if uneven amount of creation operators minus sign, see masterthesis
      return operation( true, false );

   return operation( false, false );
}

operation rot_k( idx_1p_t& idx )
{
   idx( k ) = rot_k_idx_arr[idx( k )];
   return operation( false, false );
}

operation mirror_vert( idx_1p_t& idx )
{
   mirror_mom_vert( idx( k ) );
   return operation( false, false );
}

