
/*******************************************************************************************//** @file
 *  		
 * 	file: 		symbase.h
 * 	contents:  	Types and helper functions and arrays common to the 1P as well as the 2P Symmetry code.
 * 			Contains class operation that possible operations attached to a symmetry operation
 * 			(sign change, complex conjugation). Also contains the struct ind_cpl_t that corresponds
 * 			to an index in the independent coupling vector with attached operations.
 * 			Other functions and arrays specified here define how to perform certain operations on
 * 			the frequencies, momenta and the discrete quantum numbers.
 * 
 ****************************************************************************************************/


#ifndef SYMBASE_H
#define SYMBASE_H

#include <utility>

/**
 *	Set of possible operations after symmetry operation
 */
class operation : public std::pair<bool,bool>
{
   public:
      operation(const bool& first_, const bool& second_) ///< Constructor taking two bools as argument
      {
	 (*this).first = first_ ;
	 (*this).second = second_;
      }

      operation  operator*(const operation& b)	///< Overload multiplication operator for successive application of operations
      {
	 return operation(  (*this).first xor b.first, (*this).second xor b.second ) ;
      }

   private:
};


/**
 *	Elements of vertex tensor. States index of independent coupling array and possible operations on it
 */
struct ind_cpl_t
{
   public:
      unsigned int ind;	///< Index in the vector of independent couplings.
      bool forced_zero;	///< States weather this element is explicitly forced to zero
      operation oper; 	///< Possible operations that relate two tensor elements. First bool indicates possible sign change, second one complex conjugation

      /**
       *	Default Constructor, elements default as forced to zero	
       */
      ind_cpl_t():
	 ind(0), forced_zero(true), oper(false, false)
   {}

      /**
       *	Constructor int, bool, bool	
       */
      ind_cpl_t(int ind_ , bool first_ = false , bool second_ = false ):
	 ind(ind_), forced_zero(false), oper(first_, second_)
   {}

      /**
       *		Constructor (int, operation)
       */   
      ind_cpl_t(int ind_, operation oper_):
	 ind(ind_), forced_zero(false), oper(oper_)
   {}
};

// Helper functions and necessary arrays

void freq_sign_change(int& ind, const int freq_count);	///< Change sign of signle frequency, MOVE TO FREQUENCY GRID
void flip_spin(int& ind);		///< Flip single spin index
void mom_sign_change(int& ind);		///< Change sign of single momentum, MOVE TO MOMENTUM GRID
void mirror_mom_vert(int& ind);		///< Mirror momentum at vertical axis, MOVE TO MOMENTUM GRID
void mirror_mom_diag(int& ind);		///< Mirror momentum at diagonal (bottom left to top right) axis, MOVE TO MOMENTUM GRID
void mirror_mom_pipi(int& ind);		///< Swap momentum index for the one obtained by taking (pi,pi) - k MOVE TO MOMENTUM GRID
void swap(int& a, int& b);		///< Swap two numbers 

const int rot_k_ind_arr[8] = {0, 2, 3, 4, 1, 7, 6, 5}; ///< Array that specifies how to rotate single momentum index
const int sign_change_k_ind_arr[8] = {0, 3, 4, 1, 2, 5, 6, 7}; ///< Array that specifies how to change sign of single momentum index ( corresponds to 2 rotations )
const int mirror_mom_vert_arr[8] = { 0, 4, 3, 2, 1, 5, 6, 7 }; ///< Array that specifies how to mirror single momentum index at vertical axis
const int mirror_mom_diag_arr[8] = { 0, 1, 4, 3, 2, 7, 6, 5 }; ///< Array that specifies how to mirror single momentum index at diagonal ( bottom left to top right ) axis
const int mirror_mom_pipi_arr[8] = { 6, 1, 2, 3, 4, 7, 0, 5 }; ///< Array that specifies how to calculate (pi,pi) - k for single momentum k

/**
 *	Array that specifies how to add to momenta in the language of the patch indeces
 */
const int sum_mom[8][8] = { 	{ 0, 1, 2, 3, 4, 5, 6, 7 }, { 1, 6, 7, 0, 5, 2, 3, 4 }, { 2, 7, 6, 5, 0, 1, 4, 3 }, { 3, 0, 5, 6, 7, 4, 1, 2 }, 
   { 4, 5, 0, 7, 6, 3, 2, 1 }, { 5, 2, 1, 4, 3, 0, 7, 6 }, { 6, 3, 4, 1, 2, 7, 0, 5 }, { 7, 4, 3, 2, 1, 6, 5, 0 } } ;

#endif 
