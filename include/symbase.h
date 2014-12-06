
/*******************************************************************************************//** @file
 *  		
 * 	file: 		symbase.h
 * 	contents:  	Types and helper functions common to the 1P as well as the 2P Symmetry code
 * 
 ****************************************************************************************************/


#ifndef SYMBASE_H
#define SYMBASE_H

#include <utility>

/**
 *	Set of possible operations after symmetry operation
 */
class operation :  public std::pair<bool,bool>
{
   public:
      ///< Constructor taking two bools as argument
      operation(const bool& first_, const bool& second_)
      {
	 (*this).first = first_ ;
	 (*this).second = second_;
      }

      ///< Overload multiplication operator for successive application of operations
      operation  operator*(const operation& b)
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
      bool checked; 	///< States wether element has been checked for symmetries
      operation oper; 	///< Possible operations that relate two tensor elements. First bool indicates possible sign change, second one complex conjugation

      ///< Default Constructor 
      ind_cpl_t():
	 ind(0), checked(false), oper(false, false)
   {}

      ///< Constructor int, bool, bool
      ind_cpl_t(int ind_ , bool first_ = false , bool second_ = false ):
	 ind(ind_), checked(true), oper(first_, second_)
   {}

      ///< Constructor int, operation
      ind_cpl_t(int ind_, operation oper_):
	 ind(ind_), checked(true), oper(oper_)
   {}
};

// Helper functions and necessary arrays

void freq_sign_change(int& ind);	///< Change sign of signle frequency
void mom_sign_change(int& ind);		///< Change sign of single momentum
void mirror_mom_vert(int& ind);		///< Mirror momentum at vertical axis
void mirror_mom_diag(int& ind);		///< Mirror momentum at diagonal (bottom left to top right) axis
void mirror_mom_pipi(int& ind);		///< Swap momentum index for the one obtained by taking (pi,pi) - k
void swap(int& a, int& b);		///< Swap two numbers

const int rot_k_ind_arr[8] = {0, 2, 3, 4, 1, 7, 6, 5}; ///< Array that specifies how to rotate single momentum index
const int sign_change_k_ind_arr[8] = {0, 3, 4, 1, 2, 5, 6, 7}; ///< Array that specifies how to change sign of single momentum index ( corresponds to 2 rotations )
const int mirror_mom_vert_arr[8] = { 0, 4, 3, 2, 1, 5, 6, 7 }; ///< Array that specifies how to mirror single momentum index at vertical axis
const int mirror_mom_diag_arr[8] = { 0, 1, 4, 3, 2, 7, 6, 5 }; ///< Array that specifies how to mirror single momentum index at diagonal ( bottom left to top right ) axis
const int mirror_mom_pipi_arr[8] = { 6, 1, 2, 3, 4, 7, 0, 5 }; ///< Array that specifies how to calculate (pi,pi) - k for single momentum k

#endif 
