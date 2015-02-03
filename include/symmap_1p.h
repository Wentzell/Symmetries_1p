
/*******************************************************************************************//** @file
 *  		
 * 	file: 		symmap_1p.h
 * 	contents: 	Definition of index_1p_t class, represents index object for one-particle properties. 
 * 			Definition of se_tensor	class, that allows creation of look-up arrays for the self-energy
 * 
 ****************************************************************************************************/


#ifndef SYMMAP_1P_H
#define SYMMAP_1P_H

#include <boost/multi_array.hpp>
#include <symbase.h>
#include <iostream>
/**
 *	Class defines a type that contains a set of self-energy indeces specifying a tensor element. Symmetries can act on it.	
 */
class index_1p_t  // : public inherits
{
   public:
      int w; ///< Frequency index. Correspond to Matsubara frequencys caluclated according 2\Pi/\beta(n + 1/2)

      ///< Constructor for index_1p_t
      index_1p_t(int w_) :
	 w(w_)
   {}

   private:
      friend std::ostream &operator<<(std::ostream&, const index_1p_t&); ///< Define output operator << for index_1p_t
};

/**
 *	Class representing the vertex tensor
 */
class se_tensor : public boost::multi_array<ind_cpl_t, 1>
{
   public:
      typedef boost::multi_array<ind_cpl_t, 1> super;

      ///< Allow access to elements in tensor by means of index object
      ind_cpl_t& operator()(index_1p_t& ind)
      {
	 return (*this)[ind.w];
      }

      ///< Allow access to elements in tensor by specifying all indeces
      ind_cpl_t& operator()(int w) 
      {
	 return (*this)[w];
      }

      ///< Define more convenient constructor
      se_tensor(int dim_w):
	 super(boost::extents[dim_w])
   {}

   private:
};

/**
 *	Class representing the vertex tensor
 */
class se_checked_tensor : public boost::multi_array<bool, 1>
{
   public:
      typedef boost::multi_array<bool, 1> super;

      ///< Allow access to elements in tensor by means of index object
      bool& operator()(index_1p_t& ind)
      {
	 return (*this)[ind.w];
      }

      ///< Allow access to elements in tensor by specifying all indeces
      bool& operator()(int w) 
      {
	 return (*this)[w];
      }

      ///< Define more convenient constructor
      se_checked_tensor(int dim_w):
	 super(boost::extents[dim_w])
   { std::fill( this->data(), this->data() + this->num_elements(), false); }

   private:
};
#endif 
