
/************************************************************************************************//**
 *  		
 * 	file: 		symmap_1p.h
 * 	contents: 	Definition of vertex-energy index class and clases necessary to establish
 * 			the symmetry mapping	
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
      int k;	///< Momentum patch index
      int s_in, s_out;	///< In and outgoing discrete quantum numbers. Correspond to tupels of e.g. spin, orbital ...

      ///< Constructor for index_1p_t
      index_1p_t(int w_, int k_, int s_in_, int s_out_) :
	 w(w_), k(k_), s_in(s_in_), s_out(s_out_)
   {}

   private:
      friend std::ostream &operator<<(std::ostream&, const index_1p_t&); ///< Define output operator << for index_1p_t
};

/**
 *	Class representing the vertex tensor
 */
class se_tensor : public boost::multi_array<ind_cpl_t, 4>
{
   public:
      typedef boost::multi_array<ind_cpl_t, 4> super;

      ///< Allow access to elements in tensor by means of index object
      ind_cpl_t& operator()(index_1p_t& ind)
      {
	 return (*this)[ind.w][ind.k][ind.s_in][ind.s_out];
      }

      ///< Allow access to elements in tensor by specifying all indeces
      ind_cpl_t& operator()(int w, int k, int s_in, int s_out) 
      {
	 return (*this)[w][k][s_in][s_out];
      }

      ///< Define more convenient constructor
      se_tensor(int dim_w, int dim_k, int dim_s):
	 super(boost::extents[dim_w][dim_k][dim_s][dim_s])
   {}

   private:
};

#endif 
