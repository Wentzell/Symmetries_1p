
/************************************************************************************************//**
 *  		
 * 	file: 		symmap.h
 * 	contents: 	Definition of vertex-energy index class and clases necessary to establish
 * 			the symmetry mapping	
 * 
 ****************************************************************************************************/


#ifndef SYMMAP_H
#define SYMMAP_H

#include <boost/multi_array.hpp>

/**
 *	Class defines a type that contains a set of self-energy indeces specifying a tensor element. Symmetries can act on it.	
 */
class index_1p_t  // : public inherits
{
   public:
      int w; ///< Frequency index. Correspond to Matsubara frequencys caluclated according 2\Pi/\beta(n + 1/2)
      int k;	///< Momentum patch index
      int s_in, s_out;	///< In and outgoing discrete quantum numbers. Correspond to tupels of e.g. spin, orbital ...

      ///< Constructor for index_2p_t
      index_2p_t(int w_, int k_, int s_in_, int s_out_) :
	 w(w_), k(k_), s_in(s_in_), s_out(s_out_)
   {}

   private:
      friend std::ostream &operator<<(std::ostream&, const index_1p_t&); ///< Define output operator << for index_2p_t
};



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

/**
 *	Class representing the vertex tensor
 */
class se_tensor : public boost::multi_array<ind_cpl_t, 4>
{
   public:
      typedef boost::multi_array<ind_cpl_t, 4> super;

      ///< Allow access to elements in tensor by means of index object
      ind_cpl_t& operator()(index_2p_t& ind)
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
