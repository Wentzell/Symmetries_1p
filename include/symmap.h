
/************************************************************************************************//**
 *  		
 * 	file: 		symmap.h
 * 	contents:  	Use index symmetries of a two-particle vertex tensor to determine a minimal 
 * 			set of independent elements	
 * 
 ****************************************************************************************************/


#ifndef SYMMAP_H
#define SYMMAP_H

#include <boost/multi_array.hpp>

/**
 *	Class defines a type that contains a set of vertex indeces specifying a tensor element. Symmetries can act on it.	
 */
class index_t  // : public inherits
{
   public:
      int w1_in, w2_in, w1_out; ///< In and outgoing frequency indeces. Correspond to Matsubara frequencys caluclated according 2\Pi/\beta(n + 1/2)
      int k1_in, k2_in, k1_out;	///< In and outgoing momentum patch indeces.
      int s1_in, s2_in, s1_out, s2_out;	///< In and outgoing discrete quantum numbers. Correspond to tupels of e.g. spin, orbital ...

      ///< Constructor for index_t
      index_t(int w1_in_, int w2_in_, int w1_out_, int k1_in_, int k2_in_, int k1_out_, int s1_in_, int s2_in_, int s1_out_, int s2_out_) :
	 w1_in(w1_in_), w2_in(w2_in_), w1_out(w1_out_), k1_in(k1_in_), k2_in(k2_in_), k1_out(k1_out_), s1_in(s1_in_), s2_in(s2_in_), s1_out(s1_out_), s2_out(s2_out_)
   {}

      /********************* specify possible index operations ********************/ 

   protected:

   private:

      friend std::ostream &operator<<(std::ostream&, const index_t&); ///< Define output operator << for index_t
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

   protected:

   private:

};

//operation  operator*(const operation& a, const operation& b)
//{
//return operation(  a.first xor b.first, a.second xor b.second ) ;
//}



/**
 *	Elements of vertex tensor. States index of independent coupling array and possible operations on it
 */
struct ind_cpl_t
{
   // CHANGE !!!! IMPLEMENT HERE UNSIGNED SHORT INT
   int ind;	///< Index in the vector of independent couplings.
   operation oper; ///< Possible operations that relate two tensor elements. First bool indicates possible sign change, second one complex conjugation

   ///< Constructor int, bool, bool
   ind_cpl_t(int ind_ = -1, bool first_ = false , bool second_ = false ):
      ind(ind_), oper(first_, second_)
   {}
   ///< Constructor int, operation
   ind_cpl_t(int ind_, operation oper_):
      ind(ind_), oper(oper_)
   {}
};


/**
 *	Class representing the vertex tensor
 */
class vertex_tensor : public boost::multi_array<ind_cpl_t, 10>
{
   public:

      typedef boost::multi_array<ind_cpl_t, 10> super;

      ///< Allow access to elements in tensor by means of index object
      ind_cpl_t& operator()(index_t& ind)
      {
	 return (*this)[ind.w1_in][ind.w2_in][ind.w1_out][ind.k1_in][ind.k2_in][ind.k1_out][ind.s1_in][ind.s2_in][ind.s1_out][ind.s2_out];
      }

      ///< Allow access to elements in tensor by specifying all indeces
      ind_cpl_t& operator()(int w1_in, int w2_in, int w1_out, int k1_in, int k2_in, int k1_out, int s1_in, int s2_in, int s1_out, int s2_out) 
      {
	 return (*this)[w1_in][w2_in][w1_out][k1_in][k2_in][k1_out][s1_in][s2_in][s1_out][s2_out];
      }

      ///< Define more convenient constructor
      vertex_tensor(int dim_w, int dim_k, int dim_s):
	 super(boost::extents[dim_w][dim_w][dim_w][dim_k][dim_k][dim_k][dim_s][dim_s][dim_s][dim_s])
   {}

   protected:

   private:

};

#endif 
