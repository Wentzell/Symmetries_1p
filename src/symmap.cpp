
/************************************************************************************************//**
 *  		
 * 	file: 		symmetries.cpp
 * 	contents:   	See symmetries.h
 * 
 ****************************************************************************************************/

#include <symmetries.h>


std::ostream &operator<<(std::ostream& os, const index_t& ind)
{
   os << " w1_in " << ind.w1_in << " w2_in " << ind.w2_in << " w1_out " << ind.w1_out << " k1_in " << ind.k1_in << " k2_in " << ind.k2_in << " k1_ou " << ind.k1_out << " s1_in " << ind.s1_in << " s2_in " << ind.s2_in << " s1_ou " << ind.s1_out << " s2_ou " << ind.s2_out << std::endl;
   return os;
}
