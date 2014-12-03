
/************************************************************************************************//**
 *  		
 * 	file: 		symmap.cpp
 * 	contents:   	See symmap.h
 * 
 ****************************************************************************************************/

#include <symmap.h>


std::ostream &operator<<(std::ostream& os, const index_1p_t& ind)
{
   os << " w " << ind.w << " k " << ind.k1_in << " s_in " << ind.s_in << " s_out " << ind.s_out<< std::endl;
   return os;
}
