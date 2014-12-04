
/************************************************************************************************//**
 *  		
 * 	file: 		symmap_1p.cpp
 * 	contents:   	See symmap_1p.h
 * 
 ****************************************************************************************************/

#include <symmap_1p.h>


std::ostream &operator<<(std::ostream& os, const index_1p_t& ind)
{
   os << " w " << ind.w << " k " << ind.k << " s_in " << ind.s_in << " s_out " << ind.s_out<< std::endl;
   return os;
}
