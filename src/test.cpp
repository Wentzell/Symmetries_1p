
/************************************************************************************************//**
 *  		
 * 	file: 		test.cpp
 * 	contents:   	Test source to use symmetries of vertex object
 * 
 ****************************************************************************************************/

#include <symmetries_1p.h>
#include <iostream>

using namespace std;

/**
 *	Test
 */
int main (int argc, char * argv[])
{
   // Initialize vertex tensor
   shared_ptr<se_tensor> se_ptr(new se_tensor(FREQ_COUNT_SE, PATCH_COUNT, QN_COUNT));

   // Initialize independent coupling vector
   vector<index_1p_t> ind_cpl_list;

   init_symm( se_ptr, ind_cpl_list );

   std::cout << " Number of tensor indices : " <<  FREQ_COUNT_SE * PATCH_COUNT * QN_COUNT * QN_COUNT << endl;
   std::cout << " Number of independent couplings : " << ind_cpl_list.size() << endl;
   std::cout << " Reduction factor: "	<< FREQ_COUNT_SE * PATCH_COUNT * QN_COUNT * QN_COUNT / ind_cpl_list.size() << endl;

   return 0;
}

