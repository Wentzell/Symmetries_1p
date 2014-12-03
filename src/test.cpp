
/************************************************************************************************//**
 *  		
 * 	file: 		test.cpp
 * 	contents:   	Test source to use symmetries of vertex object
 * 
 ****************************************************************************************************/

#include <symmetries.h>
#include <iostream>

using namespace std;

/**
 *	Test
 */
int main (int argc, char * argv[])
{
   // Initialize vertex tensor
   shared_ptr<se_tensor> vertex_ptr(new se_tensor(FREQ_COUNT_SE, PATCH_COUNT, QN_COUNT));

   // Initialize independent coupling vector
   vector<index_t> ind_cpl_list;

   init_symm( vertex_ptr, ind_cpl_list );

   std::cout << " Number of tensor indices : " << TENSOR_IND_COUNT << endl;
   std::cout << " Number of independent couplings : " << ind_cpl_list.size() << endl;
   std::cout << " Reduction factor: "	<< 1.0*TENSOR_IND_COUNT/ ind_cpl_list.size() << endl;

   return 0;
}

