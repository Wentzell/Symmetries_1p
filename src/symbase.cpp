
/************************************************************************************************//**
 *  		
 * 	file: 		symbase.cpp
 * 	contents:   	See symbase.h
 * 
 ****************************************************************************************************/

#include <symbase.h>
#include <const.h>

void freq_sign_change(int& ind, const int freq_count)
{
   ind = freq_count - ind - 1;
}

void mom_sign_change(int& ind)
{
   ind = sign_change_k_ind_arr[ind];
}

void mirror_mom_vert(int& ind)
{
   ind = mirror_mom_vert_arr[ind];
}

void mirror_mom_diag(int& ind)
{
   ind = mirror_mom_diag_arr[ind];
}

void mirror_mom_pipi(int& ind)
{   
   ind = mirror_mom_pipi_arr[ind]; 
}

void swap(int& a, int& b)
{
   int temp = a;
   a = b;
   b = temp;
}

