
/************************************************************************************************//**
 *  		
 * 	file: 		const.h
 * 	contents:  	Contains numerical constants
 * 
 ****************************************************************************************************/


#ifndef CONST_H
#define CONST_H

const int POS_FREQ_COUNT_VERT = 20; 	///< Amount of positive matsubara frequencies for vertrex
const int FREQ_COUNT_VERT = 2*POS_FREQ_COUNT_VERT;	///< Amount of frequencies including the negative ones
const int PATCH_COUNT = 8;		///< Amount of k-patches
const int QN_COUNT = 1;			///< Amount of possible tuples of the discrete quantum numbers

///< Amount of tensor components
const int TENSOR_IND_COUNT = FREQ_COUNT_VERT * FREQ_COUNT_VERT * FREQ_COUNT_VERT * PATCH_COUNT * PATCH_COUNT * PATCH_COUNT * QN_COUNT  * QN_COUNT * QN_COUNT * QN_COUNT;

#endif 
