/*
 *  Memory operations tests.
 *
 *  Created by Lance Putnam on 2/27/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdio.h>
#include <string.h>
#include "arr.h"
#include "mem.h"

#define SET_M1 mem::set(arrE, val(-1), Loop(lenE)); mem::set(arrO, val(-1), Loop(lenO));
#define SET_LINE arr::lineSlope1(arrE, lenE); arr::lineSlope1(arrO, lenO);
#define PRINT_EVEN for(i=0; i<lenE; i++) printf("% 2d ", arrE[i]); printf("\n");
#define PRINT_ODD  for(i=0; i<lenO; i++) printf("% 2d ", arrO[i]); printf("\n");
#define PRINT_IND  for(i=0; i<lenO; i++) printf("% 2d ", i); printf("\n");
#define PRINT_BOTH PRINT_IND PRINT_EVEN PRINT_ODD
	
#define ASSERT_GUARDS \
	if(arrE[lenE] != 100){ printf("**arrE guard modified to % 2d\n", arrE[lenE]); return 1;}\
	if(arrO[lenO] != 100){ printf("**arrO guard modified to % 2d\n", arrO[lenO]); return 1;}

using namespace gam;
using namespace gam::gen;

int main(int argc, char* argv[]){
	int i=0;

	const int lenE = 8;
	const int lenO = lenE + 1;
	const int lengthE2 = lenE>>1;
	//const int lengthO2 = lenO>>1;
	
	int arrE[lenE + 1];
	int arrO[lenO + 1];
	int tempE[lenE];
	int tempO[lenO];

	SET_M1

	// Set guard elements
	arrE[lenE] = 100;
	arrO[lenO] = 100;

	printf("\nCopy:\n");
	arr::lineSlope1(tempE, lenE);
	arr::lineSlope1(tempO, lenO);
	mem::copy(arrE, tempE, lenE);
	mem::copy(arrO, tempO, lenO);
	PRINT_BOTH ASSERT_GUARDS

	printf("\nInterleave 2:\n");
	SET_LINE
	mem::interleave2(arrE, tempE, lengthE2);
	PRINT_EVEN ASSERT_GUARDS

	printf("\nDeinterleave 2:\n");
	mem::copy(tempE, arrE, lenE);
	mem::deinterleave2(arrE, tempE, lengthE2);
	PRINT_EVEN ASSERT_GUARDS

//	printf("\nEquals:\n");
//	arr::lineSlope1(tempE, lenE);
//	arr::lineSlope1(tempO, lenO);
//	SET_LINE
//	printf("%s %s\n", 
//		mem::equal(arrE, tempE, lenE) ? "true" : "false",
//		mem::equal(arrO, tempO, lenO) ? "true" : "false"
//	);

	printf("\nExpand (2):\n");
	arr::lineSlope1(tempE, lenE);
	SET_M1
	mem::expand(arrE, tempE, lenE/2, 2UL);
	PRINT_EVEN ASSERT_GUARDS

	printf("\nKeep (stride=2, offset=0):\n");
	SET_M1
	mem::keep(arrE, lenE, 2, 0);
	mem::keep(arrO, lenO, 2, 0);
	PRINT_BOTH ASSERT_GUARDS

	printf("\nKeep (stride=3, offset=2):\n");
	SET_M1
	mem::keep(arrE, lenE, 3, 2);
	mem::keep(arrO, lenO, 3, 2);
	PRINT_BOTH ASSERT_GUARDS

	printf("\nMirror right:\n");
	SET_LINE
	mem::mirrorR(arrE, lenE);
	mem::mirrorR(arrO, lenO);
	PRINT_BOTH ASSERT_GUARDS

	printf("\nPivot around 2:\n");
	SET_LINE
	mem::pivot(arrE, lenE, 2);
	mem::pivot(arrO, lenO, 2);
	PRINT_BOTH ASSERT_GUARDS

	printf("\nReplace:\n");
	SET_M1
	mem::replace(arrE, lenE, -1, 0);
	mem::replace(arrO, lenO, -1, 0);
	PRINT_BOTH ASSERT_GUARDS
	
//	char test[] = "abcdabcd";
//	mem::replace(test, strlen(test), 'a', 'z');
//	printf("%s\n", test);

	printf("\nReverse:\n");
	SET_LINE
	mem::reverse(arrE, lenE);
	mem::reverse(arrO, lenO);
	PRINT_BOTH ASSERT_GUARDS

	printf("\nReverse every 2:\n");
	SET_LINE
	mem::reverse2(arrE, lenE);
	mem::reverse2(arrO, lenO);
	PRINT_BOTH ASSERT_GUARDS

	printf("\nRotate half:\n");
	SET_LINE
	mem::rotateH(arrE, lenE);
	mem::rotateH(arrO, lenO);
	PRINT_BOTH ASSERT_GUARDS

	printf("\nRotate left by 1:\n");
	SET_LINE
	mem::rotateL1(arrE, lenE);
	mem::rotateL1(arrO, lenO);
	PRINT_BOTH ASSERT_GUARDS
	
	printf("\nRotate left by 3:\n");
	SET_LINE
	mem::rotateL(arrE, lenE, 0);
	mem::rotateL(arrO, lenO, 0);
	PRINT_BOTH ASSERT_GUARDS
	
	printf("\nRotate right by 1:\n");
	SET_LINE
	mem::rotateR1(arrE, lenE);
	mem::rotateR1(arrO, lenO);
	PRINT_BOTH ASSERT_GUARDS

	printf("\nScale (cropping) by 2.5:\n");
	arr::lineSlope1(tempE, lenE);
	SET_M1
	mem::scaleCrop(arrE, tempE, lenE, 2.5f);
	PRINT_EVEN ASSERT_GUARDS

	printf("\nSet to 1:\n");
	mem::set(arrE, val(1), Loop(lenE));
	mem::set(arrO, val(1), Loop(lenO));
	PRINT_BOTH ASSERT_GUARDS

	printf("\nSet to 1 (stride=2, offset=0):\n");
	SET_M1
	mem::set(arrE, val(1), Loop(lenE,2));
	mem::set(arrO, val(1), Loop(lenO,2));
	PRINT_BOTH ASSERT_GUARDS

	printf("\nSet to 1 (stride=2, offset=1):\n");
	SET_M1
	mem::set(arrE, val(1), Loop(lenE,2,1));
	mem::set(arrO, val(1), Loop(lenO,2,1));
	PRINT_BOTH ASSERT_GUARDS

	printf("\nStretch (in-place):\n");
	SET_LINE
	mem::stretch(arrE, lenE, 2);
	PRINT_EVEN ASSERT_GUARDS
	
	printf("\nStretch (out-place):\n");
	arr::lineSlope1(tempE, lenE);
	mem::zero(arrE, lenE);
	mem::stretch(arrE, tempE, lenE/2, 2);
	PRINT_EVEN ASSERT_GUARDS

	printf("\nSwap (odd -line with even +line):\n");
	arr::lineSlope1(arrE, lenE, 1);
	arr::lineSlope(arrO, lenO, -1, -1);
	mem::swap(arrE, arrO, lenE);
	PRINT_BOTH ASSERT_GUARDS
	
	printf("\nTranspose2:\n");
	SET_LINE
	mem::transpose2(arrE, lenE);
	PRINT_EVEN ASSERT_GUARDS

	printf("\nZero:\n");
	SET_M1
	mem::zero(arrE, lenE);
	PRINT_EVEN ASSERT_GUARDS

	printf("\nZero (stride=2, offset=0):\n");
	SET_M1
	mem::set(arrE, val(0), Loop(lenE,2));
	PRINT_EVEN ASSERT_GUARDS

	printf("\nZero (stride=2, offset=1):\n");
	SET_M1
	mem::set(arrE, val(0), Loop(lenE,2,1));
	PRINT_EVEN ASSERT_GUARDS


	{
		using namespace gam::mem;
		using namespace gam::gen;
	
		#define PRINT(obj) printf("\t"); for(uint32_t i=0;i<obj.size();++i) printf("%c", obj[i]); printf("\n");
	
		Multi<8, char> multi;
		
		printf("\nmem::set\n");
		set(multi, "12345678", Loop(multi.size())); PRINT(multi)		
		set(multi, val('1'), Loop(multi.size())); PRINT(multi)
		set(multi, rAdd1('a'), Loop(multi.size())); PRINT(multi)
		set(multi, rAdd((char)2,'a'), Loop(multi.size())); PRINT(multi)
		set(multi, rAdd((char)2,'b'), Loop(multi.size())); PRINT(multi)
		
		set(multi, val('.'), Loop(multi.size()));
		set(multi, val('1'), Loop(multi.size(),2)); PRINT(multi)
		set(multi, val('.'), Loop(multi.size()));
		set(multi, val('1'), Loop(multi.size(),3)); PRINT(multi)
	}	

	return 0;
}

