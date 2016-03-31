//#define SET_M1 mem::set(arrE, val(-1), Loop(lenE)); mem::set(arrO, val(-1), Loop(lenO));
//#define SET_LINE arr::lineSlope1(arrE, lenE); arr::lineSlope1(arrO, lenO);
//#define PRINT_EVEN for(i=0; i<lenE; i++) printf("% 2d ", arrE[i]); printf("\n");
//#define PRINT_ODD  for(i=0; i<lenO; i++) printf("% 2d ", arrO[i]); printf("\n");
//#define PRINT_IND  for(i=0; i<lenO; i++) printf("% 2d ", i); printf("\n");
//#define PRINT_BOTH PRINT_IND PRINT_EVEN PRINT_ODD
//	
//#define ASSERT_GUARDS 
//	if(arrE[lenE] != 100){ printf("**arrE guard modified to % 2d\n", arrE[lenE]); return 1;}
//	if(arrO[lenO] != 100){ printf("**arrO guard modified to % 2d\n", arrO[lenO]); return 1;}

{
	using namespace gam::gen;

	const unsigned N=8;
	int A[N], T[N];
	
	#define ASSERT(A, a,b,c,d,e,f,g,h) assert(A[0]==a && A[1]==b && A[2]==c && A[3]==d && A[4]==e && A[5]==f && A[6]==g && A[7]==h);

	slice(T,N) = RAdd1<int>();
	slice(A,N) = 0;
	mem::deepCopy(A,T,N); ASSERT(A, 0,1,2,3,4,5,6,7)

	slice(A,N) = 0;
	slice(A,N).copy(slice(T,N) = RAdd1<int>());
	assert(mem::deepEqual(A,T,N));

	slice(T,N) = RAdd1<int>();
	slice(A,N) = 0;
	mem::deepMove(A,T,N); ASSERT(A, 0,1,2,3,4,5,6,7)
	mem::deepMove(A,A+1,4); ASSERT(A, 1,2,3,4,4,5,6,7)

	slice(A,N) = RAdd1<int>(1);
	mem::deepZero(A,N); ASSERT(A, 0,0,0,0,0,0,0,0)

	slice(A,N) = RAdd1<int>();
	mem::pivot(3, A,N); ASSERT(A, 6,5,4,3,2,1,0,7)
	mem::pivot(0, A,N); ASSERT(A, 6,7,0,1,2,3,4,5)
	mem::pivot(7, A,N); ASSERT(A, 4,3,2,1,0,7,6,5)

	slice(A,N) = RAdd1<int>();
	mem::reflectLeft(A,N,1); ASSERT(A, 7,6,5,4,4,5,6,7)
	mem::reflectLeft(A,N,2); ASSERT(A, 6,6,4,4,4,5,6,7)

	slice(A,N) = RAdd1<int>();
	mem::reflectRight(A,N,1); ASSERT(A, 0,1,2,3,3,2,1,0)
	mem::reflectRight(A,N,2); ASSERT(A, 0,1,2,3,2,2,0,0)

	slice(A,N) = RAdd1<int>();
	mem::repeat(3, A,N); ASSERT(A, 0,1,2,0,1,2,0,1)
	mem::repeat(1, A,N); ASSERT(A, 0,0,0,0,0,0,0,0)

	slice(A, N/2) = 0;
	slice(A+N/2, N/2) = 1;
	mem::replace(1, 2, A, N,1); ASSERT(A, 0,0,0,0,2,2,2,2)
	{	int v[2]={0,2}, w[2]={1,3};
		mem::replace(v,w,2, A, N,1); ASSERT(A, 1,1,1,1,3,3,3,3)
	}

	slice(A,N) = RAdd1<int>();
	mem::reverse(A,N,1); ASSERT(A, 7,6,5,4,3,2,1,0)
	mem::reverse(A,N,2); ASSERT(A, 1,6,3,4,5,2,7,0)

	slice(A,N) = RAdd1<int>();
	mem::reverse2(A,N,1); ASSERT(A, 1,0,3,2,5,4,7,6)

	slice(A,N) = RAdd1<int>();
	mem::rotateHalf(A,N,1); ASSERT(A, 4,5,6,7,0,1,2,3)
	mem::rotateHalf(A,N,2); ASSERT(A, 0,5,2,7,4,1,6,3)

	slice(A,N) = RAdd1<int>();
	mem::rotateLeft(1, A,N); ASSERT(A, 1,2,3,4,5,6,7,0)
	mem::rotateLeft(2, A,N); ASSERT(A, 3,4,5,6,7,0,1,2)
	mem::rotateLeft(8, A,N); ASSERT(A, 3,4,5,6,7,0,1,2)

	slice(A,N) = RAdd1<int>();
	mem::rotateLeft1(A,N,1); ASSERT(A, 1,2,3,4,5,6,7,0)
	mem::rotateLeft1(A,N,2); ASSERT(A, 3,2,5,4,7,6,1,0)

	slice(A,N) = RAdd1<int>();
	mem::rotateRight1(A,N,1); ASSERT(A, 7,0,1,2,3,4,5,6)
	mem::rotateRight1(A,N,2); ASSERT(A, 5,0,7,2,1,4,3,6)

	slice(A,N) = RAdd1<int>();
	assert(mem::unique(A,N));

	slice(A,N) = RAdd1<int>();
	mem::zero(A,N,1); ASSERT(A, 0,0,0,0,0,0,0,0)
	slice(A,N) = RAdd1<int>();
	mem::zero(A,N,2); ASSERT(A, 0,1,0,3,0,5,0,7)

	#define T(x) assert(mem::select(x, 0, 1, 2) == x);
	T(0) T(1) T(2)
	#undef T

	#define T(x) assert(mem::select(x, 0, 1, 2, 3) == x);
	T(0) T(1) T(2) T(3)
	#undef T

	#undef ASSERT

//	int i=0;
//
//	const int lenE = 8;
//	const int lenO = lenE + 1;
//	const int lengthE2 = lenE>>1;
//	//const int lengthO2 = lenO>>1;
//	
//	int arrE[lenE + 1];
//	int arrO[lenO + 1];
//	int tempE[lenE];
//	int tempO[lenO];
//
//	SET_M1
//
//	// Set guard elements
//	arrE[lenE] = 100;
//	arrO[lenO] = 100;
//
//
//	printf("\nInterleave 2:\n");
//	SET_LINE
//	mem::interleave2(arrE, tempE, lengthE2);
//	PRINT_EVEN ASSERT_GUARDS
//
//	printf("\nDeinterleave 2:\n");
//	mem::copy(tempE, arrE, lenE);
//	mem::deinterleave2(arrE, tempE, lengthE2);
//	PRINT_EVEN ASSERT_GUARDS
//
//	printf("\nExpand (2):\n");
//	arr::lineSlope1(tempE, lenE);
//	SET_M1
//	mem::expand(arrE, tempE, lenE/2, 2UL);
//	PRINT_EVEN ASSERT_GUARDS
//
//	printf("\nKeep (stride=2, offset=0):\n");
//	SET_M1
//	mem::keep(arrE, lenE, 2, 0);
//	mem::keep(arrO, lenO, 2, 0);
//	PRINT_BOTH ASSERT_GUARDS
//
//	printf("\nKeep (stride=3, offset=2):\n");
//	SET_M1
//	mem::keep(arrE, lenE, 3, 2);
//	mem::keep(arrO, lenO, 3, 2);
//	PRINT_BOTH ASSERT_GUARDS
//
//	printf("\nScale (cropping) by 2.5:\n");
//	arr::lineSlope1(tempE, lenE);
//	SET_M1
//	mem::scaleCrop(arrE, tempE, lenE, 2.5f);
//	PRINT_EVEN ASSERT_GUARDS
//
//	printf("\nStretch (in-place):\n");
//	SET_LINE
//	mem::stretch(arrE, lenE, 2);
//	PRINT_EVEN ASSERT_GUARDS
//	
//	printf("\nStretch (out-place):\n");
//	arr::lineSlope1(tempE, lenE);
//	mem::zero(arrE, lenE);
//	mem::stretch(arrE, tempE, lenE/2, 2);
//	PRINT_EVEN ASSERT_GUARDS
//
//	printf("\nTranspose2:\n");
//	SET_LINE
//	mem::transpose2(arrE, lenE);
//	PRINT_EVEN ASSERT_GUARDS

}
