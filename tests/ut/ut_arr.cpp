{
	const unsigned N=8;
	double A[N], B[N];

	#define IOTA(A) for(unsigned i=0;i<N;++i) A[i]=i
	#define ASSERT(A, a,b,c,d,e,f,g,h) assert(A[0]==a && A[1]==b && A[2]==c && A[3]==d && A[4]==e && A[5]==f && A[6]==g && A[7]==h);
	#define SET(A, a,b,c,d,e,f,g,h) A[0]=a; A[1]=b; A[2]=c; A[3]=d; A[4]=e; A[5]=f; A[6]=g; A[7]=h;
	#define PRINT(A) for(unsigned i=0;i<N;++i) printf("%g ", A[i]); printf("\n")

	IOTA(A);
	IOTA(B);
	assert(
		arr::dot(A,B,N) == (0*0+1*1+2*2+3*3+4*4+5*5+6*6+7*7)
	);

	IOTA(A);
	assert(
		scl::almostEqual(arr::normalize(A,N), 1./(N-1))
	);
	for(unsigned i=0; i<N; ++i)
		assert(scl::almostEqual(A[i], i/double(N-1)));

	#undef IOTA
	#undef ASSERT
	#undef SET
	#undef PRINT
}

//{
//	const unsigned lenE = 8;
//	const unsigned lenO = lenE + 1;
////	const unsigned lenE2 = lenE>>1;
////	const unsigned lenO2 = lenO>>1;
//	
//	float arrE[lenE + 1];
//	float arrO[lenO + 1];
////	float tempE[lenE];
////	float tempO[lenO];
//
//	// Set guard elements
//	arrE[lenE] = 100.f;
//	arrO[lenO] = 100.f;
//
//	printf("\nAdd flipped 1:\n");
//	ZEROS
//	//mem::set(arrE + 1, lenE, 2.f, 2);
//	arr::addFlip(arrE, lenE, 1.f);
//	PRINT_EVEN ASSERT_GUARDS
//
//	printf("\nClip to [1, -1]:\n");
//	printf("Input:\n");
//	arrE[0] = 1.f / FLT_MAX;
//	arrE[1] = -1.f / FLT_MAX;
//	arrE[2] = sqrt(-1.);
//	arrE[3] = -sqrt(-1.);	
//	arrE[4] = 0.f;
//	arrE[5] = -0.f;
//	arrE[6] = 1.f;
//	arrE[7] = -1.f;
//	PRINT_BOTH
//	
//	printf("Output:\n");
//	arr::clip1(arrE, lenE);
//	arr::clip1(arrO, lenO);
//	PRINT_BOTH
//
//	printf("\nDifferentiate:\n");
//	LINE
//	float prev = -1.f; arr::differentiate(arrE, lenE, prev);
//	      prev = -1.f; arr::differentiate(arrO, lenO, prev);
//	PRINT_BOTH ASSERT_GUARDS
//
//	printf("\nExponentiate base 2:\n");
//	LINE
//	arr::expBase(arrE, lenE, 2.f);
//	arr::expBase(arrO, lenO, 2.f);
//	PRINT_BOTH ASSERT_GUARDS
//
//	printf("\nLinear to dB (thresh = -12):\n");
//	printf("Output:\n");
//	arr::linToDB(arrE, lenE, -12.f);
//	PRINT_EVEN
//
//	printf("\nMirror (dp):\n");
//	LINE
//	arr::mirror_dp(arrE, lenE);
//	arr::mirror_dp(arrO, lenO);
//	PRINT_BOTH ASSERT_GUARDS
//
//	printf("\nMirror (dq):\n");
//	LINE
//	arr::mirror_dq(arrE, lenE);
//	arr::mirror_dq(arrO, lenO);
//	PRINT_BOTH ASSERT_GUARDS
//
//	printf("\nMultiply 1s by Bartlett:\n");
//	ONES
//	arr::mulBartlett(arrE, lenE);
//	arr::mulBartlett(arrO, lenO);
//	PRINT_BOTH
//
//	printf("\nNormalize:\n");
//	printf("Input:\n");
//	arrE[0] = 0.f;
//	arrE[1] = 0.1f;
//	arrE[2] = -0.1f;
//	arrE[3] = 0.2f;	
//	arrE[4] = -0.2f;
//	arrE[5] = 0.3f;
//	arrE[6] = -0.3f;
//	arrE[7] = -0.4f;
//	PRINT_EVEN
//	float norm = arr::normalize(arrE, lenE);
//	printf("Output (norm = %f):\n", norm);
//	PRINT_EVEN
//	
//	return 0;
//}

