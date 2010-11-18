{
	{
		using namespace acc;

		assert(None::map(0,2,1) == 0);
		assert(Wrap::map(0,2,1) == 2);
		assert(Clip::map(0,2,1) == 1);

		assert(None::map(2,2,1) == 2);
		assert(Wrap::map(2,2,1) == 2);
		assert(Clip::map(2,2,1) == 2);

		assert(None::map(3,2,1) == 3);
		assert(Wrap::map(3,2,1) == 1);
		assert(Clip::map(3,2,1) == 2);
	}
	
	const uint32_t N=8;
	double A[N], B[N];
	Slice<double> a(A,N), b(B,N);

	#define ASSERT(A, a,b,c,d,e,f,g,h) assert(A[0]==a && A[1]==b && A[2]==c && A[3]==d && A[4]==e && A[5]==f && A[6]==g && A[7]==h)
	#define SET(A, a,b,c,d,e,f,g,h) A[0]=a; A[1]=b; A[2]=c; A[3]=d; A[4]=e; A[5]=f; A[6]=g; A[7]=h
	#define PRINT(A) for(uint32_t i=0;i<N;++i) printf("%g ", A[i]); printf("\n")

//	{
//		Indexer ind(8);
//		for(int i=0; i<8; ++i){ A[i]=ind(i,8); }	ASSERT(A, 0,1,2,3,4,5,6,7);
//	}
//
//	{
//		IndexerInt ind(8);
//		for(int i=0; i<8; ++i){ A[i]=ind(i,8); }	ASSERT(A, 0,1,2,3,4,5,6,7);
//
//		ind.stride(2);
//		for(int i=0; i<8; ++i){ A[i]=ind(i,8); }	ASSERT(A, 0,2,4,6,8,10,12,14);
//
//		ind.stride(-1);
//		for(int i=0; i<8; ++i){ A[i]=ind(i,8); }	ASSERT(A, 0,-1,-2,-3,-4,-5,-6,-7);
//
//		ind.stride(-2);
//		for(int i=0; i<8; ++i){ A[i]=ind(i,8); }	ASSERT(A, 0,-2,-4,-6,-8,-10,-12,-14);
//	}
//
//	{
//		IndexerReal ind(8);
//		for(int i=0; i<8; ++i){ A[i]=ind(i,8); }	ASSERT(A, 0,1,2,3,4,5,6,7);
//
//		ind.stride(2);
//		for(int i=0; i<8; ++i){ A[i]=ind(i,8); }	ASSERT(A, 0,2,4,6,0,2,4,6);
//
//		ind.stride(-1);
//		for(int i=0; i<8; ++i){ A[i]=ind(i,8); }	ASSERT(A, 0,7,6,5,4,3,2,1);
//
//		//ind.stride(4. + 1/2.);
//		ind.strides(4, 2);
//		for(int i=0; i<8; ++i){ A[i]=ind(i,8); }	ASSERT(A, 0,4,1,5,2,6,3,7);
//
//		//ind.stride(2. + 1/4.);
//		ind.strides(2, 4);
//		for(int i=0; i<8; ++i){ A[i]=ind(i,8); }	ASSERT(A, 0,2,4,6,1,3,5,7);
//
//		ind.stride(0.5);
//		for(int i=0; i<8; ++i){ A[i]=ind(i,8); }	ASSERT(A, 0,0,1,1,2,2,3,3);
//
//		// overstrided
//		ind.stride(1. + 32.);
//		for(int i=0; i<8; ++i){ A[i]=ind(i,8); }	ASSERT(A, 0,1,2,3,4,5,6,7);
//
//		ind.stride(-1. + 32.);
//		for(int i=0; i<8; ++i){ A[i]=ind(i,8); }	ASSERT(A, 0,7,6,5,4,3,2,1);
//			
//		//for(int32_t i=1; i<32678*; ++i){ /*printf("%d ", i);*/ assert(1 == int32_t(i * (1./i) + 1e-8)); }
//		//ind.strides(2, 2, 2); 
//		//ind.stride(2 + 2./7); for(int i=0; i<8; ++i){ printf("%d ", ind(i,8)); } printf("\n");		
//	}


	// Basic slice functionality

	// Random access
	SET(A,0,1,2,3,4,5,6,7);	for(int i=0;i<8;++i) assert(a[i] == i);
	
	// Setting elements
	a.set();						assert(a == double());
	a.set(1.);						assert(a == 1.);
	a = 0.;							assert(a == 0.);
	a(4,2) = 1.;					assert(a(4,2) == 1.);
	a(4,2,1) = 2.;					assert(a(4,2,1) == 2.);
	a(4,-2,-1) = 3.;				assert(a(4,-2,-1) == 3.);
	a(4,-2,-2) = 4.;				assert(a(4,-2,-2) == 4.);
	
	// Copying
	a=0.; b=1.;
	a.copy(b);						assert(a == 1.);

	SET(B,1,2,3,4,5,6,7,8);
	a.copy(b(8,-1,-1));				ASSERT(A,8,7,6,5,4,3,2,1);
	a = 0.;
	a.copy(b.reversed());			ASSERT(A,8,7,6,5,4,3,2,1);
	a.swap(b);						ASSERT(A, 1,2,3,4,5,6,7,8);
									ASSERT(B, 8,7,6,5,4,3,2,1);
	a(4).swap(a(4,1,4).reverse());	ASSERT(A, 8,7,6,5,4,3,2,1);
	
	// Arithmetic operations
	a = 0.;
	a += 1.;						assert(a == 1.);
	a += (b=1.);					assert(a == 2.);
	a += b(4);						ASSERT(A, 3,3,3,3,2,2,2,2);

	a = 2.; b = 2.;
	a += b;							assert(a == 4.);
	a *= b;							assert(a == 8.);
	a += a;							assert(a == 16.);
	a *= a;							assert(a == 256.);

	a = 1.;
	assert(a.sum() == 8);

	{
	double A[4];
	double B[4] = {1,2,3,4};
	Slice<double> a(A,4), b(B,4);

	a.copy(b);			// A = {1,2,3,4}
	b = 0.;				// B = {0,0,0,0}
	a.swap(b);			// A = {0,0,0,0}, B = {1,2,3,4}
	}

	
//	SET(B,1,2,3,4,5,6,7,8)
//	a = 0.;
//	
////	0 1 2 3 4 5 6 7
////	0 4 1 5 2 6 3 7
//	
//	for(int i=0; i<N; ++i){
//		int j = scl::wrap((2. + 1/4.)*i, 8.);
//		A[i] = B[j];
//	}
//	
//	//ASSERT(A,1,5,2,6,3,7,4,8)
	#undef ASSERT
	#undef SET
	#undef PRINT
}

/*
		{
		
			typedef std::map<uint32_t, std::vector<double> > MyMap;
			MyMap m;
		
			int K=1024*32;
			for(int j=0; j<K; ++j){
				double s = double(j)/K * 8.;
				ind.stride(s);
				for(int i=0; i<8; ++i){ A[i] = ind(i,8); }
				
				if(mem::unique(A,N)){
				
					uint32_t h =0;
					for(int i=0; i<8; ++i) h += pow(10,7-i)*A[i];
				
					m[h].push_back(s);
					//for(int i=0; i<8; ++i){ printf("%g ", A[i]); } printf("    %g\n", s);
				}
			}

			MyMap::iterator it;
			// show content:
			printf("%d unique sequences\n\n", m.size());
			for( it=m.begin() ; it != m.end(); it++ ){
				uint32_t h=(*it).first;
				
				printf("0%d  ", h);
				
				const std::vector<double>& v = (*it).second;
				for(unsigned i=0; i<v.size(); ++i){
					//printf("%g ", v[i]);
				}
				
				printf("\n");
			}
		}
		
*/


