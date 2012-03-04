{
	using namespace gam::gen;

	const int N = 8;
	double A[N];
	
	#define ASSERT(A, a,b,c,d,e,f,g,h) assert(A[0]==a && A[1]==b && A[2]==c && A[3]==d && A[4]==e && A[5]==f && A[6]==g && A[7]==h);

	slice(A,N) = Val<>(1);
	ASSERT(A, 1,1,1,1,1,1,1,1)
	
	slice(A,N) = gen::Impulse<>();
	ASSERT(A, 1,0,0,0,0,0,0,0)

	slice(A,N) = RCos<>(1./N, 1);
	for(int i=0; i<N; ++i) assert(A[i] = cos(M_2PI*i/N));
	slice(A,N) = RCos<>(2./N, 0.5);
	for(int i=0; i<N; ++i) assert(A[i] = 0.5*cos(2*M_2PI*i/N));

	slice(A,N) = RSin<>(1./N, 1);
	for(int i=0; i<N; ++i) assert(A[i] = cos(M_2PI*i/N));
	slice(A,N) = RSin<>(2./N, 0.5);
	for(int i=0; i<N; ++i) assert(A[i] = 0.5*cos(2*M_2PI*i/N));

	{
		double f = 1./8;
		RSin<> g(f, 1);
		
//		for(int i=0; i<8; ++i){
//			double r = g();
//			printf("% 6.3f, % 6.3f, % 6.3f\n", r, g.phase(), g.amp());
//		}

		for(int i=0; i<8; ++i){
			g();
			assert(scl::almostEqual(g.freq(), f));
			assert(scl::almostEqual(g.amp(), 1));
			//assert(scl::almostEqual(g.phase(), 0) || scl::almostEqual(g.phase(), 1));
		}
	}

	slice(A,N) = RAdd<>(1,0);
	ASSERT(A, 0,1,2,3,4,5,6,7)
	slice(A,N) = RAdd<>(2,1);
	ASSERT(A, 1,3,5,7,9,11,13,15)

	slice(A,N) = RAdd1<>(0);
	ASSERT(A, 0,1,2,3,4,5,6,7)
	slice(A,N) = RAdd1<>(10);
	ASSERT(A, 10,11,12,13,14,15,16,17)

	slice(A,N) = RAddN<1>(0);
	ASSERT(A, 0,1,2,3,4,5,6,7)
	slice(A,N) = RAddN<2>(1);
	ASSERT(A, 1,3,5,7,9,11,13,15)

	slice(A,N) = RAddWrap<>(1,0,4);
	ASSERT(A, 0,1,2,3,0,1,2,3)
	slice(A,N) = RAddWrap<>(3,0,8);
	ASSERT(A, 0,3,6,1,4,7,2,5)
	slice(A,N) = RAddWrap<>(3,10,18,10);
	ASSERT(A, 10,13,16,11,14,17,12,15)
	slice(A,N) = RMul<>(2,1);
	ASSERT(A, 1,2,4,8,16,32,64,128)

	slice(A,N) = RMulAdd<>(2,1,0);
	ASSERT(A, 0,1,3,7,15,31,63,127)
	
	#undef ASSERT
}
