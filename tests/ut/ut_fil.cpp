{
	const int N = 8;
	double ds = 1./N;
	double eps = 1e-12;
	double A[N];
	
	#define RESET for(int i=0;i<N;++i) A[i]=i+1;

	#define ASSERT(A, a,b,c,d,e,f,g,h) assert(A[0]==a && A[1]==b && A[2]==c && A[3]==d && A[4]==e && A[5]==f && A[6]==g && A[7]==h);

	
	//fil::LCCD<0,1,double,double> lccd;


	{
		Slice<double> s(A,N);
		s = 1.;
		s.filter(fil::Integrator<>());
		ASSERT(A, 1,2,3,4,5,6,7,8);
		
		s.filter(fil::Delay1<>());
		ASSERT(A, 0,1,2,3,4,5,6,7);

		s.filter(fil::Delay2<>(6,7));
		ASSERT(A, 6,7,0,1,2,3,4,5);
		
		RESET
		s.filter(fil::Hold<>(4));
		ASSERT(A, 1,1,1,1,5,5,5,5);
		
	}


	fil::Reson<double> res(ds, 1, 0, 1);
	
	assert(res.freq() == ds);
	assert(res.decay() == 1);
	assert(res.ahead() == 1);
	assert(scl::almostEqual(res.behind().i, sin(-2*ds*M_2PI)));
	for(int i=0; i<N; ++i) assert(scl::abs(sin(i*ds*M_2PI) - res().i) < eps);

	res.set(0, 1, 0, 1);
	res.decay(0.5, N);
	for(int i=0; i<N; ++i) res();
	assert(scl::almostEqual(res.r, 0.5) );

	res.set(ds, 0, 0, 1);
	for(int i=0; i<N; ++i) assert(scl::abs(sin(i*ds*M_2PI) - res(i?0:1).i) < eps);

//	fil::TransferFunc fr;
//	//fr.addX(0.1, 0).addY( 0.9, 1);	// one-pole lo-pass
//	//fr.addX(0.1, 0).addY(-0.9, 1);	// one-pole hi-pass
//	//fr.addX(0.5, 0).addX( 0.5, 1);	// one-zero lo-pass
//	//fr.addX(0.5, 0).addX(-0.5, 1);	// one-zero hi-pass
//	//fr.addX(0.25,0).addX(-0.5, 1).addY(0.5, 1).gain(2); // 1st order all-pass
//	fr.addX(0.25,0).addX(-0.5, 4).addY(0.5, 4).gain(2); // 4th order all-pass comb
//	//fr.addX(0.5, 0).addX( 0.5, 8);	// 4-notch lo-pass comb
//	//fr.addX(0.5, 0).addX(-0.5, 8);	// 4-notch hi-pass comb
//	//fr.addX(0.2, 0).addY( 0.8, 8);	// 4-peak lo-pass comb
//	//fr.addX(0.2, 0).addY(-0.8, 8);	// 4-peak hi-pass comb
//
//	for(int i=0; i<32; ++i){
//		float f = float(i)/32 * 0.5;
//		Complexd::Polar p = fr(f);
//		printf("% 6.2f % 6.2f ", p.m, p.p*M_1_PI);
//		printPlot(p.m);
//		printPlot(p.p*M_1_PI);
//		printf("\n");
//	}

//	{
//		//fil::FunctionObject1<float, scl::clipS<float> > f(1);
//		fil::FunctionObject1<float, scl::clipS<float> > f(1);
//		assert(f(10) == 1); assert(f(-10) == -1);
//		
//		fil::Function2<float, float, float> f2(scl::clipS<float>, 10, 1);
//		
//		assert(f2() == 1);
//		assert(f2(-10) == -1);
//	}
	#undef RESET
	#undef ASSERT
}
