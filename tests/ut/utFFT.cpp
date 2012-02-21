{
	typedef Complex<double> complex;

	struct F{
		static bool aeq(double a, double b, double eps=1e-8){
			return scl::abs(a-b) <= eps;
		}
		static bool aeq(const complex& a, const complex& b, double eps=1e-8){
			return aeq(a.r, b.r, eps) && aeq(a.i, b.i, eps);
		}
	};

	const int N=16;
	complex in[N];
	complex ac[N];
	double ar2[N+2];
	double * ar = ar2+1;
	CFFT<double> cfft(N);
	RFFT<double> rfft(N);

	for(int i=0; i<N; ++i){
		double phs = double(i)/N * M_2PI;
		in[i] = 1. + complex().fromPhase(phs) + complex().fromPhase(2*phs);
		in[i] += (i&1?-1:1); // Nyquist
		ac[i] = in[i];
		ar[i] = ac[i].r;
	}
	
	// test forward transform
	cfft.forward(ac);
	rfft.forward(ar);

//	for(int i=0; i<N; ++i){
//		printf("[%3d] % f % f\n", i, ac[i].r, ac[i].i);
//	}
//	for(int i=0; i<N; ++i){
//		printf("[%3d] % f\n", i, ar[i]);
//	}

	assert(F::aeq(ac[0], complex(1,0)));
	assert(F::aeq(ac[1], complex(1,0)));
	assert(F::aeq(ac[2], complex(1,0)));
	assert(F::aeq(ac[N/2], complex(1,0)));
	assert(F::aeq(ar[0], 1));
	assert(F::aeq(ar[1], 1));
	assert(F::aeq(ar[3], 1));
	assert(F::aeq(ar[N-1], 1));

	// test inverse transform
	cfft.inverse(ac);
	rfft.inverse(ar);

	for(int i=0; i<N; ++i){
		assert(F::aeq(ac[i], in[i]));
		assert(F::aeq(ar[i], in[i].r));
	}
	
	
	// test complex output buffer real-to-complex
	rfft.forward(ar2, true, true);

	assert(F::aeq(ar2[0], 1));
	assert(F::aeq(ar2[1], 0));
	assert(F::aeq(ar2[2], 1));
	assert(F::aeq(ar2[4], 1));
	assert(F::aeq(ar2[N], 1));
	assert(F::aeq(ar2[N+1], 0));

	rfft.inverse(ar2, true);
	for(int i=0; i<N; ++i){
		assert(F::aeq(ar2[i+1], in[i].r));
	}

//	for(int i=0; i<N; ++i) printf("[%3d] % f % f\n", i, ac[i].r, ac[i].i);
//	
//	printf("\n");
////	for(int i=1; i<N-1; i+=2) printf("[%3d] % f % f\n", i, ar[i], ar[i+1]);
//	for(int i=0; i<N; i+=2) printf("% f\n% f ", ar[i], ar[i+1]);

//	{
//		for(int i=0; i<N; ++i){
//			//buf[i+1] = i?0:1;
//			buf[i+1] = 1 + cos(float(i)/N * M_2PI) + (i&1?-1:1);
//			//buf[i+1] = (i&1?-1:1);
//		}
//		fft.forward(buf, true, true);
//		for(int i=0; i<N+2; i+=2) printf("%6.3f %6.3f\n", buf[i], buf[i+1]);
//		fft.inverse(buf, true, true);
//		for(int i=0; i<N; ++i) printf("%g\n", buf[i+1]);
//	}
}
