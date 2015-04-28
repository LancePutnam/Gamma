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
		
		//in[i] = (i==0) ? 1 : 0;
		
		ac[i] = in[i];
		ar[i] = ac[i].r;
	}
	
	// test forward transform
	cfft.forward(ac);
	rfft.forward(ar);

	//for(int i=0; i<N; ++i) printf("[%3d] % f % f\n", i, ac[i].r, ac[i].i); printf("\n");
	//for(int i=0; i<N; ++i) printf("[%3d] % f\n", i, ar[i]); printf("\n");

	assert(F::aeq(ac[  0], complex(1,0)));
	assert(F::aeq(ac[  1], complex(1,0)));
	assert(F::aeq(ac[  2], complex(1,0)));
	assert(F::aeq(ac[N/2], complex(1,0)));
	assert(F::aeq(ar[  0], 1.0));
	assert(F::aeq(ar[  1], 0.5));
	assert(F::aeq(ar[  3], 0.5));
	assert(F::aeq(ar[N-1], 1.0));

	// test inverse transform
	cfft.inverse(ac);
	rfft.inverse(ar);

	//for(int i=0; i<N; ++i) printf("[%3d] % f\n", i, ar[i]); printf("\n");

	for(int i=0; i<N; ++i){
		assert(F::aeq(ac[i], in[i]));
		assert(F::aeq(ar[i], in[i].r));
	}
	
	
	// test complex output buffer real-to-complex
	rfft.forward(ar2, true);

	assert(F::aeq(ar2[  0], 1.0));
	assert(F::aeq(ar2[  1], 0. ));
	assert(F::aeq(ar2[  2], 0.5));
	assert(F::aeq(ar2[  3], 0. ));
	assert(F::aeq(ar2[  4], 0.5));
	assert(F::aeq(ar2[  N], 1. ));
	assert(F::aeq(ar2[N+1], 0.0));

	rfft.inverse(ar2, true);
	for(int i=0; i<N; ++i){
		assert(F::aeq(ar2[i+1], in[i].r));
	}
}


{
	for(int j=0; j<3; ++j){
		const int N = 32;
		SpectralType specType[] = {COMPLEX, MAG_PHASE, MAG_FREQ};

		STFT stft(N, N/4, 0, HANN, specType[j]);

		// configure to produce most precise output
		stft.inverseWindowing(false);
		stft.precise(true);

		float t=0; // output

		for(int i=0; i<N*4; ++i){
			float p = float(i)/N * 2*M_PI;
			float s = cos(p);
			if(stft(s)){
				//printf("\n");
			}

			// output should match input after N-1 samples
			if(i>=N){
				//printf("[%3d] %f %f\n", i, s,t);
				assert(near(s,t, 2e-6));
			}

			t = stft();
			//printf("[%3d] ", i); printPlot(t); printf("\n");
		}
	}
}

