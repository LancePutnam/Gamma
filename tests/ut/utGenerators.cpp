{
	const int N = 32;
	RFFT<float> fft(N);
	float buf[N+2];
	std::complex<float> * fdbuf = (std::complex<float> *)buf;
	float * tdbuf = buf+1;
	Sync::master().spu(1);

	/*
	NOTE: Amplitudes of oscillators are multiplied by 2 before being fed to
	FFT so that magnitudes fall in a more intuitive range. The real-to-complex
	FFT actually reports magnitudes that are half of what you would expect.
	*/

	{
		Sine<> g;
		g.freq(1./N);
		for(int i=0;i<N;++i) tdbuf[i]=g()*2;
		//for(int i=0;i<N;++i) printf("% f\n", tdbuf[i]);
		fft.forward(buf, true);
		//for(int i=0;i<N/2+1;++i) printf("[%2d] % f\t% f\n", i, norm(fdbuf[i]), arg(fdbuf[i]));
		assert(near(abs(fdbuf[1]), 1.f));
		assert(near(arg(fdbuf[1]),-M_PI/2));
		
		g.phase(0.25);
		for(int i=0;i<N;++i) tdbuf[i]=g()*2;
		fft.forward(buf, true);
		assert(near(abs(fdbuf[1]), 1.f));
		assert(near(arg(fdbuf[1]), 0.f));		
	}

	{
		double eps=1e-3;
		DSF<> g;
		g.freq(1./N);
		g.harmonics(3);
		g.ampRatio(0.5);
		//g.freqRatio(2);
		for(int i=0;i<N;++i) tdbuf[i]=g()*2;
		fft.forward(buf, true);
		//for(int i=0;i<N/2+1;++i) printf("[%2d] % f\t% f\n", i, abs(fdbuf[i]), arg(fdbuf[i]));
		assert(near(abs(fdbuf[1]), 1, eps));
		assert(near(abs(fdbuf[2]), 0.5, eps));
		assert(near(abs(fdbuf[3]), 0.5*0.5, eps));
		assert(near(abs(fdbuf[4]), 0, eps));
		
		g.freqRatio(2);
		for(int i=0;i<N;++i) tdbuf[i]=g()*2;
		fft.forward(buf, true);
		//for(int i=0;i<N/2+1;++i) printf("[%2d] % f\t% f\n", i, abs(fdbuf[i]), arg(fdbuf[i]));
		assert(near(abs(fdbuf[1]), 1, eps));
		assert(near(abs(fdbuf[3]), 0.5, eps));
		assert(near(abs(fdbuf[5]), 0.5*0.5, eps));
		assert(near(abs(fdbuf[2]), 0, eps));
		assert(near(abs(fdbuf[4]), 0, eps));
		assert(near(abs(fdbuf[6]), 0, eps));
		assert(near(abs(fdbuf[7]), 0, eps));
	}

	{
		{ Player<> p; }
		{ Player<> p; Player<> q(p); }
		{ Array<float> a; Player<> p(a, 1); }
		//{ Player<> p("path/to/soundfile.wav"); }
		
		{
			const double SR = 44100;

			Array<float> a(N);
			Player<float, ipl::Trunc, tap::Clip> p(a, SR);
			
			Sync::master().spu(SR);
			
			assert(p.min() == 0);
			assert(p.max() == N);
			assert(p.rate() == 1);
			assert(p.sampleRate() == SR);
			assert(p.channels() == 1);
			
			assert(p.pos() == 0);
			
			p.advance();
			assert(p.pos() == 1);
			
			for(int i=0; i<N; ++i) p.advance();

			assert(p.pos() == N);			
			
			p.reset();
			assert(p.pos() == 0);
		}
	}
}



