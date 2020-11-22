{
	const int N = 32;
	RFFT<float> fft(N);
	float buf[N+2];
	std::complex<float> * fdbuf = (std::complex<float> *)buf;
	float * tdbuf = buf+1;
	gam::sampleRate(1);

	{
		Accum<> g;
		g.freq(1./4);
		g.phase(0);
		assert(g.nextPhase() == 0);
		assert(g.nextPhase() == 1UL<<30);
		assert(g.nextPhase() == 1UL<<31);

		g.phaseMax();
		assert(g() != 0);
		assert(g() == 0);
	}

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
		{
			SamplePlayer<> p;
			assert(!p.valid());
		}

		{ SamplePlayer<> p; SamplePlayer<> q(p); }
		{ Array<float> a; SamplePlayer<> p(a, 1); }
		//{ SamplePlayer<> p("path/to/soundfile.wav"); }
		
		{
			const double SR = 44100;

			Array<float> a(N);
			SamplePlayer<float, ipl::Trunc, tap::Clip> p(a, SR);
			
			gam::sampleRate(SR);
			
			assert(p.min() == 0);
			assert(p.max() == N);
			assert(p.rate() == 1);
			assert(p.frameRate() == SR);
			assert(p.channels() == 1);
			
			assert(p.pos() == 0);
			
			p.advance();
			assert(p.pos() == 1);
			
			for(int i=0; i<N; ++i) p.advance();

			assert(p.pos() == (N-1));			
			
			p.reset();
			assert(p.pos() == 0);
		}

		{ // Multichannel
			SamplePlayer<float> p;
			static const int frames = 3;
			static const int chans = 2;
			float SR = 1;
			gam::sampleRate(SR);
			float bufI[frames*chans] = {1.,10., 2.,20., 3.,30.};
			p.buffer(bufI, frames, SR, chans, true /*interleaved*/);
			for(int i=1; i<=frames; ++i){
				assert(p.read(0) == i);
				assert(p.read(1) == i*10);
				p.advance();
			}
			//printf("%g %g\n", p.read(0), p.read(1));
			float bufD[frames*chans] = {1.,2.,3., 10.,20.,30.};
			p.buffer(bufD, frames, SR, chans, false /*interleaved*/);
			for(int i=1; i<=frames; ++i){
				assert(p.read(0) == i);
				assert(p.read(1) == i*10);
				p.advance();
			}
		}
	}
}



