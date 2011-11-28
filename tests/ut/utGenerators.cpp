{
	{
		{ Player<> p; }
		{ Player<> p; Player<> q(p); }
		{ Array<float> a; Player<> p(a, 1); }
		//{ Player<> p("path/to/soundfile.wav"); }
		
		{
			const int N = 128;
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



