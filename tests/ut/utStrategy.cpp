{
	static const int N = 16;
	typedef float sample_t;
	sample_t src[N];
	for(int i=0; i<N; ++i) src[i] = i;

	{	using namespace gam::ipl;
		assert(scl::almostEqual( Trunc<sample_t>()(src, 1, 0.5, N-1, 0), sample_t(1.)));
		assert(scl::almostEqual( Round<sample_t>()(src, 1, 0.5, N-1, 0), sample_t(2.)));
		assert(scl::almostEqual(Linear<sample_t>()(src, 1, 0.5, N-1, 0), sample_t(1.5)));
	}

	{	using namespace gam::iplSeq;
		{	Base<3,float> I;
			I.push(1);		assert(1 == I.val());
			I.push(2);		assert(2 == I.val());
			I.push(3);		assert(3 == I.val());
			I.push(4);		assert(4 == I.val());
		}
		{	None<float> I;
			I.push(1);
			I.push(2);
			assert(2 == I(0));
			assert(2 == I(0.5));
		}
		{	Trunc<float> I;
			I.push(1);
			I.push(2);
			assert(1 == I(0));
			assert(1 == I(0.5));
		}
		{	Round<float> I;
			I.push(1);
			I.push(2);
			assert(1 == I(0));
			assert(2 == I(0.5));
		}
		{	Linear<float> I;
			I.push(1);
			I.push(2);
			assert(1   == I(0));
			assert(scl::almostEqual(sample_t(1.5), I(0.5)));
		}
	}

	{	using namespace gam::phsInc;
		assert(incClip(     0.f, 1.f, 500.f,0.f) == 1.f);
		assert(incClip(-900.f, 1.f, 500.f,0.f) == 0.f); // below min
		{ auto x=incClip(900.f, 1.f, 500.f,0.f); assert(499.f < x && x < 500.f); } // above max
	}
}
