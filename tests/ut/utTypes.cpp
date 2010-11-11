{
	typedef Complex<double> Complexd;
	typedef Polar<double> Polard;
	typedef Quat<double> Quatd;
	typedef Vec3<double> Vec3d;

	// Constants
	
	#define T(x, y) assert(x == y);
	T(MaskSign< float>(), 0x80000000)
	T(MaskSign<double>(), UINT64_C(0x8000000000000000))
	#undef T


	//{ int x=0; printf("%d\n", (x-9) % 8); }
	//{ unsigned x=0; printf("%ud\n", (x-1)); }
	//{ ValWrap<int> x(7); x=-1; x++; x+=14; x=(x+1)*8; printf("%f\n", x.fraction()); }

	{	Bits<> b;
		#define T(x, y) assert(x == y);
		b.enable(1<<0);				T(b.enabled(1<<0), 1<<0)
		b.enable(1<<1);				T(b.enabled(1<<1), 1<<1)
		b.toggle(1<<0);				T(b.enabled(1<<0), 0) T(b.enabled(1<<1), 1<<1)
		b.set(1<<1, false);			T(b.enabled(1<<1), 0)
		b.set(b.mask(0,1), true);	T(b.enabled(1<<0), 1<<0) T(b.enabled(1<<1), 1<<1)
		b.disable(1<<0);			T(b.enabled(1<<0), 0) T(b.enabled(1<<1), 1<<1)
		b.zero();					T(b(), 0)
		#undef T
	}

	{	Complexd c(0,0);
		#define T(x, y) assert(x == y);
		T(c, Complexd(0,0))
		c.fromPolar(1, 0.2);	T(c, Polard(0.2))
		c.fromPhase(2.3);		T(c, Polard(2.3))
		T((c != Complexd(0,0)), true)
		T(c.conj(), Complexd(c.r, -c.i))
		#undef T

		#define T(x, y) assert(scl::almostEqual(x,y,2));
		c.normalize();			T(c.norm(), 1)
		double p=0.1; c(1,0); c *= Polard(1, p); T(c.arg(), p)

		c.fromPolar(4,0.2);
		T(c.sqrt().norm(), 2)
		T(c.sqrt().arg(), 0.1)
		#undef T
	}

	{	Quatd q(0,0,0,0);
		#define T(x, y) assert(x == y);
		T(q, Quatd(0,0,0,0))
		T(q.conj(), Quatd(q.r, -q.i, -q.j, -q.k))
		#undef T
	}

	{	Multi<3, double> v = {{0,0,0}};
		assert(v == 0);
		assert(v != 1);
		v = 1; assert(v == 1);
	}
	
	{	ShiftBuffer<2, double> v;
		assert(v == 0);
		v(1);		assert(v[0]==1); assert(v[1]==0);
		v(2);		assert(v[0]==2); assert(v[1]==1);
	}
	
	{	Vec3d v(0,0,0);
		assert(v == 0);
		v  = 1;		assert(v == 1); assert(v == Vec3d(1,1,1));
		v += 1;		assert(v == 2);
		v -= 1;		assert(v == 1);
		v *= 2;		assert(v == 2);
		v /= 2;		assert(v == 1);
		v(1,2,3).normalize(); assert(scl::almostEqual(v.norm(),1));
		v(1,2,3);	assert(v.sgn() == v.normalize());
	}
}
