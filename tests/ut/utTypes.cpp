{
	typedef Complex<double> Complexd;
	typedef Polar<double> Polard;
	typedef Vec<3,double> Vec3d;

	// Constants
	
	#define T(x, y) assert(x == y);
	T(MaskSign< float>(), 0x80000000UL)
	T(MaskSign<double>(), 0x8000000000000000ULL)
	#undef T


	//{ int x=0; printf("%d\n", (x-9) % 8); }
	//{ unsigned x=0; printf("%ud\n", (x-1)); }
	//{ ValWrap<int> x(7); x=-1; x++; x+=14; x=(x+1)*8; printf("%f\n", x.fraction()); }

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

	{	Multi<3, double> v = {{0,0,0}};
		assert(v == 0);
		assert(v != 1);
		v = 1; assert(v == 1);
	}
	
	{	typedef Vec<3,double> Vec3;
		Vec3 v(0,0,0);
		assert(v == 0);
		v  = 1;						assert(v == 1); assert(v == Vec3(1,1,1));
		v += 1;						assert(v == 2);
		v -= 1;						assert(v == 1);
		v *= 2;						assert(v == 2);
		v /= 2;						assert(v == 1);
		v.set(1,2,3).normalize();	assert(scl::almostEqual(v.mag(),1));
		v.set(1,2,3);				assert(v.sgn() == v.normalize());
	}
}
