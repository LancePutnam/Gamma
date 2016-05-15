{
	using namespace gam::scl;

	// Constants
	assert(MaskSign< float>() == 0x80000000UL);
	assert(MaskSign<double>() == 0x8000000000000000ULL);

	// Complex
	{
		typedef Complex<double> Complex;
		typedef Polar<double> Polar;

		// Memory
		Complex c(1,2);
			assert(c == Complex(1,2));
			assert(c != Complex(1,0));
		c = Complex(2,3);
			assert(c == Complex(2,3));

		// Complex operations
		assert(almostEqual(Complex(1,3).conj(), Complex(1,-3)));
		assert(almostEqual(Complex(3,4).mag(), 5));
		assert(almostEqual(Complex( 1, 1).arg(),  45*M_PI/180));
		assert(almostEqual(Complex(-1, 1).arg(), 135*M_PI/180));
		assert(almostEqual(Complex(-1,-1).arg(),-135*M_PI/180));
		assert(almostEqual(Complex( 1,-1).arg(), -45*M_PI/180));
		assert(almostEqual(Complex(1.2,-0.7).normalize().mag(), 1));

		// Polar conversion
		c.fromPolar(1, 0.2);
			assert(almostEqual(c, Polar(0.2)));
		c.fromPhase(2.3);
			assert(almostEqual(c, Polar(2.3)));

		// Basic arithmetic
		c = Complex(1,1);
		Complex b(1,2);
			assert(almostEqual(b+c, Complex( 2, 3)));
			assert(almostEqual(b+1.,Complex( 2, 2)));
			assert(almostEqual(b-c, Complex( 0, 1)));
			assert(almostEqual(b-1.,Complex( 0, 2)));
		double p=0.1;
		c = Complex(1, 0);
		c *= Polar(1, p);
			assert(almostEqual(c.arg(), p));
			assert(almostEqual(b*2.,Complex( 2, 4)));
			assert(almostEqual(2.*b,Complex( 2, 4)));
			assert(almostEqual(b/2.,Complex(0.5,1)));
			assert(almostEqual(2./b,Complex(1 * 2./5, -2 * 2./5)));

		c.fromPolar(4, 0.2);
			assert(almostEqual(sqrt(c).norm(), 2));
			assert(almostEqual(sqrt(c).arg(), 0.1));
	}

	// Vec
	{	typedef Vec<3,double> Vec3;
		Vec3 v(0,0,0);
		assert(v == 0);
		v  = 1;						assert(v == 1); assert(v == Vec3(1,1,1));
		v += 1;						assert(v == 2);
		v -= 1;						assert(v == 1);
		v *= 2;						assert(v == 2);
		v /= 2;						assert(v == 1);
		v.set(1,2,3).normalize();	assert(scl::almostEqual(v.mag(),1));
		v.set(1,2,3);				assert(v.normalized() == v.normalize());
	}
}
