{
	const double pinf = INFINITY;		// + infinity
	const double ninf = -INFINITY;		// - infinity
//	const double nan  = pinf * 0.;		// not-a-number

	#define T(x, y) assert(scl::abs(x) == y);
	T(1, 1) T(-1, 1) T(1.f, 1.f) T(-1.f, 1.f) T(1., 1.) T(-1., 1.)
	#undef T

	#define T(x, y) assert(scl::almostEqual(x,y,10));
	T(0.f, 0.f) T(1.f, 0.9999999f) T(1.f, 1.000001f)
	T(0.0, 0.0) T(1.0, 0.999999999999999) T(1.0, 1.000000000000001)
	#undef T

	#define T(x, y) assert(scl::ceil(x) == y);
	T(0., 0.)	T( 1., 1.) T( 1.2, 2.) T( 1.8, 2.) T( 1000.1, 1001.)
				T(-1.,-1.) T(-1.2,-1.) T(-1.8,-1.) T(-1000.1,-1000.)
	#undef T

	#define T(x, y) assert(scl::ceilEven(x) == y);
	T(0, 0) T(1, 2) T(2, 2) T(3, 4) T(1001, 1002)
	#undef T

	#define T(x, y) assert(scl::ceilPow2(x) == y);
	T(0, 0) T(1, 1) T(2, 2) T(3, 4)
	T(500, 512) T(999, 1024)
	#undef T

	#define T(x, y) assert(scl::clip(x) == y);
	T(0., 0.) T(0.5, 0.5) T(1., 1.) T(1.2, 1.) T(-0.5, 0.) T(pinf, 1.) T(ninf, 0.)
	#undef T

	#define T(x, y) assert(scl::clipS(x) == y);
	T(0., 0.) T(0.5, 0.5) T(1., 1.) T(1.2, 1.) T(-0.5, -0.5) T(-1., -1) T(-1.2, -1.)
	#undef T

	#define T(x) assert(scl::abs(scl::cosP3(x) - cos(x*M_2PI)) < 0.018);
	T(0.) T(0.1) T(0.2) T(0.3) T(0.4) T(0.5)
	#undef T
	
	#define T(x, y) assert(scl::factorial12(x) == y);
	T(0, 1) T(1, 1) T(2, 2*1) T(3, 3*2*1) T(4, 4*3*2*1)
	T(5, 5*4*3*2*1) T(6, 6*5*4*3*2*1) T(7, 7*6*5*4*3*2*1) T(8, 8*7*6*5*4*3*2*1)
	T(9, 9*8*7*6*5*4*3*2*1) T(10, 10*9*8*7*6*5*4*3*2*1) T(11, 11*10*9*8*7*6*5*4*3*2*1)
	T(12, 12*11*10*9*8*7*6*5*4*3*2*1)
	#undef T

	#define T(x, y) assert(scl::floor(x) == y);
	T(0., 0.)	T( 1., 1.) T( 1.2, 1.) T( 1.8, 1.) T( 1000.1, 1000.)
				T(-1.,-1.) T(-1.2,-2.) T(-1.8,-2.) T(-1000.1,-1001.)
	#undef T

	#define T(x, y) assert(scl::floorPow2(x) == y);
	T(0, 1) T(1, 1) T(2, 2) T(3, 2)
	T(513, 512) T(1090, 1024)
	#undef T

	#define T(x, y) assert(scl::almostEqual(scl::fold(x), y));
	T(0., 0.) T(0.5, 0.5) T(1., 1.) T(1.2, 0.8) T(-0.2, 0.2)
	T(2.2, 0.2) T(3.2, 0.8) T(4.2, 0.2) T(5.2, 0.8)
	#undef T
	
	#define T(x,y,r) assert(scl::gcd(x,y) == r);
	T(7,7,7) T(7,4,1) T(8,4,4)
	#undef T

	#define T(x) assert(scl::abs(scl::invSqrt<1>(x) - 1./sqrt(x)) < 0.002);
	T(0.5f) T(1.f) T(4.f) T(8.f) T(1111.f)
	T(0.50) T(4.0) T(8.0) T(1111.0)
	#undef T

	#define T(x,y,r) assert(scl::lcm(x,y)==r);
	T(7,3,21) T(8,4,8) T(3,1,3) //T(0,0,0)
	#undef T

	#define T(x) assert(scl::log2(1<<x) == x);
	T(0) T(1) T(2) T(3) T(4) T(29) T(30) T(31)
	#undef T

	#define T(x) assert(scl::almostEqual(scl::log2Fast(x), log2f(x), 3800000));
	T(1) T(2) T(8) T(128) T(1.2) T(1.4) T(0.2)
	#undef T

	#define T(x,y) assert(scl::almostEqual(scl::mapDepth(x, 0.5), y));
	T(-1., 0.5) T(0., 0.75) T(1., 1.)
	#undef T

	#define T(x,y) assert(scl::almostEqual(scl::mapLin(x, -1., 1., 0., 1.), y));
	T(-1., 0.) T(0., 0.5) T(1., 1.)
	#undef T

	#define T(x) assert(scl::pow2(x) == x*x);
	T(0) T(1) T(2) T(3) T(-1) T(-2) T(-3)
	#undef T

	#define T(x) assert(scl::pow2S(x) == x*scl::abs(x));
	T(0) T(1) T(2) T(3) T(-1) T(-2) T(-3)
	#undef T

	#define T(x) assert(scl::pow3(x) == x*x*x);
	T(0) T(1) T(2) T(3) T(-1) T(-2) T(-3)
	#undef T

	#define T(x) assert(scl::pow3Abs(x) == scl::abs(x*x*x));
	T(0) T(1) T(2) T(3) T(-1) T(-2) T(-3)
	#undef T

	#define T(x) assert(scl::pow4(x) == x*x*x*x);
	T(0) T(1) T(2) T(3) T(-1) T(-2) T(-3)
	#undef T

	#define T(x) assert(scl::pow5(x) == x*x*x*x*x);
	T(0) T(1) T(2) T(3) T(-1) T(-2) T(-3)
	#undef T
	
	#define T(x) assert(scl::pow6(x) == x*x*x*x*x*x);
	T(0) T(1) T(2) T(3) T(-1) T(-2) T(-3)
	#undef T

	#define T(x) assert(scl::pow8(x) == x*x*x*x*x*x*x*x);
	T(0) T(1) T(2) T(3) T(-1) T(-2) T(-3)
	#undef T

	#define T(x) assert(scl::pow16(x) == x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x);
	T(0) T(1) T(2) T(3) T(-1) T(-2) T(-3)
	#undef T

	#define T(x) assert(scl::almostEqual(scl::pow64(x), x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x, 10));
	T(0.) T(1.) T(1.01) T(1.02) T(-1.) T(-1.01) T(-1.02)
	#undef T
	
	#define T(x,y,r) assert(scl::remainder(x,y) == r);
	T(7,7,0) T(7,1,0) T(7,4,3) T(7,3,1) T(14,3,2)
	#undef T

	#define T(x,y) assert(scl::round(x) == y);
	T(0.f, 0.f) T(0.2f, 0.f) T(0.8f, 1.f) T(-0.2f, 0.f) T(-0.8f,-1.f) T(0.5f, 0.f) T(-0.5f, 0.f)
	T(0.0, 0.0) T(0.20, 0.0) T(0.80, 1.0) T(-0.20, 0.0) T(-0.80,-1.0) T(0.50, 0.0) T(-0.50, 0.0)
	#undef T

	#define T(x) assert(scl::abs(scl::sqrt<1>(x) - sqrt(x)) < 0.008);
	T(0.f) T(1.f) T(4.f) T(8.f) T(1111.f)
	T(0.0) T(1.0) T(4.0) T(8.0) T(1111.0)
	#undef T

	#define T(x,y) assert(scl::sumOfSquares(x) == y);
	T(1., 1.) T(2., 1*1+2*2) T(3., 1*1+2*2+3*3) T(4., 1*1+2*2+3*3+4*4) T(5., 1*1+2*2+3*3+4*4+5*5)
	#undef T

	#define T(x,y) assert(scl::trunc(x) == y);
	T(0.f, 0.f) T(0.2f, 0.f) T(0.8f, 0.f) T(-0.2f, 0.f) T(-0.8f, 0.f) T(0.5f, 0.f) T(-0.5f, 0.f)
	T(0.0, 0.0) T(0.20, 0.0) T(0.80, 0.0) T(-0.20, 0.0) T(-0.80, 0.0) T(0.50, 0.0) T(-0.50, 0.0)
	#undef T

	#define T(x, y) assert(scl::almostEqual(scl::wrap(x, 1., -1.), y));
	T(0., 0.)	T( 0.5, 0.5) T( 1.,-1.) T( 1.2,-0.8) T( 2.2, 0.2)
				T(-0.5,-0.5) T(-1.,-1.) T(-1.2, 0.8) T(-2.2,-0.2)
	#undef T
}
