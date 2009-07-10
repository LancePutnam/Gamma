#include <assert.h>
#include "Gamma.h"
#include "Access.h"
#include "Conversion.h"
#include "SmartObject.h"
#include "Types.h"
#include "Thread.h"
#include <map>

using namespace gam;

THREAD_FUNCTION(threadFunc){
	*(int *)user = 1; return NULL;
}


struct TestSmartObject : public SmartObject<TestSmartObject>{};


int main(int argc, char* argv[]){

	// Unit tests are ordered from the least to most dependent functions/objects
	// in order to catch errors in base functionality.


	// Memory management
	{
		auto TestSmartObject a;
		static TestSmartObject s;
		TestSmartObject * d = new TestSmartObject;
		
		assert(!a.dynamicAlloc());
		assert(!s.dynamicAlloc());
		assert(d->dynamicAlloc());

//		struct A{
//			
//			A& operator*(){ printf("*\n"); return *this; }
//			const A& operator*() const { printf("*\n"); return *this; }
//			
//			//const A& operator*() const { printf("* const\n"); return *this; }
//			A* operator&(){ printf("&\n"); return this; }
//			
//			int count;
//		};
//		
//		A as;
//		
//		A * ap = &as;
//		A& ar = *ap;
//		
//		std::list<double> testList;
//		std::map<float, float> testMap;
//		printf("%d\n", sizeof(testList));
//		printf("%d\n", sizeof(testMap));

	}


	// Constants
	
	#define T(x, y) assert(x == y);
	T(MaskSign< float>(), 0x80000000)
	T(MaskSign<double>(), UINT64_C(0x8000000000000000))
	#undef T
	


	// Types

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
		c.fromPolar(1, 0.2);	T(c, Complexd::Polar(1, 0.2))
		c.fromPhase(2.3);		T(c, Complexd::Polar(1, 2.3))
		T(c != Complexd(0,0), true)
		T(c.conj(), Complexd(c.r, -c.i))
		#undef T

		#define T(x, y) assert(scl::almostEqual(x,y,2));
		c.normalize();			T(c.norm(), 1)
		double p=0.1; c(1,0); c *= Complexd::Polar(1, p); T(c.arg(), p)
		#undef T
	}

	{	Quatd q(0,0,0,0);
		#define T(x, y) assert(x == y);
		T(q, Quatd(0,0,0,0))
		T(q.conj(), Quatd(q.r, -q.i, -q.j, -q.k))
		#undef T
	}
	
	
	// Containers
	{
		typedef int t;
		typedef Array<t> array_t;
		array_t * a = new array_t(16);
		array_t * b = new array_t(*a);

		(*a)[0] = 123;
		assert((*a)[0] == 123);
		assert(a->elems() == b->elems());
		assert(array_t::references(a->elems()) == 2);

		delete a;
		assert(array_t::references(b->elems()) == 1);
		
		array_t * c = new array_t(b->elems(), b->size());
		assert(array_t::references(b->elems()) == 2);

		delete b;
		assert(array_t::references(c->elems()) == 1);
		
		t * elemsC = c->elems();
		delete c;
		assert(array_t::references(elemsC) == 0);
		
		a = new array_t(16);
		b = new array_t(*a);
		
		b->own();
		assert(a->elems() != b->elems());
		assert(array_t::references(a->elems()) == 1);
		assert(array_t::references(b->elems()) == 1);
		
		t * elemsA = a->elems();
		t * elemsB = b->elems();
		a->source(*b);
		assert(a->elems() == b->elems());
		assert(array_t::references(elemsA) == 0);
		assert(array_t::references(elemsB) == 2);		
	}
	

	// Conversion

	#define T(x, y) assert(bitsToUInt(x) == y);
	T("", 0) T("0", 0) T("1", 1) T("10", 2) T("11", 3)
	T("0001", 1) T("11111111111111111111111111111111", 0xffffffff)
	T("11111111111111111111111111111110", 0xfffffffe)
	#undef T

	{	Twiddle<float> t(0);
		#define T(x, y) assert(x == y);
		T(t.i, 0) T(t.u, 0) T(t.f, 0)
		t.f = 0.5; T(t.u, bitsToUInt("001111110")<<23)
		t.f = 1.0; T(t.u, bitsToUInt("001111111")<<23)
		t.f = 2.0; T(t.u, bitsToUInt("010000000")<<23)
		t.f =-0.5; T(t.u, bitsToUInt("101111110")<<23)
		t.f =-1.0; T(t.u, bitsToUInt("101111111")<<23)
		t.f =-2.0; T(t.u, bitsToUInt("110000000")<<23)
		#undef T		
	}

	#define T(x, y) assert(castIntRound(x) == y);
	T( 0.0, 0)	T( 0.2, 0) T( 1.0, 1) T( 1.2, 1) T( 1.5, 2) T( 1.8, 2)
				T(-0.2, 0) T(-1.0,-1) T(-1.2,-1) T(-1.5,-2) T(-1.8,-2)
	#undef T

	#define T(x, y) assert(castIntTrunc(x) == y);
	T( 0.0, 0)	T( 0.2, 0) T( 1.0, 1) T( 1.2, 1) T( 1.5, 1) T( 1.8, 1)
				T(-0.2, 0) T(-1.0,-1) T(-1.2,-1) T(-1.5,-1) T(-1.8,-1)
	#undef T

	#define T(x, y) assert(floatExponent(x) == y);
	T(0.125, 124) T(0.25, 125) T(0.50, 126) T(1.00, 127) T(2.00, 128) T(4.00, 129)
	T(0.249, 124) T(0.49, 125) T(0.99, 126) T(1.99, 127) T(3.99, 128) T(7.99, 129)
	T(-0.125, 124) T(-0.25, 125) T(-0.50, 126) T(-1.00, 127) T(-2.00, 128) T(-4.00, 129)
	#undef T

	#define T(x, y) assert(floatMantissa(x) == y);
	T(0.1250, 0.0) T(0.250, 0.0) T(0.50, 0.0) T(1.00, 0.0) T(2.00, 0.0) T(4.00, 0.0)
	T(0.1875, 0.5) T(0.375, 0.5) T(0.75, 0.5) T(1.50, 0.5) T(3.00, 0.5) T(6.00, 0.5)
	#undef T

	#define T(x, y) assert(floatToInt(x) == y);
	T( 0.0, 0)	T( 0.2, 0) T( 1.0, 1) T( 1.2, 1) T( 1.5, 1) T( 1.8, 1)
				T(-0.2, 0) T(-1.0,-1) T(-1.2,-1) T(-1.5,-1) T(-1.8,-1)
	#undef T

	#define T(x, y) assert(floatToUInt(x) == y);
	T( 0.0, 0)	T( 0.2, 0) T( 1.0, 1) T( 1.2, 1) T( 1.5, 1) T( 1.8, 1)
	T(-0.0, 0)	T(-0.2, 0) T(-1.0, 1) T(-1.2, 1) T(-1.5, 1) T(-1.8, 1) 
	#undef T

	#define T(x, y) assert(intToUnit(int16_t(x)) == y);
	T(0, 0) T(-32768, -1) T(32767, 32767./32768)
	#undef T	

	#define T(x, y) assert(uintToUnit<float>(x) == y);
	T(0, 0.0) T(1<<9, 1.1920928955078125e-07) T(0xffffffff, 0.99999988079071044921875)
	 T(1<<29, 0.125) T(1<<30, 0.25) T(1<<31, 0.5)
	#undef T

	

	// Scalar

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
	T(0., 0.) T(0.5, 0.5) T(1., 1.) T(1.2, 1.) T(-0.5, 0.)
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

	#define T(x) assert(scl::abs(scl::invSqrt<1>(x) - 1./sqrt(x)) < 0.002);
	T(0.5f) T(1.f) T(4.f) T(8.f) T(1111.f)
	T(0.50) T(4.0) T(8.0) T(1111.0)
	#undef T

	#define T(x) assert(scl::item(x, 0, 1, 2) == x);
	T(0) T(1) T(2)
	#undef T

	#define T(x) assert(scl::item(x, 0, 1, 2, 3) == x);
	T(0) T(1) T(2) T(3)
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


	// Thread
	{	int x=0;
		Thread t(threadFunc, &x);
		t.wait();
		assert(x == 1);
	}

	return 0;
}
