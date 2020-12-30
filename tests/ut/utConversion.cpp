{
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
	
	#define T(x, y) assert(bits(x) == y);
	T("0",0) T("1",1) T("01",1) T("10",2) T("1111", 15) T("1...", 8)
	#undef T

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
	
	{
	int32_t i;
	#define T(x) assert((split(x, i) == (x - int32_t(x))) && (i == int32_t(x)));
	T(0.f) T(0.1f) T(0.9f) T(1.f) T(1.1f) T(10.7f) T(1000.4512331f) T(16777216.f) T(16777216.5f) T(16777217.f)
	#undef T
	}

	#define T(x, y) assert(uintToUnit<float>(x) == y);
	T(0, 0.0) T(1<<9, 1.1920928955078125e-07) T(0xffffffff, 0.99999988079071044921875)
	 T(1<<29, 0.125) T(1<<30, 0.25) T(1<<31, 0.5)
	#undef T

	#define T(x, y) assert(unitToUInt8(x) == y);
	T(0, 0) T(1./4, 64) T(1./2, 128) T(3./4, 192)
	#undef T

	// from Sample.h
	for(int i=-127; i<=127; ++i){
		auto v = sampleTo<float>(char(i));
		auto j = sampleTo<char>(v);
		//printf("%d -> %g -> %d\n", i,v,j);
		assert(("char-float roundtrip", j==i));
	}
	for(int i=-127; i<=127; ++i){
		auto v = sampleTo<short>(char(i));
		auto j = sampleTo<char>(v);
		//printf("%d -> %d -> %d\n", i,v,j);
		assert(("char-short roundtrip", j==i));
	}
	for(int i=-32767; i<=32767; ++i){
		auto v = sampleTo<float>(short(i));
		auto j = sampleTo<short>(v);
		assert(("short-float roundtrip", j==i));
	}
}
