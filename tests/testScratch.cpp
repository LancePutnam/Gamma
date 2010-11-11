#include "Gamma/Gamma.h"
//#include "../doc/scratch.h"
#include "Gamma/Effects.h"

#define ONES	mem::set(buf, gen::Val<>(1), Loop(size));
#define ZEROS	mem::zero(buf, size);
#define PRINT	for(uint32_t i=0; i<size; i++) printf("% 5.2f ", buf[i]); printf("\n");

int main(int argc, char** argv){
/*
	using namespace gam;

	const uint32_t size = 8;
	float buf[size];
	Sync::master().spu(size);
	
	{	using namespace mem;
		using namespace scl;

		float temp[16];
		arr::lineSlope1(temp, 16);
		//call(pow3<float>, temp, temp, 16);

		#define MANY(name, o0, i0, len)\
			for(int i=0; i<len; ++i) o0[i] = name(i0[i]);
			
		MANY(pow3, temp, temp, 16);
		
		arr::print(temp, 16);
	}
	
	#define DO(fnc)\
		printf("\n"#fnc":\n"); fnc;\
		for(int i=0; i<multi.size(); ++i) printf("% 6.2f ", multi[i]); printf("\n");
		
	Val2<float> multi(1, 2);
	DO(multi)
	DO(multi *= 2)
	DO(multi *= multi)
	DO(multi += multi)
	DO(multi = multi + multi)
	DO(multi += 1)
	DO(multi = multi + 1)
	DO(multi = -multi)
	DO(multi = 0)

	gen::RAdd< Val2<float> > addGen( Val2<float>(0.1, 0.2), 0 );
	
	for(int i=0; i<16; ++i){
		Val2<float> v = addGen();
		printf("%f %f\n", v[0], v[1]);
	}
	
	#undef DO
	
	{	
		#define DO(fnc) printf("\n"#fnc"\n\ti "); PRINT fnc; printf("\to "); PRINT
		
		using namespace scratch;
		using namespace gam::arr;
		using namespace gam::gen;
		using namespace gam::mem;

		ONES  DO( add (buf, buf, size) )
		ONES  DO( add (buf, val(1), size) )
		ONES  DO( add (buf, buf, Loop(size, 2)) )
		ONES  DO( add (buf, val(1), Loop(size, 2)) )
		ZEROS DO( add (buf, rAdd(1), size) )
		ZEROS DO( add	(buf, rAdd(1), Loop(size, 2)) )
		ONES  DO( invert (buf, buf, size) )
		ONES  DO( invert (buf, buf, Loop(size,2)) )
		ONES  DO( invert (buf, buf, Loop(size,2,1)) )
		ONES  DO( mul (buf, val(2), size) )
		ONES  DO( mul (buf, buf, val(2), size) )
		ONES  DO( sub (buf, buf, size) )
		ONES  DO( sub (buf, val(3), size) )
		ONES  DO( sub (buf, val(3), buf, size) )


		mem::set(buf, RAdd<>(4./size, -2), Loop(size));
		DO( filter(buf, Clip(1.f), size) )

		mem::set(buf, RAdd<>(5./size, -2), Loop(size));
		DO( filter(buf, Wrap(1.f), size) )
		
		mem::set(buf, RAdd<>(4./size, -2), Loop(size));
		DO( filter(buf, Pow2, size) )
		
		
		
		ONES  DO( add(buf, buf, size))
		
		Indices indices(size); indices << 1 << 2 << 4 << 5;
		ONES  DO( add(buf, buf, indices))
		ZEROS DO( set(buf, rAdd(1), size))
		ZEROS DO( set(buf, NoiseWhite<>(), size))
	}
	
	#undef DO
*/
/*
	{	using namespace scratch;
	
		printf("\nDFT:\n");
		const uint32_t n=32;
		float in[n], out[n+2];
	
		tbl::cosine(in, n); //for(int i=0; i<n; i++) printf("% 6.3f\n", in[i]);
		dftR2C(out, in, n);
	
		for(int i=0; i<n+2; i+=2) printf("\t[%2d] % 6.3f % 6.3f\n", i>>1, out[i], out[i+1]);
	}
*/	
	return 0;
}








