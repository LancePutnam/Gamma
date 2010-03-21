#include <assert.h>
#include <stdio.h>
#include <iostream>

#include "Gamma/arr.h"
#include "Gamma/gen.h"
#include "Gamma/mem.h"
#include "Gamma/Containers.h"

using namespace gam;
using namespace gam::gen;
using namespace std;

int main(int argc, char* argv[]){


	#define PRINT(a) printf("%p: ", &a); for(uint32_t i=0; i<a.size(); i++) printf("%d ", a[i]); printf("\n");
	{
		Array<int> array(8);
		PRINT(array)
		for(uint32_t i=0; i<array.size(); i++) array[i]=i;
		PRINT(array)
		
		Array<int> array2(array);
		PRINT(array2)
		array2.own();
		mem::zero(array.elems(), array.size());
		PRINT(array2)
		
		ArrayPow2<int> arrayPow2(8);
		PRINT(arrayPow2)
		for(uint32_t i=0; i<arrayPow2.size(); i++) arrayPow2[i]=i;
		PRINT(arrayPow2)		
	}
	#undef PRINT
	

//	double pos = 0;
//	for(unsigned j=1; j<34; j++){
//		printf("%2i ", (unsigned)pos);
//		pos = scl::wrap((pos + 4.25), 16.);
//	} printf("\n");

	int x = 1;
	DelayN<int> delayN(4);
	
	for(uint32_t i=0; i<delayN.size() * 2; i++){
		int input = x++;
		int d = delayN(input);
		printf("(%d) -> %d [", input, d);
		for(uint32_t j=0; j<delayN.size(); j++){
			printf("%d ", delayN[j]);
		}
		printf("]\n");
	}


	printf("\nDoubleRing:\n");
	DoubleRing<int> doubleRing(4);
	
	for(int i=0; i<8; ++i){
		doubleRing(i+1);
		
		printf("\t[");
		for(uint32_t j=0; j<doubleRing.size(); ++j) printf("%d ", doubleRing[j]);
		
		printf("] [");
		int * buf = doubleRing.copy();
		for(uint32_t j=0; j<doubleRing.size(); ++j) printf("%d ", buf[j]);
		
		printf("]\n");
	}

	
	printf("\nMulti:\n");
	{	Multi<5,char> multi = {{'m','u','l','t','i'}};
		for(uint32_t i=0; i<multi.size(); ++i) cout<<multi[i]; cout<<endl;}	
	{	Multi<5,char> multi = {{'z'}};
		for(uint32_t i=0; i<multi.size(); ++i) cout<<multi[i]; cout<<endl;}
	
	
	{ printf("\nRing:\n");
		const int N=4;
		Ring<int> ring(N);
		bool r=false;
		
		mem::zero(ring.elems(), ring.size());
		//r = mem::equal(ring, Val<int>(0), N);
		r = slice(&ring[0], N) == Val<int>(0);
		printf("\tzero(): %s\n", r ? "pass" : "fail");
		
		for(int i=0; i<N*2; ++i){
			ring(i+1);
			printf("\t"); for(int j=0; j<N; ++j) printf("%d ", ring[j]);
			printf("  "); for(int j=0; j<N; ++j) printf("%d ", ring.atPrev(j));
			printf("\n");
		}
	}
	
	{ printf("\nVec:\n");
		const int N=4;
		Vec<N, int> v1;
		Vec<N, int> v2;
		
		#define DO(op)\
			printf("\n%s\n", #op);\
			printf("\ti: v1:"); for(int i=0; i<N; ++i) printf("%2d", v1[i]);\
			printf("  v2:"); for(int i=0; i<N; ++i) printf("%2d", v2[i]); printf("\n");\
			op;\
			printf("\to: v1:"); for(int i=0; i<N; ++i) printf("%2d", v1[i]);\
			printf("  v2:"); for(int i=0; i<N; ++i) printf("%2d", v2[i]); printf("\n");
		
		DO(v1 = 3; v2 = 2)
		DO(v1 += v2)
		DO(v1 -= v2)
		DO(v1 *= v2)
		DO(v1 /= v2)
		DO(v1 += 1)
		DO(v1 -= 1)
		DO(v1 *= 2)
		DO(v1 /= 2)
		DO(v1 = v1 * 2 + v2)
		DO(v1 = -v2)
				
		#undef DO
	}
//	#define PRINT(obj, fnc) for(int i=0; i<obj.size(); ++i) cout << fnc; cout << endl;
//	
//	printf("\nSeq:\n");
//	
//	// forward
//	{	Seq<8,char> seq("12345678"); PRINT(seq, seq()) }
//
//	// backward
//	{	Seq<8,char, gen::RAddN<-1> > seq("12345678"); seq.tap() = 8; PRINT(seq, seq()) }
//
//	// stutter
//	{	Seq<8,char,gen::RAdd<float> > seq("12345678");
//		seq.tap() = -0.5; seq.tap().add = 0.5;
//		PRINT(seq, seq()) }
//
//	
//	{
//		printf("\nComplex\n");
//		Complex<double> c(Complex<double>::Polar(1, M_PI));
//		
//		for(int i=0; i<32; ++i){
//			for(int i=0; i<2; ++i) scl::printPlot(c[i], 16); printf("\n");
//			c *= c.normalize();
//		}
//	}
//
//
//	{
//		printf("\nQuat\n");
//		Quat<double> q(1,2/3.,1/3.,0), qr(13,1,1,1);
//		q.normalize();
//		qr.normalize();
//		
//		for(int i=0; i<64; ++i){
//			for(int i=0; i<4; ++i) scl::printPlot(q[i], 16); printf("\n");
//			q *= qr;
//		}
//	}
	
	return 0;
}

