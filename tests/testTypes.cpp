#include <stdio.h>
#include "Types.h"

using namespace gam;

int main(int argc, char* argv[]){


	printf("\nMulti:\n");
	{	Multi<5,char> multi = {'m','u','l','t','i'};
		for(int i=0; i<multi.size(); ++i) printf("%c", multi[i]); printf("\n");}
	{	Multi<5,char> multi = {'z'};
		for(int i=0; i<multi.size(); ++i) printf("%c", multi[i]); printf("\n");}

	
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
	
	
	{
//		//Elem2<int> elem2 = {1,2};
//		Elem2<int> elem2(1,2);
//		printf("Elem2<int>: sizeof=%d B, size=%d\n", sizeof(elem2), elem2.size());
//		printf("x,y\t{%d,%d}\n", elem2.x, elem2.y);
//		printf("r,i\t{%d,%d}\n", elem2.r, elem2.i);
//		printf("[0],[1]\t{%d,%d}\n", elem2[0], elem2[1]);
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

