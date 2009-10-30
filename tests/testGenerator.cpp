#include "Gamma.h"
#include "Access.h"
#include "Constants.h"
//#include "Visual.h"

//#define GEN(fnc) slice(arr,N) = fnc; printf("\n"#fnc":\n"); PRINT_TABLE
//
//#define PRINT_TABLE \
//	for(int i=0; i<N; i++){\
//		float v = arr[i];\
//		printf("[%2d] % 5.3f  ", i, v);\
//		printPlot(v, 32); printf("\n");\
//	}

using namespace gam;
using namespace gam::gen;

int main(int argc, char* argv[]){
	
	const int N = 8;
	double A[N];
	
	#define ASSERT(A, a,b,c,d,e,f,g,h) assert(A[0]==a && A[1]==b && A[2]==c && A[3]==d && A[4]==e && A[5]==f && A[6]==g && A[7]==h);
	
	slice(A,N) = Val<>(1);
	ASSERT(A, 1,1,1,1,1,1,1,1)
	
	slice(A,N) = Impulse<>();
	ASSERT(A, 1,0,0,0,0,0,0,0)

	slice(A,N) = RCos<>(1./N, 1);
	for(int i=0; i<N; ++i) assert(A[i] = cos(M_2PI*i/N));
	slice(A,N) = RCos<>(2./N, 0.5);
	for(int i=0; i<N; ++i) assert(A[i] = 0.5*cos(2*M_2PI*i/N));

	slice(A,N) = RSin<>(1./N, 1);
	for(int i=0; i<N; ++i) assert(A[i] = cos(M_2PI*i/N));
	slice(A,N) = RSin<>(2./N, 0.5);
	for(int i=0; i<N; ++i) assert(A[i] = 0.5*cos(2*M_2PI*i/N));

	slice(A,N) = RAdd<>(1,0);
	ASSERT(A, 0,1,2,3,4,5,6,7)
	slice(A,N) = RAdd<>(2,1);
	ASSERT(A, 1,3,5,7,9,11,13,15)

	slice(A,N) = RAdd1<>(0);
	ASSERT(A, 0,1,2,3,4,5,6,7)
	slice(A,N) = RAdd1<>(10);
	ASSERT(A, 10,11,12,13,14,15,16,17)

	slice(A,N) = RAddN<1>(0);
	ASSERT(A, 0,1,2,3,4,5,6,7)
	slice(A,N) = RAddN<2>(1);
	ASSERT(A, 1,3,5,7,9,11,13,15)

	slice(A,N) = RAddWrap<>(1,0,4);
	ASSERT(A, 0,1,2,3,0,1,2,3)
	slice(A,N) = RAddWrap<>(3,0,8);
	ASSERT(A, 0,3,6,1,4,7,2,5)
	slice(A,N) = RAddWrap<>(3,10,18,10);
	ASSERT(A, 10,13,16,11,14,17,12,15)
	slice(A,N) = RMul<>(2,1);
	ASSERT(A, 1,2,4,8,16,32,64,128)

	slice(A,N) = RMulAdd<>(2,1,0);
	ASSERT(A, 0,1,3,7,15,31,63,127)


//	{ using namespace gam::gen;
//
//		GEN(rAdd(1./N))
//		GEN(rAdd(-1./N, 1.))
//		GEN(val(1.))
//		GEN(rMul(pow(0.1, 1./N)))
//		GEN(rMul(pow(1./0.1, 1./N), 0.1))
//		GEN(rMulAdd(0.9, 0.1, -1.))
//		GEN(rMulAdd(0.9, -0.1, 1.))
//		GEN(nyquist(1))
//		GEN(recip(1.))
//		GEN(rCos<>(1./N))
//		GEN(rCos<>(1./N, 0.5))
//		GEN(Sin<>(M_2PI/N, 0))
//		GEN(Sin<>(M_2PI/N, M_PI_2, 0.5))
//		GEN(rSin(1./N))
//		GEN(rSin(1./N, 0.25))
//		GEN(rSin(1./N, 0., 0.5))
////		GEN(RSin2<>(1./N))
////		GEN(RSin2<>(1./N, 0.25))
////		GEN(RSin2<>(1./N, 0.25, 0.9))
//		
//		
////		{	printf("\n");
////			//int x3=-3, x2=-2, x1=-1; // 0,1,2,3, ...
////			//int x3=-2, x2=-1, x1=-1; // 0,0,1,1,2,2, ...
////			int x3=-5, x2=-2, x1=-3; //  ...
////			for(int i=0; i<64; i++){
////				int x0 = x1+x2-x3; x3=x2; x2=x1; x1=x0;
////				printf("\t[%2d] % d\n", i, x0);
////			}
////		}
//			
//			
//		
////		for(int k=0; k<16; k++){
////		for(int j=0; j<16; j++){
////		float c2 = k/15.*6-3;
////		float c1 = j/15.*6-3;
////		printf("\nc1 = %f, c2 = %f, c1/c2 = %f:", c1, c2, c1/c2);
////		GEN(Test<>(cos(-3*ang), cos(-2*ang)*0.9, cos(-ang)*0.6, 1, c2, c1))
////		}}
//		
//		//GEN(RSin2<>(-0.125, 0.125, 1.5, -0.9))
//	
//		//GEN(Test<>(cos(-3*ang), cos(-2*ang)*0.9, cos(-ang)*0.6, 1, -2.7, 2.7))
//	
//	/*
//		//GEN(TSin< Sin<> >( Sin<>(M_2PI/N, 0), 0 ) )
//		//GEN(TSin< Add<> >( Add<>(0.01, M_2PI/N), 0 ) )
//	
//		//GEN( TMulOp< Sin<>, Const<> >( Sin<>(M_2PI/N, 0), Const<>(0.5) ) )
//		
//		mem::set(arr, TMulOp< Sin<>, Val<> >( Sin<>(M_2PI/N, 0), Val<>(0.5), Loop(N) )); printf("\n"":\n"); PRINT_TABLE
//
//
//		Sin<> gen1(M_2PI/N, 0);
//		Val<> gen2(0.5);
//		mem::set(arr, mul(gen1, gen2), Loop(N)); printf("\n"":\n"); PRINT_TABLE
//
//		printf("\n");
////		mem::zero(arr, N);
////		for(int i=1; i < N>>1; ++i){
////			arr::add(arr, N, Sinusoid<>(M_2PI*i/N, M_PI_2));
////		}
////		PRINT_TABLE
//	*/
//	
//		double maxError = 0.00001;
////		{	bool result=false;
////			RSin<> iRSin(1./N);
////			for(int i=0; i<N; ++i){
////				float v = iRSin();
////				float a = sin((M_2PI/N)*i);
////				float e = scl::abs(v-a);
////				if(e>maxError) goto done;
////			}
////			result=true;
////			done: printf("RSin: %s\n", result ? "pass" : "fail");
////		}
//
//		#define TEST(Class, arg, truth)\
//		{	bool result=true;\
//			Class<> i##Class(arg);\
//			for(int i=0; i<N; ++i){\
//				float v = i##Class();\
//				float a = truth;\
//				float e = scl::abs(v-a);\
//				if(e>maxError){ result=false; break; };\
//			}\
//			printf("%8s: %s\n", #Class, result ? "pass" : "fail");\
//		}
//		
//		printf("\nGenerator Unit Tests\n");
//		TEST(Val,  1, 1)
//		TEST(Impulse, 1, i==0)
//		TEST(Nyquist, 1, scl::even(i)?1:-1)
//		TEST(Recip, 1, 1./(i+1))
//		TEST(RAdd, 1, i)
//		TEST(RAdd1, 0, i)
//		TEST(RMul, 0.9, pow(0.9,i))
//		TEST(RCos, 1./N, cos((M_2PI/N)*i))
//		TEST(RSin, 1./N, sin((M_2PI/N)*i))
//		TEST(Saw,  1./N, float(i)/N)
//	}
	
	
	
	
	return 0;
}

