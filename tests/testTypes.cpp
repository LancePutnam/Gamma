#include <stdio.h>
#include "Gamma/Types.h"
#include "Gamma/ipl.h"
#include "Gamma/scl.h"

using namespace gam;


template <uint32_t N, class T, class S>
struct TVec : public Multi<N,T> {

	//typedef Vec<N,T> V;

	TVec(const T& v=T(0)){ (*this) = v; }

	#define DO for(uint32_t i=0; i<N; ++i)
	S& operator +=(const S& v){ DO (*this)[i] += v[i]; return *this; }

	T dot() const { return dot(*this); }
	T dot(const S& v) const { T r=(T)0; DO r+=(*this)[i]*v[i]; return r; }
	T norm() const { return sqrt(dot()); }

	S sgn() const { S(*this) /= norm(); }

	#undef DO
};



template <class T>
struct TVec3 : public TVec<3, T, TVec3<T> >{

	//typedef TVec<3, T, TVec3<T> > base;
	typedef TVec<3, T, TVec3<T> > V;
	//using base::base();
	TVec3(const V& v){ (*this)=v; }
	TVec3(const T& v=T(0)){ (*this)(v,v,v); }
	TVec3(const T& v1, const T& v2, const T& v3=T(0)){ (*this)(v1,v2,v3); }
	
	TVec3& operator()(const T& v1, const T& v2, const T& v3){ TVec3& t=*this; t[0]=v1; t[1]=v2; t[2]=v3; return t; }
};


int main(int argc, char* argv[]){

	{
		Vec3<float> v(1,0,0), u(0,1,0);

		//v = ipl::linear(0.5, v, u*0.5);
	}


	printf("\nMulti:\n");
	{	Multi<5,char> multi = {{'m','u','l','t','i'}};
		for(uint32_t i=0; i<multi.size(); ++i) printf("%c", multi[i]); printf("\n");}
	{	Multi<5,char> multi = {{'z'}};
		for(uint32_t i=0; i<multi.size(); ++i) printf("%c", multi[i]); printf("\n");}

	
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
		//Elem2<int> elem2 = {1,2};
//		Elem2<int> elem2(1,2);
//		printf("Elem2<int>: sizeof=%d B, size=%d\n", sizeof(elem2), elem2.size());
//		printf("x,y\t{%d,%d}\n", elem2.x, elem2.y);
//		printf("r,i\t{%d,%d}\n", elem2.r, elem2.i);
//		printf("[0],[1]\t{%d,%d}\n", elem2[0], elem2[1]);

		
//		GComplex<float> vec2(1,2);
//		printf("Thing: sizeof=%d B, size=%d\n", sizeof(vec2), vec2.size());
//		vec2 += GComplex<float>(5,5);
//		vec2 += 1;
//		vec2 -= 1;
//		vec2 -= GComplex<float>(5,5);
//		vec2 *= GComplex<float>(2,2);
//		vec2 *= 3;
//		printf("x,y\t{%f,%f}\n", vec2.x, vec2.y);
//		printf("r,i\t{%f,%f}\n", vec2.r, vec2.i);
//		printf("[0],[1]\t{%f,%f}\n", vec2[0], vec2[1]);
//		//sss;sss
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
	{
		printf("\nComplex\n");
		Complex<double> c(Polar<double>(0.5, 0.1));
		
		printf("real:  %f imag: %f\n", c.r, c.i);
		printf("norm:  %f\n", c.norm());
		printf("norm2: %f\n", c.norm2());
		printf("arg:   %f\n", c.arg());
		printf("sgn:   %f, %f\n", c.sgn().r, c.sgn().i);
		printf("recip: %f, %f\n", c.recip().r, c.recip().i);
	}
	
//	{
//		//Vec3f v1, v2;
//		//Vec3f t = (v1-v2).cross(v2);
//		Vector3<float> sub;
//		sub(1,2,3);
//		sub += sub;
//		printf("%f\n", sub.doSomething().noop(4));
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

