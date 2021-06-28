#ifndef GAMMA_CONSTANTS_H_INC
#define GAMMA_CONSTANTS_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

	File Description: 
	Definitions of numerical constants.
*/


#include "Gamma/pstdint.h"

namespace gam{

#define DEFCONST(N, vf, vd)\
	template <class T> struct N;\
	template<> struct N< float>{ operator uint32_t() const { return vf##UL; } };\
	template<> struct N<double>{ operator uint64_t() const { return vd##ULL; } };

	/*		sign +	expo +	frac =	total
	float	1		 8		23		32
	double	1		11		52		64
	*/

	DEFCONST(MaskExpo, 0x7F800000, 0x7FF0000000000000) // IEEE-754 floating-point exponent bit mask
	DEFCONST(MaskFrac, 0x007FFFFF, 0x000FFFFFFFFFFFFF) // IEEE-754 floating-point fraction bit mask
	DEFCONST(MaskSign, 0x80000000, 0x8000000000000000) // IEEE-754 floating-point sign bit mask
	DEFCONST(Expo1_2 , 0x3F000000, 0x3FE0000000000000) // IEEE-754 floating-point [0.5, 1) exponent interval
	DEFCONST(Expo1   , 0x3F800000, 0x3FF0000000000000) // IEEE-754 floating-point [ 1, 2) exponent interval
	DEFCONST(Expo2   , 0x40000000, 0x4000000000000000) // IEEE-754 floating-point [ 2, 4) exponent interval
	DEFCONST(Expo4   , 0x40800000, 0x4010000000000000) // IEEE-754 floating-point [ 4, 8) exponent interval
	DEFCONST(ExpoNeg1, 0xBF800000, 0xBFF0000000000000) // IEEE-754 floating-point [-1,-2) exponent interval
#undef DEFCONST

const float justUnder1f = 0.99999997f; //0x3f7fffff;


/// 0 if little-endian, 1 if big-endian.
#if defined(__ppc__)
	const int endian = 1;
#else
	const int endian = 0;
#endif

inline bool isLittleEndian(){ return endian==0; }

namespace{
	const double roundMagic = 6755399441055744.; // 2^52 * 1.5
	
	template<class T> constexpr T roundEps();
	template<> inline constexpr float  roundEps<float >() { return 0.499999925f; }
	template<> inline constexpr double roundEps<double>() { return 0.499999985; }
}


// constant macros

// const double pi = 4*atan(1.0);

// math.h sometimes does not include these ANSI C macros.
#ifndef M_E
#define M_E			2.71828182845904523536028747135266250
#endif
#ifndef M_LOG2E
#define M_LOG2E		1.44269504088896340735992468100189214
#endif
#ifndef M_LOG10E
#define M_LOG10E	0.434294481903251827651128918916605082
#endif
#ifndef M_LN2
#define M_LN2		0.693147180559945309417232121458176568
#endif
#ifndef M_LN10
#define M_LN10		2.30258509299404568401799145468436421
#endif
#ifndef M_PI
#define M_PI		3.14159265358979323846264338327950288
#endif
#ifndef M_PI_2
#define M_PI_2		1.57079632679489661923132169163975144
#endif
#ifndef M_PI_4
#define M_PI_4		0.785398163397448309615660845819875721
#endif
#ifndef M_1_PI
#define M_1_PI		0.318309886183790671537767526745028724
#endif
#ifndef M_2_PI
#define M_2_PI		0.636619772367581343075535053490057448
#endif
#ifndef M_2_SQRTPI
#define M_2_SQRTPI	1.12837916709551257389615890312154517
#endif
#ifndef M_SQRT2
#define M_SQRT2		1.41421356237309504880168872420969808
#endif
#ifndef M_SQRT1_2
#define M_SQRT1_2	0.707106781186547524400844362104849039
#endif

// Some other useful constants
#ifndef M_2PI
#define M_2PI		6.283185307179586231941716828464095101		// 2pi
#endif
#ifndef M_4PI
#define M_4PI		12.566370614359172463937643765552465425		// 4pi
#endif
#ifndef M_1_2PI
#define M_1_2PI		0.159154943091895345554011992339482617		// 1/(2pi)
#endif
#ifndef M_3PI_2
#define M_3PI_2		4.712388980384689673996945202816277742		// 3pi/2
#endif
#ifndef M_3PI_4
#define M_3PI_4		2.356194490192343282632028017564707056		// 3pi/4
#endif
#ifndef M_LN001		
#define M_LN001		-6.90775527898								// ln(0.001)
#endif

// Function for getting values
// printf("%.36Lf\n", 4.0 * M_PI);

} // gam::

#endif
