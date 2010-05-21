#ifndef GAMMA_FILTER_INC_H
#define GAMMA_FILTER_INC_H

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Gamma/Sync.h"
#include "Gamma/Types.h"

namespace gam{


/// First-order IIR section (1 pole/ 1 zero)

/// y[n] = ci0*x[n] + ci1*x[n-1] + co1*y[n-1]
///
template <class Tv, class Tp>
class IIR1{
public:

	IIR1(): d1(0), ci0(0), ci1(0), co1(0){}

	/// Filter input sample
	Tv operator()(const Tv& i0){
		Tv d0 = i0 + d1*co1;
		Tv o0 = d0*ci0 + d1*ci1;
		d1=d0;
		return o0;
	}

	/// Set filter coefficients
	void coefs(Tp vi0, Tp vi1, Tp vo1){ ci0=vi0; ci1=vi1; co1=vo1; }
	
	/// Zero internal delays
	void zero(){ d1=0; }

protected:
	Tv d1;
	Tp ci0, ci1, co1;
};



/// Second-order IIR section (2 poles, 2 zeros)

/// y[n] = ci0*x[n] + ci1*x[n-1] + ci2*x[n-2] + co1*y[n-1] + co2*y[n-2]
///
template <class Tv=gam::real, class Tp=gam::real>
class IIR2{
public:

	IIR2(): d1(0), d2(0), ci0(0), ci1(0), ci2(0), co1(0), co2(0){}

	/// Filter input sample
	Tv operator()(const Tv& i0){	// direct form II
		Tv d0 = i0 + d1*co1 + d2*co2;
		Tv o0 = d0*ci0 + d1*ci1 + d2*ci2;
		d2=d1; d1=d0;
		return o0;
	}

	/// Get filter coefficients
	Tp * coefs(){ return c; }

	/// Set filter coefficients
	void coefs(Tp vi0, Tp vi1, Tp vi2, Tp vo1, Tp vo2){
		ci0=vi0; ci1=vi1; ci2=vi2; co1=vo1; co2=vo2;
	}
	
	/// Zero internal delays
	void zero(){ d1=0; d2=0; }

protected:
	Tv d1, d2;
	union{ struct{ Tp ci0, ci1, ci2, co1, co2; }; Tp c[5]; };
};



template <int N, class Tv=gam::real, class Tp=gam::real>
class IIRN{
public:

	Tv d[N];	///< delay buffer
	Tp a[N];	///< feedforward coefficients
	Tp b[N];	///< feedback coefficients

	/// Filter input sample
	Tv operator()(const Tv& i0){	// direct form II
		d[0] = i0;
		Tv o0 = 0;
		for(int i=N-1; i>=1; --i){
			d[0] += d[i]*b[i];	// sum feedback
			o0   += d[i]*a[i];	// sum feedforward
			d[i] = d[i-1];		// shift delay
		}
		o0 += d[0]*a[0];

		return o0*b[0];
	}
	
	/// Zero internal delays
	void zero(){ for(int i=0; i<N; ++i) d[i]=Tv(0); }
	
};

/// General-purpose multi-staged IIR filter

// The process for designing an IIR is:
//		1. pre-warp center frequency
//		2. determine location of poles and zeros (roots) on s-plane
//		3. apply bilinear transform on roots to go from s-plane to z-plane
//		4. convert z-plane roots to filter coefficients
//
// s-plane characteristics
//	+ im	+ frequencies
//	- im	- frequencies
//	- re	growing amplitudes
//	+ re	decaying amplitudes

template <class Tv=gam::real, class Tp=gam::real>
class IIRSeries{
public:

	/// @param[in]	ord		order of filter
	IIRSeries(uint32_t ord)
	:	mOrder(ord), mPoles(0), mIIR1(0), mIIR2(0)
	{
		mPoles = new complex[numPoles()];
		if(odd()) mIIR1 = new IIR1<Tv, Tp>;
		mIIR2 = new IIR2<Tv, Tp>[numIIR2()];

//		mPoleAngle.fromPhase(M_PI/order());
//		mPoleStart.fromPhase(-M_PI_2 + M_PI_2/order());

		// precompute s-plane pole angles in (-,+) quarter
		complex p0, p1;							// position, angle
		p0.fromPhase(-M_PI_2 + M_PI_2/order());	// start near (0, 1)
		p1.fromPhase(M_PI/order());				// go clockwise towards (-1, 0)

		for(uint32_t i=0; i<numPoles(); ++i){
			mPoles[i] = p0; p0 *= p1;
		}
	}

	virtual ~IIRSeries(){
		delete [] mPoles; mPoles=0;
		if(mIIR1) delete mIIR1;
		delete [] mIIR2; mIIR2=0;
	}
	
	/// Filter input sample
	Tv operator()(const Tv& i0) {
		Tv o0 = i0;
		for(uint32_t i=0; i<numIIR2(); ++i) o0 = mIIR2[i](o0);
		if(odd()) o0 = (*mIIR1)(o0);
		return o0;
	}

protected:
	typedef Complex<double> complex;

	uint32_t mOrder;		// order of the filter (# pole/zero pairs)
	complex * mPoles;		// array of non-conjugated poles (others mirrored around imaginary axis)

	IIR1<Tv, Tp> * mIIR1;	// pointer to first order stage for odd order filters
	IIR2<Tv, Tp> * mIIR2;	// array of second order stages
	
	uint32_t numPoles() const { return (order()+1)>>1; }	// # non-conjugated poles
	uint32_t numIIR2() const { return order()>>1; }			// # 2nd-order IIRs
	uint32_t odd() const { return order()&1; }				// is filter odd order?
	uint32_t order() const { return mOrder; }				// filter order
	complex& pole1(){ return mPoles[(order()-1)>>1]; }

	// Convert z-plane pole to first-order IIR coefs
	void poleToCoef1(){
		if(odd()){
			double t = pole1().r;
			t = (t - 1.)/(t + 1.);
			mIIR1->coefs(1, 1, -t);
		}
	}
	
	// Conformal map from s-plane to z-plane.
	// The transform is p = (1+p)/(1-p), a special type of Mobius transform.
	template <class T>
	static void bilinear(Complex<T>& p){
		Complex<T> den(T(1) - p.r, -p.i);
		p.r += Tp(1);
		p /= den;
	}

	// Convert a z-plane pole real component to IIR1 low-pass coefs
	static void convertLP(IIR1<Tv,Tp>& f, Tp pr){
		pr = (pr - Tp(1))/(pr + Tp(1));			
		Tp ci = (Tp(1) + pr) * Tp(0.5);
		f.coefs(ci, ci, -pr);
	}

	// Convert a z-plane pole to IIR2 low-pass coefs
	static void convertLP(IIR2<Tv,Tp>& f, const Complex<Tp>& p){	
		Tp co2 = Tp(1)/p.dot();						// recip of mag squared
		Tp co1 = Tp(-2)*p.r*co2;
		Tp ci  = (Tp(1) + co1 + co2) * Tp(0.25);	// input gain compensation
		f.coefs(ci, ci*Tp(2), ci, -co1, -co2);
	}
};


#define INHERIT_IIRSERIES\
	typedef IIRSeries<Tv,Tp> super;\
	typedef Complex<double> complex;\
	using super::mIIR1;\
	using super::mIIR2;\
	using super::mPoles;\
	using super::bilinear;\
	using super::numIIR2;\
	using super::numPoles;\
	using super::order;\
	using super::odd;\
	using super::convertLP;


/// Butterworth filter

/// Maximally flat in the pass- and stop-bands.
///
template <class Tv=gam::real, class Tp=gam::real>
class IIRButter: public IIRSeries<Tv,Tp>{
public:

	IIRButter(Tp frq, uint32_t order=2): super(order){
		freq(frq);
	}

	/// Set cutoff frequency (0, 0.5)
	void freq(Tp v){

		v = tan(M_PI*v);	// pre-warp frequency

		// convert IIR2s
		uint32_t j=0;
		for(; j<numIIR2(); ++j){

			// get s-plane pole
			Complex<Tp> p(Tp(mPoles[j].r) * v, Tp(mPoles[j].i) * v);

			bilinear(p);				// bilinear xform
			convertLP(mIIR2[j], p);		// z-plane pole to coefs
		}
		
		// convert IIR1
		if(odd()) convertLP(*mIIR1, Tp(mPoles[j].r) * v);
	}

protected:
	INHERIT_IIRSERIES;
};




/// Chebyshev filter

/// A Chebyshev filter is a special case of a Butterworth filter with a narrower
/// transition band but ripple in the passband.

// TODO: odd orders are unity gain, but even orders are not...
//
template <class Tv=gam::real, class Tp=gam::real>
class IIRCheby: public IIRSeries<Tv,Tp>{
public:

	IIRCheby(Tp frq, Tp rip, uint32_t order=2)
	:	super(order), mWarpR(1.), mWarpI(1.){
		set(frq, rip);
	}

	/// Set cutoff frequency (0, 0.5)
	void freq(Tp v){

		v = tan(M_PI*v);		// pre-warp frequency
		
		Tp mr = mWarpR * v;		// bring poles closer to jw axis
		Tp mi = mWarpI * v;	

		// convert IIR2s
		uint32_t j=0;
		for(; j<numIIR2(); ++j){

			// get s-plane pole
			Complex<Tp> p(Tp(mPoles[j].r) * mr, Tp(mPoles[j].i) * mi);

			bilinear(p);				// bilinear xform
			convertLP(mIIR2[j], p);		// z-plane pole to coefs
		}
		
		// convert IIR1
		if(odd()) convertLP(*mIIR1, Tp(mPoles[j].r) * mr);
	}

	/// Set cutoff frequency (0, 0.5) and passband ripple (dB)
	void set(Tp frq, Tp rip){
		ripple(rip); freq(frq);
	}

protected:
	INHERIT_IIRSERIES;
	
	Tp mWarpR, mWarpI;	// multiplicands to warp circle to ellipse

	void ripple(Tp v){
		v = pow(10., v/10.) - 1.;
		v = 1./pow(v, 1./order());
		Tp asinh = log(v + sqrt(1. + v*v));
		Tp v0 = asinh/(double)order();
	
		mWarpR = sinh(v0);
		mWarpI = cosh(v0);
	}
};



//template <class Tv=gam::real, class Tp=gam::real>
//class IIRCheby: public IIRSeries<Tv,Tp>{
//public:
//
//	IIRCheby(Tp frq, Tp rip, uint32_t order=2): super(order){
//		set(frq, rip);
//	}
//
//
//	// frq - desired cutoff frequency (0, 0.5)
//	// rip = passband ripple in dB
//	void set(Tp frq, Tp rip){
//
//		rip = pow(10., rip/10.) - 1.;
//		rip = pow(rip, 1./order());
//
//		// pre-warped angular frequency
//		frq = tan(M_PI*frq);
//
//		getPoles(frq, rip);
//		bilinear();
//		
//		//mGain = 1;
//		poleToCoef1();
//		poleToCoef2();
//	}
//	
//
//protected:
//
//	INHERIT_IIRSERIES;
//
//
//	// The poles form an ellipse in the complex plane.
//	void getPoles(double w, double rip) {
//
//		double x = 1./rip;
//		double asinh = log(x + sqrt(1. + x*x));
//		double v0 = asinh/(double)order();
//		double sm = sinh(v0);
//		double cm = cosh(v0);
//
//		complex& cr = mPoleAngle;
//		complex cv  = mPoleStart * w;
//		complex * p = &pole1();
//
//		(*p)(-sm*cv.r, cm*cv.i);
//		for(uint32_t j=1; j<numPoles(); ++j){
//			cv *= cr;
//			(*--p)(-sm*cv.r, cm*cv.i);
//		}
//	}
//
//	// Convert roots to 2nd-order IIR coefs
//	void poleToCoef2(){
//		for(int j=0; j<numIIR2(); j++){
//			double rr = mPoles[j].r;		// root real comp.
//			double rm = mPoles[j].dot();	// root mag squared
//			mIIR2[j].coefs(1,2,1, 2.*rr, -rm);
//			//mGain *= (rm - 2.*rr + 1.) / (4.*rm); 
//		}
//	}
//	
//};


#undef INHERIT_IIRSERIES

/*

//! Do bilinear transformation
// n2 = half the order
// r = complex roots
void bilinear(long n2, complex<double>* r){
	double td;
	int j;

	const double a = 2.; // butter, ellip
	const double a = 1.; // cheby

	if(odd){
		// For s-domain zero (or pole) at infinity we
		// get z-domain zero (or pole) at z = -1
		double tmp = (r[0].real() - a)/(a + r[0].real());
		iir_1->set_coeff(-tmp);
		gain *= 0.5*(1+tmp);
	}
	
	// cheby, butter
	for(j=odd;j<n2;j++){                       
		td = a*a - 2*a*r[j].real() + magsq(r[j]);
		r[j] = complex<double>(
			a*a - magsq(r[j]),
			2.0*a*r[j].imag()
		)/td;
	}

	// ellip
	for(j=odd; j<n2; ++j){
		r[j] = (1. + r[j])/(1. - r[j]);
	}
}


void get_coeff(long n2){
	for(int j=odd; j<n2; j++){
	
		double rr = roots[j].r;		// root real comp.
		double rm = roots[j].dot(); // root mag squared
		double zr = zeros[j].r;		// zero real comp.
		double zm = zeros[j].dot(); // zero mag squared
	
		// cheby
		iir[j-odd].set_coeff(-2*rr, rm);
		
		// butter
		iir[j-odd].set_coeff(-2*rr/rm, 1./rm);

		// ellip
		iir[j-odd].set_a(-2*rr, rm);
		iir[j-odd].set_b(-2*zr, zm);
		
		// all
		gain *= (rm - 2*rr + 1.) / (4*rm); 
	}
}
*/




// Rotation amount

// This object is used to construct complex rotation filters. If the bandwidth
// is 0, it effectively acts as a quadrature oscillator. Non-zero bandwidth
// amounts can be used for decaying complex sinusoids or complex two-pole 
// band-pass filters.
template <class Tv=double, class Ts=Synced>
class Rotation : public Complex<Tv>, Ts{
public:

	Rotation(Tv freq=440, Tv width=100)
	:	mFreq(freq), mWidth(width)
	{	Ts::initSynced(); }
	
	/// Set 60 dB decay length
	void decay(const Tv& v){
		//width(-M_LN001/(v * M_PI));
		width(Tv(2.198806796637603 /* -ln(0.001)/pi */)/v);
	}
	
	void width(const Tv& v){
		mWidth = v;
		mDecay = ::exp(-M_PI * v * Ts::ups());
		freq(mFreq);
	}

	void freq(const Tv& v){
		mFreq = v;
		Tv phaseInc = v * Ts::ups() * (Tv)M_2PI;
		this->fromPolar(mDecay, phaseInc);
	}
	
	Tv gain() const {
		return ((Tv)1 - mDecay*mDecay) /** (*this).i*/;
	}
	
	virtual void onResync(double r){ width(mWidth); }

protected:
	Tv mFreq, mWidth, mDecay;
};



/*
 -------------------------------------- A collection of IIR filters ----

        R in the filters means resonance, steepness and narrowness.

        *** Halfband lowpass ***

          b1  =  0.641339     b10 = -0.0227141
          b2  = -3.02936
          b3  =  1.65298      a0a12 = 0.008097
          b4  = -3.4186       a1a11 = 0.048141
          b5  =  1.50021      a2a10 = 0.159244
          b6  = -1.73656      a3a9  = 0.365604
          b7  =  0.554138     a4a8  = 0.636780
          b8  = -0.371742     a5a7  = 0.876793
          b9  =  0.0671787    a6    = 0.973529

          output(t) = a0a12*(input(t  ) + input(t-12))
                    + a1a11*(input(t-1) + input(t-11))
                    + a2a10*(input(t-2) + input(t-10))
                    + a3a9* (input(t-3) + input(t-9 ))
                    + a4a8* (input(t-4) + input(t-8 ))
                    + a5a7* (input(t-5) + input(t-7 )) + a6*input(t-6)
                    + b1*output(t-1) + b2*output(t-2) + b3*output(t-3)
                    + b4*output(t-4) + b5*output(t-5) + b6*output(t-6)
                    + b7*output(t-7) + b8*output(t-8) + b9*output(t-9)
                    + b10*output(t-10)

        *** Fast bandpass ***

          freq = passband center frequency
          r = 0..1, but not 1

          fx = cos(2*pi*freq/samplerate)
          a0 = (1-r)*sqrt(r*(r-4*fx^2+2)+1)
          b1 = 2 * fx * r
          b2 = -r^2

          output(t) = a0*input(t)
                    + b1*output(t-1) + b2*output(t-2)



        -------------------------------------- A collection of FIR filters ----

        *** 90 degree phase shift ***

          Coefficients from: (1-cos(t*pi))/(t*pi)
          Limit t->0 = 0

          Coefficients at even t values are zero.
          Antisymmetrical around t=0.

        *** Phase shift of any angle ***

          a = angle

          Coefficients from: (sin(t*pi-a)+sin(a))/(t*pi)
          Coefficients at even t values, not including t=0, are zero.
          Limit t->0 = cos(a)

        *** Lowpass ***

          f = cutoff frequency / (samplerate/2)

          Coefficients from: sin(t*f*pi)/(t*pi)
          Limit t->0 = f

          Symmetrical around t=0.

        *** Halfband lowpass ***

          Coefficients from: 0.5*sin(t*pi/2)/(t*pi/2)
          Limit t->0 = 0.5

          Coefficients at even t values, not including t=0, are zero.
          Symmetrical around t=0.

        *** Highpass ***

          f = cutoff frequency / (samplerate/2)

          Coefficients from: (sin(t*pi)-sin(t*f*pi))/(t*pi)
          Limit t->0 = 1-f

          Symmetrical around t=0.

        *** Bandpass ***

          f1 = lower frequency limit / (samplerate/2)
          f2 = higher frequency limit / (samplerate/2)

          Coefficients from: sin(t*f2*pi)/(t*pi) - sin(t*f1*pi)/(t*pi)
          Limit t->0 = f2-f1

          Symmetrical around t=0.

        *** Bandpass magnitude-ramp ***

          x1 = lower limit frequency / (samplerate/2)
          x2 = higher limit frequency / (samplerate/2)
          y1 = magnitude at lower limit frequency
          y2 = magnitude at higher limit frequency

          c = (y1-y2)/(x1-x2)
          d = (x1*y2-x2*y1)/(x1-x2)

          Coefficients from:
          ((d+c*x2)*sin(x2*t*pi)+(-d-c*x1)*sin(x1*t*pi)+
          (c*cos(x2*t*pi)-c*cos(x1*t*pi))/(t*pi))/(t*pi)
          Limit t->0 = (x2-x1)*(y1+y2)/2

          Symmetrical around t=0.

        *** Negative frequency removal ***

          Complex coefficients from: sin(t*pi)/(t*pi) + i*(cos(t*pi)-1)/(t*pi)
          Limit t->0 = 1

          Real part symmetrical around t=0.
          Imaginary part antisymmetrical around t=0.

        *** Combined negative frequency removal and bandpass ***

          f1 = lower frequency limit / (samplerate/2)
          f2 = higher frequency limit / (samplerate/2)

          Complex coefficients from:
          (sin(f2*t*pi)-sin(f1*t*pi))/(t*pi) +
          i*(cos(f2*t*pi)-cos(f1*t*pi))/(t*pi)
          Limit t->0 = f2-f1

          Real part symmetrical around t=0.
          Imaginary part antisymmetrical around t=0.
*/






} // gam::
#endif
