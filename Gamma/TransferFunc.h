#ifndef GAMMA_TRANSFER_FUNC_H_INC
#define GAMMA_TRANSFER_FUNC_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <cmath>
#include <vector>
#include <complex>

namespace gam{

/// Transfer function of an arbitrary difference equation

/// http://en.wikipedia.org/wiki/Transfer_function
///
class TransferFunc {
public:
	typedef std::complex<double> Complex;

	struct DelayUnit{
		/// param[in] c		weighting coefficient
		/// param[in] d		delay in samples
		DelayUnit(double c_=0, double d_=0): c(c_), d(d_){}

		DelayUnit& weight(double v){ c=v; return *this; }

		DelayUnit& delay(double v){ d=v; return *this; }

		// H(z) = c z^d
		Complex response(Complex z){
			return c * pow(z,-d);
		}
		
		/// \param[in] f	frequency, in radians
		Complex response(double f){
			double phs = f*d;
			return Complex(c*cos(phs), c*sin(phs));
			// return response(Complex(cos(f), sin(f)));
		}

		double c, d;
	};


	/// \param[in] gain		overall filter gain
	TransferFunc(double gain=1): mGain(gain){}


	/// Set number of feedforward (x) and feedback (y) components
	TransferFunc& resize(int nx, int ny){
		mx.resize(nx);
		my.resize(ny);
		return *this;
	}

	DelayUnit& x(int idx){ return mx[idx]; }

	DelayUnit& y(int idx){ return my[idx]; }

	/// Add feedforward sample delay
	TransferFunc& addX(double c, double d){
		mx.push_back(DelayUnit(c,d)); return *this; }

	/// Add feedback sample delay
	TransferFunc& addY(double c, double d){
		my.push_back(DelayUnit(c,d)); return *this; }

	/// Clear sample delays
	TransferFunc& clear(){ mx.clear(); my.clear(); return *this; }

	/// Set overall gain factor
	TransferFunc& gain(double v){ mGain=v; return *this; }
	
	/// Returns frequency response at unit frequency [-0.5, 0.5]
	Complex operator()(double f){
		Complex X(0,0), Y(1,0);
		f *= M_2PI;
		for(unsigned i=0; i<mx.size(); ++i) X += mx[i].response(f);
		for(unsigned i=0; i<my.size(); ++i) Y -= my[i].response(f);
		return (X/Y) * mGain; // H(z) = Y(z)/X(z)
	}

	/// Returns frequency response at coordinate on z-plane
	Complex operator()(Complex z){
		Complex X(0,0), Y(1,0);
		for(unsigned i=0; i<mx.size(); ++i) X += mx[i].response(z);
		for(unsigned i=0; i<my.size(); ++i) Y -= my[i].response(z);
		return (X/Y) * mGain; // H(z) = Y(z)/X(z)
	}

protected:
	std::vector<DelayUnit> mx, my;
	double mGain;
};

} // gam::

#endif
