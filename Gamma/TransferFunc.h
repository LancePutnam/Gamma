#ifndef GAMMA_TRANSFER_FUNC_H_INC
#define GAMMA_TRANSFER_FUNC_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <math.h>
#include <vector>
#include <complex>

namespace gam{


/// Transfer function of an arbitrary difference equation
    
/// http://en.wikipedia.org/wiki/Transfer_function
///
class TransferFunc {
public:
	typedef std::complex<double> Complex;

	/// \param[in] gain		overall filter gain
	TransferFunc(double gain=1): mGain(gain){}

	/// Add feedforward sample delay
	TransferFunc& addX(double c, double d){ mX.push_back(DelayUnit(c,d)); return *this; }

	/// Add feedback sample delay
	TransferFunc& addY(double c, double d){ mY.push_back(DelayUnit(c,d)); return *this; }
	
	/// Clear sample delays
	TransferFunc& clear(){ mX.clear(); mY.clear(); return *this; }
	
	/// Set overall gain factor
	TransferFunc& gain(double v){ mGain=v; return *this; }
	
	/// Returns frequency response at unit frequency [-0.5, 0.5]
	Complex operator()(double f){
		Complex X(0,0), Y(1,0);
		f *= M_2PI;
		for(uint32_t i=0; i<mX.size(); ++i) X += mX[i].response(f);
		for(uint32_t i=0; i<mY.size(); ++i) Y -= mY[i].response(f);
		return X/Y * mGain; // H(z) = Y(z)/X(z)
	}

protected:

	struct DelayUnit{
		/// param[in] c		weighting coefficient
		/// param[in] d		delay in samples
		DelayUnit(double c_, double d_): c(c_), d(d_){}

		Complex response(double f){
			double phs = f*d;
			return Complex(c*::cos(phs), c*::sin(phs));
		}

		double c, d;
	};

	std::vector<DelayUnit> mX, mY;
	double mGain;
};

} // gam::

#endif
