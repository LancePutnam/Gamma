#ifndef GAMMA_RAMPED_H_INC
#define GAMMA_RAMPED_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

	File Description: 
	Value smoother
*/

namespace gam{

/// Applies linear smoothing to value updates
template <class T>
class Ramped {
public:

	Ramped(){}
	Ramped(T v){ *this = v; }

	/// Set new target value

	/// @param[in] v		new target value
	/// @param[in] steps	number of samples to ramp over
	Ramped& target(const T& v, float steps){
		mTarget = v;
		mInc = (mTarget - mVal)/steps;
		return *this;
	}

	/// Get target value
	const T& target() const { return mTarget; }

	/// Set current and target values
	Ramped& operator= (const T& v){ mVal=mTarget=v; return *this; }

	/// Get smoothed value
	operator T() const { return mVal; }

	/// Update value; \returns smoothed value
	T operator()(){
		mVal += mInc;
		if((mInc>T(0) && mVal>mTarget) || (mInc<T(0) && mVal<mTarget))
			mVal = mTarget;
		return mVal;
	}

private:
	T mVal = T(0);
	T mInc = T(0);
	T mTarget = T(0);
};

} // gam::

#endif // include guard
