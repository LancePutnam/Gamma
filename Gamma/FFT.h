#ifndef GAMMA_FFT_H_INC
#define GAMMA_FFT_H_INC


namespace gam{


/// One-dimensional complex fast Fourier transform

/// The complex sequence format for forward and inverse transforms is
/// [r0, i0, r1, i1, ... , r(n-1), i(n-1)]
template <class T>
class CFFT{
public:

	typedef T value_type;

	/// \param[in] size		size of complex input sequence; 
	///						most efficient when a product of small primes
	CFFT(int size=0);
	
	~CFFT();

	/// Get size of transform
	int size() const;

	/// Perform forward transform in-place
	
	/// \param[in,out] buf		input/output buffer
	/// \param[in] normalize	whether to scale magnitudes by 1/N
	/// \param[in] nrmGain		gain to apply if normalizing
	void forward(T * buf, bool normalize=true, T nrmGain=1.);

	template <template <class> class ComplexType>
	void forward(ComplexType<T> * buf, bool normalize=true, T nrmGain=1.){
		forward((T*)buf, normalize, nrmGain);
	}

	/// Perform inverse transform in-place
	void inverse(T * buf);

	template <template <class> class ComplexType>
	void inverse(ComplexType<T> * buf){ inverse((T*)buf); }

	/// Set size of transform
	void resize(int n);

private:
	class Impl; Impl * mImpl;
};



/// One-dimensional real-to-complex fast Fourier transform

/// The complex sequence format for forward and inverse transforms is
/// [2r0, r1, i1, ... , r(n/2-1), i(n/2-1), 2r(n/2)]
template <class T>
class RFFT{
public:

	typedef T value_type;

	/// \param[in] size		size of real input sequence; 
	///						most efficient when a product of small primes
	RFFT(int size=0);
	
	~RFFT();

	/// Get size of transform
	int size() const;
	
	/// Perform real-to-complex forward transform in-place
	
	/// \param[in,out]	buf			input is real sequence, output is complex sequence
	/// \param[in]		complexBuf	If true, then 
	///									input is  [ *, x0, x1, x2, ..., x(n),   *] and
	///									output is [r0,  0, r1, i1, ..., r(n/2), 0].
	///								If false, then 
	///									input is  [x0, x1, x2, ..., x(n)  ] and
	///									output is [r0, r1, i1, ..., r(n/2)].
	/// \param[in]		normalize	whether to scale magnitudes by 1/N
	/// \param[in]		nrmGain		gain to apply if normalizing
	void forward(T * buf, bool complexBuf=false, bool normalize=true, T nrmGain=1.);
	
	/// Perform complex-to-real inverse transform in-place

	/// \param[in,out]	buf			input is complex sequence, output is real sequence
	/// \param[in]		complexBuf	If true, then 
	///									input is  [r0,  0, r1, i1, ..., r(n/2), 0] and
	///									output is [ *, x0, x1, x2, ..., x(n),   *].
	///								If false, then 
	///									input is  [r0, r1, i1, ..., r(n/2)]  and 
	///									output is [x0, x1, x2, ..., x(n)  ].
	void inverse(T * buf, bool complexBuf=false);

	/// Set size of transform
	void resize(int n);
	

private:
	class Impl; Impl * mImpl;
};

} // gam::

#endif
