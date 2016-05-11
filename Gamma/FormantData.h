#ifndef GAMMA_FORMANTDATA_H_INC
#define GAMMA_FORMANTDATA_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

	File Description: 
	Frequency/amplitude data of formants.
*/

#include <cmath>

namespace gam{

/// Formant data of vowel sounds

/// Peterson, G.E., and H.L. Barney, "Control Methods Used in a Study of 
///	the Vowels," Journal of the Acoustical Society of America, vol. 24, 1952.
class Vowel{
public:

	/// Voice types
	enum Voice{
		MAN=0,		/**< */
		WOMAN=1,	/**< */
		CHILD=2,	/**< */
		NUM_VOICES
	};

	/// Phonemes
	enum Phoneme{
		HEED=0,		/**< */
		HID=1,		/**< */
		HEAD=2,		/**< */
		HAD=3,		/**< */
		HOD=4,		/**< */
		HAWED=5,	/**< */
		HOOD=6,		/**< */
		WHOD=7,		/**< */
		HUD=8,		/**< */
		HEARD=9,	/**< */
		NUM_PHONEMES
	};

	/// Get amplitude of filter
	static float amp(Voice v, Phoneme p, int i){
		return pow(10., dB(v,p,i)/20.);
	}

	/// Get amplitude, in decibels, of filter
	static float dB(Voice v, Phoneme p, int i){
		static const float voweldB[NUM_PHONEMES][3] = {
			{-4, -24, -28},
			{-3, -23, -27},
			{-2, -17, -24},
			{-1, -12, -22},
			{-1, - 5, -28},
			{ 0, - 7, -34},
			{-1, -12, -34},
			{-3, -19, -43},
			{-1, -10, -27},
			{-5, -15, -20}
		};
		return voweldB[p][i];
	}

	/// Get center frequency of filter
	static float freq(Voice v, Phoneme p, int i){
		static const float vowelHz[NUM_VOICES][NUM_PHONEMES][3] = {
			{	{270, 2290, 3010},
				{390, 1990, 2550},
				{530, 1840, 2480},
				{660, 1720, 2410},
				{730, 1090, 2440},
				{570,  840, 2410},
				{440, 1020, 2240},
				{300,  870, 2240},
				{640, 1190, 2390},
				{490, 1350, 1690}
			},
			{	{310, 2790, 3310},
				{430, 2480, 3070},
				{610, 2330, 2990},
				{860, 2050, 2850},
				{850, 1220, 2810},
				{590,  920, 2710},
				{470, 1160, 2680},
				{370,  950, 2670},
				{760, 1400, 2780},
				{500, 1640, 1960}
			},
			{	{370, 3200, 3730},
				{530, 2730, 3600},
				{690, 2610, 3570},
				{1010,2320, 3320},
				{1030,1370, 3170},
				{680, 1060, 3180},
				{560, 1410, 3310},
				{430, 1170, 3260},
				{850, 1590, 3360},
				{560, 1820, 2160}
			},
		};
		return vowelHz[v][p][i];
	}

	/// Get fundamental frequency of voice
	static float pitch(Voice v, Phoneme p, int i){
		static const float vowelFund[NUM_VOICES][NUM_PHONEMES] = {
			{136, 135, 130, 127, 124, 129, 137, 141, 130, 133},
			{235, 232, 223, 210, 212, 216, 232, 231, 221, 218},
			{272, 269, 260, 251, 256, 263, 276, 274, 261, 261},
		};
		return vowelFund[v][p];
	}

};

} // gam::

#endif
