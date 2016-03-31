/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Algorithmic / Note-based functions
	Description:	This demonstrates various functions for mapping notes onto
					numerical values and absolute frequencies.
*/

#include <stdio.h>
#include "Gamma/scl.h"
#include "Gamma/rnd.h"
using namespace gam;

int main(){

	/* 
	scl::freq returns a frequency, in Hz, from a note string consisting of a 
	pitch class letter followed by an (optional) accidental and lastly an 
	octave number where 4 is the middle octave.
	*/

	// Middle C
	printf("\"c 4\" = %g Hz\n", scl::freq("c 4"));
	printf("\"C4\"  = %g Hz\n", scl::freq("C4"));

	// Middle C sharp
	printf("\"c+4\" = %g Hz\n", scl::freq("c+4"));
	printf("\"C#4\" = %g Hz\n", scl::freq("C#4"));
	
	// Middle C flat
	printf("\"c-4\" = %g Hz\n", scl::freq("c-4"));
	printf("\"Cb4\" = %g Hz\n", scl::freq("Cb4"));
	
	// Low A
	printf("\"a3\"  = %g Hz\n", scl::freq("a3"));
	
	// The C major scale
	printf("\nFrequencies of the C-major scale (middle octave):\n");
	printf("%g, %g, %g, %g, %g, %g, %g, %g\n",
		scl::freq("c4"), scl::freq("d4"), scl::freq("e4"), scl::freq("f4"),
		scl::freq("g4"), scl::freq("a4"), scl::freq("b4"), scl::freq("c5")
	);
	
	
	/* 
	We can use scl::nearest to get the closest note number within a specified
	pitch class set. This is useful for mapping arbitrary values onto notes in a
	musical scale.
	The first argument is the note number to match, the second is a string
	containing the scale's half-step intervals, and the third argument specifies
	the modulus or divisions per octave (default is 12).
	*/
	
	// Major scale
	printf("%g\n", scl::nearest(4.1, "2212221"));
	
	// Minor scale
	printf("%g\n", scl::nearest(4.1, "2122212"));
	
	// Major pentatonic scale
	printf("%g\n", scl::nearest( 1.1, "22323"));
	printf("%g\n", scl::nearest( 4.1, "22323"));
	printf("%g\n", scl::nearest(10.7, "22323"));
	
	
	/*
	scl::ratioET computes frequency ratios within a specified 
	equal-temperament tuning system. The first argument is the note number, the
	second is the number of divisions per octave (default is 12) and the third
	is the base multiplier of the octave (default is 2).
	*/
	
	// Ratio of 1 semitone in 12-TET
	printf("1 semitone in 12-TET: %g\n", scl::ratioET(1));

	// Ratio of 7 semitones (a fifth) in 12-TET
	printf("7 semitone in 12-TET: %g\n", scl::ratioET(7));

	// Ratio of 117 cents in 12-TET
	printf("117 cents in 12-TET:  %g\n", scl::ratioET(1.17));
	
	// Ratio of 1 semitone in 7-TET
	printf("1 semitone in 7-TET:  %g\n", scl::ratioET(1, 7));

	// Ratio of 1 semitone in Bohlen-Pierce scale
	printf("1 semitone in Bohlen-Pierce scale: %g\n", scl::ratioET(1, 13, 3));


	/*
	The functions above can be combined in order to get specific frequency
	values of notes in a scale.
	*/

	printf("\nRandom frequency values in C-major pentatonic scale:\n");
	for(int i=0; i<10; ++i){
		// First, produce a random note number
		double v = rnd::uni(12.);
		
		// Next, get the nearest note number in the pentatonic scale
		v = scl::nearest(v, "22323");
		
		// Finally, compute a frequency in Hz from the note number
		v = scl::ratioET(v) * scl::freq("c4");

		// Everything in one line...
		//double v = scl::ratioET(scl::nearest(rnd::uni(12.), "22323")) * scl::freq("c4");

		printf("%g Hz\n", v);
	}
}
