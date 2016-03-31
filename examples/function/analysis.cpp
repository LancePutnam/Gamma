/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Function / Analysis
	Description:	Demonstration of various analysis functions.
*/

#include <stdio.h>
#include <string.h>
#include "Gamma/Gamma.h"
#include "Gamma/Access.h"
#include "Gamma/Print.h"

using namespace gam;
using namespace gam::gen;

int main(){

	const unsigned size = 32;
	float table[size];
	//unsigned indices[size];

	for(unsigned i=0; i<size; ++i){
		float phs = float(i)/size;
		table[i] = scl::abs(cos(2 * 2*M_PI*phs) * exp(-phs));
	}
	
	// print out function
	for(unsigned i=0; i<size; i++){
		printf("[%4d]\t% 7.4f ", i, table[i]);
		printPlot(table[i]);
		printf("\n");
	}


	unsigned featureI;
	float featureF;
	
	printf("\n");

	featureF = arr::dot(table, table, size);
	printf("Dot Product:       %f\n", featureF);

	unsigned min, max;
	arr::extrema(table, size, min, max);
	printf("Extrema:           [%f, %f]\n", table[min], table[max]);

	float slope, inter;
	arr::lineFit(table, size, slope, inter);
	printf("Linear fit:        inter = %f, slope = %f\n", inter, slope);

	featureI = arr::indexOfMax(table, size);
	printf("Max:               [%2d] %f\n", featureI, table[featureI]);

	unsigned peaks[size>>1];
	featureI = arr::maxima(peaks, table, size);
	printf("Maxima (%d):\n", featureI);
	for(unsigned long i=0; i<featureI; i++) printf("\t[%2d] %f\n", peaks[i], table[peaks[i]]);

	featureF = arr::mean(table, size);
	printf("Mean:              %f\n", featureF);

	featureF = arr::meanNorm(table, size);
	printf("Mean Norm:         %f\n", featureF);
	
	featureF = arr::meanWeightedIndex(table, size);
	printf("Mean Weight Index: %f\n", featureF);

	featureF = arr::rms(table, size);
	printf("RMS:               %f\n", featureF);

	featureF = arr::sum(table, size);
	printf("Sum:               %f\n", featureF);

	featureF = arr::variance(table, size);
	printf("Variance:          %f\n", featureF);

	featureI = arr::within(table, size, 0.5f);
	printf("Within 0.5:        %d\n", featureI);

	featureI = arr::zeroCount(table, size);
	printf("Zeros:             %d\n", featureI);
	
	featureI = arr::zeroCross(table, size, table[size-1]);
	printf("Zero Crossings:    %d\n", featureI);
	
	featureI = arr::zeroCrossN(table, size, table[size-1]);
	printf("Zero Crossings(-): %d\n", featureI);
/*
	printf("\nCluster:\n");
	arr::lineSlope1(indices, size);
	unsigned numIndices = size;
	arr::cluster(table, indices, numIndices, 0.1f);
	mem::print(table, indices, numIndices, "% f");

	printf("\nHistogram [0, 2):\n");
	mem::zero(indices, size);
	float step = 2.f / (float)size;
	arr::histogram(table, size, indices, size, 1.f / step, -1.f / step);
	for(unsigned long i=0; i<size; i++)
		printf("[%4.2f, %4.2f) %d\n", (float)i * step, (float)(i+1) * step, indices[i]);
	
	printf("\nSort - Insertion:\n");
	arr::lineSlope1(indices, size);
	arr::sortInsertion(table, indices, size);
	mem::print(table, indices, size, "% f");

	printf("\nSort - Quick:\n");
	arr::lineSlope1(indices, size);
	arr::sortQuick(table, indices, 0, size);
	mem::print(table, indices, size, "% f");
*/	
	return 0;
}
