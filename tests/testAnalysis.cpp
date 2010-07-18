#include <stdio.h>
#include <string.h>

#include "Gamma/arr.h"
#include "Gamma/gen.h"
#include "Gamma/tbl.h"
#include "Gamma/Access.h"
#include "Gamma/Print.h"

using namespace gam;
using namespace gam::gen;

int main(int argc, char* argv[]){

	const uint32_t size = 32;
	float table[size];
	//uint32_t indices[size];

	tbl::sinusoid(table, size, 0.f, 2.f);
	slice(table, size) += val(1);
	//arr::lineSlope1(table, size); arr::add(table, size, -8.f);
	//arr::mul(table, size, 0.9999f);
	slice(table, size) = rAdd(1./size, 0.);
	
	// print out function
	for(uint32_t i=0; i<size; i++){
		printf("[%4d]\t% 7.4f ", i, table[i]);
		printPlot(table[i]);
		printf("\n");
	}


	uint32_t featureI;
	float featureF;
	
	printf("\n");

	featureF = arr::dot(table, table, size);
	printf("Dot Product:       %f\n", featureF);

	uint32_t min, max;
	arr::extrema(table, size, min, max);
	printf("Extrema:           [%f, %f]\n", table[min], table[max]);

	float slope, inter;
	arr::fitLine(table, size, slope, inter);
	printf("Linear fit:        inter = %f, slope = %f\n", inter, slope);

	featureI = arr::max(table, size);
	printf("Max:               [%2d] %f\n", featureI, table[featureI]);

	uint32_t peaks[size>>1];
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
	uint32_t numIndices = size;
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
