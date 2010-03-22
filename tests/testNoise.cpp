#include "Gamma/Noise.h"
#include "Gamma/scl.h"
#include "Gamma/Print.h"

using namespace gam;

int main(int argc, char* argv[]){

	const int num=128;

	printf("\nWhite / Pink / Brown:\n");
	NoiseWhite<> white;
	NoisePink<> pink;
	NoiseBrown<> brown;
	
	for(int i=0; i<num; i++){
		float vw = white();
		float vp = pink();
		float vb = brown();
		
		printf("% 6.4f ", vw);
		printPlot(vw, 20, true, "o");
		
		printf("% 6.4f ", vp);
		printPlot(vp, 20, true, "o");
		
		printf("% 6.4f ", vb);
		printPlot(vb, 20, true, "o");
		
		printf("\n");
	}
	
//	printf("\nGeneric noise:\n");
//	Noise<> noise;
//	
//	for(int i=0; i<num; i++){
//		float v = noise(RndUni<>());
//		printf("% 6.4f ", v);
//		scl::printPlot(v, 20, "o", true);		
//	}
	

	printf("\nRandom seeding test:\n");
	NoiseWhite<> nw1, nw2;
	printf("White: v1:%d v2:%d %s\n", nw1.rng.val, nw2.rng.val, nw1.rng.val != nw2.rng.val ? "success" : "fail");

	NoisePink<> np1, np2;
	printf("Pink : v1:%d v2:%d %s\n", np1.rng.val, np2.rng.val, np1.rng.val != np2.rng.val ? "success" : "fail");
	
	NoiseBrown<> nb1, nb2;
	printf("Brown: v1:%d v2:%d %s\n", nb1.rng.val, nb2.rng.val, nb1.rng.val != nb2.rng.val ? "success" : "fail");
}
