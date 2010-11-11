#include <stdio.h>
#include "Gamma/Timer.h"

using namespace gam;

double error(double m, double g){ return (m - g) / g * 100.; }

int main(){

	Timer timer;

	printf("Current time: %lld nsec\n", timeNow());

	printf("\nSleep function test:\n");
	double wait = 0.0015;	// sleep interval (sec)
	int intervals = 4;		// number of interval doublings to test
	int trials = 100;		// number of timing measurements per interval

	for(int j=0; j<intervals; ++j){
		printf("\nDesired sleep time: % 5.3f sec\n", wait);
		printf("Trials: %d\n", trials);
		double min = wait * 2.;
		double avg = 0.;
		for(int i=0; i<trials; ++i){
			sleepSec(0.01);
			timer.start();
			sleepSec(wait);
			timer.stop();
			double dur = timer.elapsedSec();
			if(dur < min && dur > 0.) min = dur;
			avg += dur;
		}
		avg /= (double)trials;
		printf("Min: % 8.6f (% 5.2g%% error)\n", min, error(min, wait));
		printf("Avg: % 8.6f (% 5.2g%% error)\n", avg, error(avg, wait));
		wait *= 2.;
	}
	

	printf("\nPeriodic timing test:\n");
	wait = 1./60.;			// timing interval (sec)
	double dur = 10.;		// total duration of test (sec)
	
	nsec_t nwait = (nsec_t)(wait * 1e9);
	trials = (int)(dur / wait);
	printf("\nTimer period: %f sec\n", wait);
	printf("Total duration: %5.3f sec\n", dur);
	
	nsec_t target = timeNow() + nwait;
	timer.start();
	for(int i=0; i<trials; ++i){

		// do some work...
		if(i % (trials/10) == 0){
			printf("%d ", 10 - (i * 10) / trials); fflush(stdout);
		}

		sleepUntil(target);
		target += nwait;
	}
	timer.stop();
	
	double measDur = timer.elapsedSec();
	printf("\nMeasured time: % 8.6f (% 8.6f%% error)\n", measDur, error(measDur, dur));

	return 0;
}
