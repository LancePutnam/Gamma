#include <stdio.h>
#include "Particle.h"

using namespace gam;
using namespace gam::gen;

int main(int argc, char* argv[]){

	PointParticles<> p(4, 3, 3);	// 4 particles, 3 dimensions, 3 rates

	printf("\nInitial values...\n");
	p.print();

	printf("\nSet positions...\n");
	mem::set(p.pos(0), p.size(), 1.f);
	mem::set(p.pos(1), p.size(), 2.f);
	mem::set(p.pos(2), p.size(), 3.f);
	p.print();

	printf("\nSet velocities...\n");
	mem::set(p.vel(0), p.size(), 0.1f);
	mem::set(p.vel(1), p.size(), 0.2f);
	mem::set(p.vel(2), p.size(), 0.3f);
	p.print();

	printf("\nSet accelerations...\n");
	mem::set(p.acc(0), p.size(), 0.01f);
	mem::set(p.acc(1), p.size(), 0.02f);
	mem::set(p.acc(2), p.size(), 0.03f);
	p.print();

	printf("\nRun update...\n");
	p();
	p.print();

	printf("\nRun update...\n");
	p();
	p.print();
	
	return 0;

}
