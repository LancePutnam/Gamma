{
	struct TestObserver : DomainObserver{

		double checkSPU;

		TestObserver(): checkSPU(0){}

		void onDomainChange(double r){
			checkSPU = spu();
		}

	};

	Domain domA, domB;
	TestObserver obs1, obs2, obs3;

	domA.spu(10);

	assert( 0 == obs1.checkSPU);
	assert( 0 == obs2.checkSPU);
	assert( 0 == obs3.checkSPU);

	// Check observer notification upon attaching to a new subject
	domA << obs1 << obs2 << obs3;

	assert(10 == obs1.checkSPU);
	assert(10 == obs2.checkSPU);
	assert(10 == obs3.checkSPU);

	// Check observer notification upon changing subject
	domA.spu(100);

	assert(100 == obs1.checkSPU);
	assert(100 == obs2.checkSPU);
	assert(100 == obs3.checkSPU);

	// Check change of subjects
	domB.spu(200);
	domB << obs1 << obs2;
	domA.spu(10);

	assert(200 == obs1.checkSPU);
	assert(200 == obs2.checkSPU);
	assert(10 == obs3.checkSPU);
}
