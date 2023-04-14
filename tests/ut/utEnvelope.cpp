{

	// test double precision
	{
		const int N = 80000;
		const double eps = 1e-4;
		Curve<double,double> c(N, 0, 0.0,-0.5);
		assert(!c.done());
		assert(c.end() == -0.5);
		
		for(int i=0; i<N; ++i){
			double v = c();
//			printf("%g ?= %g\n", v, double(i)/N*c.end());
			assert(aeq(v, double(i)/N*c.end(), eps));
			assert(v == c.value());
		}
		
		c(); c();  // do a couple more iterations to reach endpoint
		assert(c.done());

		// test reset
		c.reset();
		for(int i=0; i<N; ++i){
			double v = c();
			assert(aeq(v, double(i)/N*c.end(), eps));
		}
	}

	
	// test various start and end points
	{
		const float points[4][2] = {{0,1}, {1,0}, {-1,1}, {1,-1}};
		const int N = 8000;
		const float eps = 1e-1;

		for(int j=0; j<4; ++j){
			float start = points[j][0];
			float end   = points[j][1];
			Curve<> c(N, 0, start,end);
			assert(!c.done());
			assert(c.end() == end);
			
			for(int i=0; i<N; ++i){
				double v = c();
				//printf("%g ?= %g\n", v, 100 - float(i)/N*50);
				assert(aeq(v, (float(i)/N*(end-start) + start), eps));
			}
			
			c(); // do one more iteration to reach endpoint
			assert(c.done());
		}
	}
	

	{
		Env<3> e;
		assert(e.size() == 3);
		assert(!e.done());
		assert(!e.released());
		
		e.levels(0, -0.2, 0.5, 0.1);
		e.maxLevel(2);
		assert(aeq(e.levels()[0], 0.0f));
		assert(aeq(e.levels()[1],-0.8f));
		assert(aeq(e.levels()[2], 2.0f));
		assert(aeq(e.levels()[3], 0.4f));
		
		e.lengths(1,2,3);
		assert(aeq(e.totalLength(), 6.f));
		
		e.totalLength(2);
		assert(aeq(e.lengths()[0], 1.f/6 * 2));
		assert(aeq(e.lengths()[1], 2.f/6 * 2));
		assert(aeq(e.lengths()[2], 3.f/6 * 2));

		e.lengths(1,3,1);
		e.totalLength(10, 1);
		assert(aeq(e.lengths()[0], 1.f));
		assert(aeq(e.lengths()[1], 8.f));
		assert(aeq(e.lengths()[2], 1.f));		
	}
	
	{
		//typedef Vec<2,double> vec;
		#define vec Vec<2,double> 
		Env<2, vec> e;
		
		assert(!e.done());
		
		e.levels(vec(0, -0.1), vec(0.5, 1), vec(0.1, 0));
		e.maxLevel(vec(2,3));
		assert(aeq(e.levels()[0][0], 0.0));
		assert(aeq(e.levels()[1][0], 2.0));
		assert(aeq(e.levels()[2][0], 0.4));
		assert(aeq(e.levels()[0][1],-0.3));
		assert(aeq(e.levels()[1][1], 3.0));
		assert(aeq(e.levels()[2][1], 0.0));
	}
}
