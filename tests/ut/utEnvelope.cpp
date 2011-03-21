{
	{
		const int N = 80000;
		const double eps = 1e-4;
		Curve<double> c(N, 0,-0.5);
		assert(!c.done());
		assert(c.end() == -0.5);
		
		for(int i=0; i<N; ++i){
			double v = c();
//			printf("%g ?= %g\n", v, double(i)/N*c.end());
			assert(scl::abs(v - double(i)/N*c.end()) < eps);
			assert(v == c.value());
		}
		
		c(); // cannot hit end point exactly
		assert(c.done());

		c.reset();
		for(int i=0; i<N; ++i){
			double v = c();
			assert(scl::abs(v - double(i)/N*c.end()) < eps);
			assert(v == c.value());
		}
	}

	{
		const int N = 80000;
		const float eps = 1e-1;
		Curve<float> c(N, 0,-0.5);
		assert(!c.done());
		assert(c.end() == -0.5);
		
		for(int i=0; i<N; ++i){
			double v = c();
//			printf("%g ?= %g\n", v, float(i)/N*c.end());
			assert(scl::abs(v - float(i)/N*c.end()) < eps);
			assert(v == c.value());
		}
		
		c(); // cannot hit end point exactly
		assert(c.done());
	}
}
