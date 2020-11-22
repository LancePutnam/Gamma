{
	{	// Test onResize
		struct TestArray : public Array<int>{
			TestArray(): resized(false){}

			void onResize(){
				resized = true;
			}
			
			bool resized;
		};
	
		TestArray a;
		assert(!a.resized);	// default ctor should not resize

		a.resized = false;	// onResize should be called whenever the size changes
		a.resize(16);
		assert(a.resized);
		a.resized = false;	
		a.resize(8);
		assert(a.resized);
		
		a.resized = false;	// onResize should NOT be called since the size didn't change
		a.resize(8);
		assert(!a.resized);
	}

	{
		const int N=16;
		typedef int t;
		typedef Array<t> array_t;
		auto * a = new array_t(N, 123);
		auto * b = new array_t;
		b->source(*a);

		for(unsigned i=0; i<a->size(); ++i) assert((*a)[i] == 123);
		assert(a->elems() == b->elems());
		assert(a->size() == b->size());
		assert(array_t::references(a->elems()) == 2);

		a->clear();
		assert(a->size() == 0);
		assert(a->elems() == 0);

		delete a;
		assert(array_t::references(b->elems()) == 1);

		{ // setting from raw pointer should not manage
			auto * c = new array_t(b->elems(), b->size());
			assert(array_t::references(b->elems()) == 1);
		}

		delete b;
		assert(array_t::references(b->elems()) == 0);

		a = new array_t(N, 123);
		b = new array_t(*a);
		
		b->own();
		assert(a->elems() != b->elems());
		assert(array_t::references(a->elems()) == 1);
		assert(array_t::references(b->elems()) == 1);
		for(unsigned i=0; i<a->size(); ++i) assert((*b)[i] == 123);
		
		t * elemsA = a->elems();
		t * elemsB = b->elems();
		a->source(*b);
		assert(a->elems() == b->elems());
		assert(array_t::references(elemsA) == 0);
		assert(array_t::references(elemsB) == 2);
	}


	{
		DelayN<int> d(2);
	
		//for(int i=0; i<8; ++i) printf("%d\n", d(i+1));
	
		assert(d(1) == int());
		assert(d(2) == int());
		assert(d(3) == 1);
		assert(d(4) == 2);
		assert(d(5) == 3);
	}
	
//	{ Array<t> a(N); }
//	{ ArrayPow2<t> a(N); }
//	{ Ring<t> a(N); }
//	{ DoubleRing<t> a(N); a.copy(); }
}
