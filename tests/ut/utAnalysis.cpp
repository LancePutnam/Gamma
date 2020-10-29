{
	{
		PCounter pcounter(4);
		assert(pcounter.count() == 0);
		assert(pcounter() == false);
		assert(pcounter() == false);
		assert(pcounter() == false);
		assert(pcounter.cycled() == false);
		assert(pcounter() == true);
		assert(pcounter.cycled() == true);
		assert(pcounter() == false);
		assert(pcounter.cycled() == false);
	}
	{
		MaxAbs<> maxAbs(4);
		maxAbs(0);
		maxAbs(-1);
		maxAbs(-2);
		maxAbs(1);
		assert(maxAbs.value() == 2);
		maxAbs(4);
		assert(maxAbs.value() == 2);
	}
	{
		ZeroCross<> zeroCross(0);
		assert(zeroCross( 0) ==  0);
		assert(zeroCross( 1) ==  1);
		assert(zeroCross(.5) ==  0);
		assert(zeroCross(-1) == -1);
		assert(zeroCross( 0) ==  0);
		assert(zeroCross( 1) ==  1);
	}
	{
		ZeroCrossRate<> zcr(4);
		zcr( 1);
		zcr(-1);
		zcr( 1);
		zcr(-1);
		assert(zcr.value() ==  1.);
		zcr(-1);
		zcr( 1);
		zcr( 1);
		zcr(-1);
		assert(zcr.value() ==  0.5);
	}

}
