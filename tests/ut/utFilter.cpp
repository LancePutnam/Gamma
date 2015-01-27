{
{
	Biquad<float, float, Domain1> fil;
		
}

{
	MovingAvg<> fil(4);
	assert(near(fil(1), 0.25));
	assert(near(fil(1), 0.5 ));
	assert(near(fil(1), 0.75));
	assert(near(fil(1), 1.  ));
	assert(near(fil(1), 1.  ));
}

}
