{
{
	Biquad<float, float, Domain1> fil;
		
}

{
	MovingAvg<> fil(4);
	assert(aeq(fil(1), 0.25));
	assert(aeq(fil(1), 0.50));
	assert(aeq(fil(1), 0.75));
	assert(aeq(fil(1), 1.00));
	assert(aeq(fil(1), 1.00));
}

}
