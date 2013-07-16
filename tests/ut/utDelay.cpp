{
	Delay<float, ipl::Trunc, Domain1> delay(4);
	//assert(delay.size() == 4);

	delay.zero();

	// Test normal filtering (delaying) operation
	assert(0 == delay(1));
	assert(0 == delay(2));
	assert(0 == delay(3));
	assert(0 == delay(4));
	assert(1 == delay(5));
	assert(2 == delay(6));

	assert(6 == delay.read(1));
	assert(5 == delay.read(2));
	assert(4 == delay.read(3));
	assert(3 == delay.read(4));
}
