{
	{
		using namespace acc;

		assert(None::map(0,2,1) == 0);
		assert(Wrap::map(0,2,1) == 2);
		assert(Clip::map(0,2,1) == 1);

		assert(None::map(2,2,1) == 2);
		assert(Wrap::map(2,2,1) == 2);
		assert(Clip::map(2,2,1) == 2);

		assert(None::map(3,2,1) == 3);
		assert(Wrap::map(3,2,1) == 1);
		assert(Clip::map(3,2,1) == 2);
	}
}
