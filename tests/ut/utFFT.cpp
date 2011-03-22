{
	const int N=16;

	Complex<double> ac[N];
	double ar[N];

	CFFT<double> cfft(N);
	RFFT<double> rfft(N);

	fil::Reson<double> r1(1./N), r2(2./N);
	for(int i=0; i<N; ++i){ ac[i] = r1()+r2()+1; ar[i] = ac[i].r; }

//	for(int i=0; i<N; ++i){
//		printf("[%3d] % f % f\n", i, ac[i].r, ac[i].i);
//	}
	
	cfft.forward((double*)ac);
	rfft.forward(ar);

//	for(int i=0; i<N; ++i) printf("[%3d] % f % f\n", i, ac[i].r, ac[i].i);
//	
//	printf("\n");
////	for(int i=1; i<N-1; i+=2) printf("[%3d] % f % f\n", i, ar[i], ar[i+1]);
//	for(int i=0; i<N; i+=2) printf("% f\n% f ", ar[i], ar[i+1]);

//	{
//		for(int i=0; i<N; ++i){
//			//buf[i+1] = i?0:1;
//			buf[i+1] = 1 + cos(float(i)/N * M_2PI) + (i&1?-1:1);
//			//buf[i+1] = (i&1?-1:1);
//		}
//		fft.forward(buf, true, true);
//		for(int i=0; i<N+2; i+=2) printf("%6.3f %6.3f\n", buf[i], buf[i+1]);
//		fft.inverse(buf, true, true);
//		for(int i=0; i<N; ++i) printf("%g\n", buf[i+1]);
//	}
}
