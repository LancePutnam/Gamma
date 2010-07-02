#include <assert.h>
#include "Gamma/Gamma.h"
#include "Gamma/Access.h"
#include "Gamma/Conversion.h"
#include "Gamma/File.h"
#include "Gamma/Strategy.h"
#include "Gamma/Types.h"
#include "Gamma/Thread.h"
#include "Gamma/UnitMapper.h"
#include <map>

using namespace gam;

THREAD_FUNCTION(threadFunc){
	*(int *)user = 1; return NULL;
}


int main(int argc, char* argv[]){

	// File I/O
	{
		const char * path = "test.txt";
		File f(path, "w");
		assert(f.open());
		
		char buf[] = {'H','e','l','l','o',' ','W','o','r','l','d','!'};
		assert(f.write(buf, 1, sizeof(buf)) == sizeof(buf));
		f.close();
		assert(!f.opened());
		
		assert(File::exists(path));
		
		f.mode("r");
		assert(f.open());
		
		char * bufR = f.readAll();
		for(int i=0; i<f.size(); ++i) assert(buf[i] == bufR[i]);
		f.close();
	}
	
	// FunctionTable
	{
		const int N=4;
		FunctionTable<double, ipl::Linear, acc::Wrap> ft(N);
		for(int i=0; i<N; ++i) ft[i]=i;
		
		assert(ft(0./N) == 0);
		assert(ft(1./N) == 1);
		
		assert(ft(0.5/N) == 0.5);
		assert(ft(1.5/N) == 1.5);
	}
	

	// Thread
	{	int x=0;
		Thread t(threadFunc, &x);
		t.wait();
		assert(x == 1);
	}

	return 0;
}
