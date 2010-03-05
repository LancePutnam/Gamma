#include <assert.h>

#include "ut.h"

int main(int argc, char* argv[]){

	// Unit tests are ordered from the least to most dependent functions/objects
	// in order to catch errors in base functionality.

	utTypes();
	utConversion();
	utContainers();
	utAccess();

	ut_fil();
	ut_gen();
	ut_mem();
	ut_scl();	
	ut_ipl();
	
	ut_arr();

	return 0;
}
