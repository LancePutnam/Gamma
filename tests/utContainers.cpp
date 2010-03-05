#include "ut.h"

int utContainers(){

	typedef int t;
	typedef Array<t> array_t;
	array_t * a = new array_t(16, 123);
	array_t * b = new array_t(*a);

	for(uint32_t i=0; i<a->size(); ++i) assert((*a)[i] == 123);
	assert(a->elems() == b->elems());
	assert(a->size() == b->size());
	assert(array_t::references(a->elems()) == 2);

	delete a;
	assert(array_t::references(b->elems()) == 1);
	
	array_t * c = new array_t(b->elems(), b->size());
	assert(array_t::references(b->elems()) == 2);

	delete b;
	assert(array_t::references(c->elems()) == 1);
	
	t * elemsC = c->elems();
	delete c;
	assert(array_t::references(elemsC) == 0);
	
	a = new array_t(16, 123);
	b = new array_t(*a);
	
	b->own();
	assert(a->elems() != b->elems());
	assert(array_t::references(a->elems()) == 1);
	assert(array_t::references(b->elems()) == 1);
	for(uint32_t i=0; i<a->size(); ++i) assert((*b)[i] == 123);
	
	t * elemsA = a->elems();
	t * elemsB = b->elems();
	a->source(*b);
	assert(a->elems() == b->elems());
	assert(array_t::references(elemsA) == 0);
	assert(array_t::references(elemsB) == 2);
	
	return 0;	
}
