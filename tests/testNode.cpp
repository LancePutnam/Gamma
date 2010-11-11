#include "Gamma/Node.h"

using namespace gam;

class MyNode2 : public Node2<MyNode2>{
public:
	int v;
};

// Using Node4 as a mix-in.
class MyNode4 : public Node4<MyNode4>{
public:
	int v;
	
	// Override Node4's print()
	void print(char * append){ printf("%d", v); printf("%s",append); }
};



int main(int argc, char ** argv){

	const int numNode2 = 4;
	MyNode2 n2[numNode2];
	
	#define PRINT_N2 \
		for(int i=0; i<numNode2; i++){\
			printf("[%p %p %p]\n", n2[i].nodeL, &n2[i], n2[i].nodeR);\
		}

	printf("\nInsert to left:\n");
	for(int i=numNode2-1; i>0; i--) n2[i-1].nodeInsertL(n2[i]);
	PRINT_N2
	
	for(int i=0; i<numNode2; i++) n2[i].nodeRemove();

	printf("\nInsert to right:\n");
	for(int i=0; i<numNode2-1; i++) n2[i+1].nodeInsertR(n2[i]);
	PRINT_N2

	printf("\nRemove n2[1]:\n");
	n2[1].nodeRemove();
	PRINT_N2
	
	printf("\nTraverse to right (head start):\n");
	MyNode2 * n2t = &n2[0];
	while(n2t){
		printf("%p ", n2t);
		n2t = n2t->nodeR;
	}

	printf("\n\nTraverse to left (tail start):\n");
	n2t = &n2[numNode2-1];
	while(n2t){
		printf("%p ", n2t);
		n2t = n2t->nodeL;
	}	



	printf("\n\n");
	MyNode4 * n1a = new MyNode4();
	MyNode4 * n2a = new MyNode4();
	MyNode4 * n2b = new MyNode4();
	MyNode4 * n3a = new MyNode4();

	n1a->add(n2a);
	n1a->add(n2b);
	n2a->add(n3a);
	
	// change member directly
	n1a->v = 4;
	n2a->v = 3;
	n2b->v = 2;
	n3a->v = 1;
	
	n1a->printDescendents("\n");
	
	// change member using iterator
	MyNode4 * node = n1a;
	int v = 1;
	
	do{
		node->v = v++;
		node = node->next(n1a);
	}while(node);
	
	n1a->printDescendents("\n");	

	n2b->setAsFirstChild();
	n1a->printDescendents("\n");

	n2b->setAsLastChild();
	n1a->printDescendents("\n");
	
	n2a->remove();
	n1a->printDescendents("\n");

	return 0;
}

