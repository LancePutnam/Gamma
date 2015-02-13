#ifndef GAMMA_NODE_H_INC
#define GAMMA_NODE_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <stdio.h>

namespace gam{

/// Doubly-linked node
template <class T>
class Node2{
public:
	Node2();
	Node2(bool zeroLinks);
	//virtual ~Node2();
	~Node2();
	
	T * nodeL;					///< Pointer to left node
	T * nodeR;					///< Pointer to right node
	
	void nodeInsertL(T& node);	///< Insert myself to left of node
	void nodeInsertR(T& node);	///< Insert myself to right of node
	void nodeRemove();			///< Remove myself from linked list

	/// Returns leftmost link
	const Node2<T>& leftmost() const {
		Node2<T> * t = nodeL;
		Node2<T> * n = nodeL;
		while(t){ n = t; t = t->nodeL; }
		return n ? *n : *this;
	}

	/// Returns whether node is linked to other nodes
	bool linked() const;

	void print(const char * append = "\n", FILE * fp = stdout) const;
	void printAll(const char * append = "\n", FILE * fp = stdout) const;
};


/// Triply-linked node
template <class T>
class Node3{
public:

	T * parent;		///< Parent node
	T * child;		///< Child node
	T * sibling;	///< Right sibling

	Node3()
	:	parent(0), child(0), sibling(0)
	{}
	
//	const T * parent() const { return mParent; }
//	const T * child() const { return mChild; }
//	const T * sibling() const { return mSibling; }
//
//	T * parent(){ return mParent; }
//	T * child(){ return mChild; }
//	T * sibling(){ return mSibling; }

	/// Add node as my first child
	void addFirstChild(T * newChild){
		newChild->removeFromParent();
		newChild->parent = self();
		newChild->sibling = child;
		child = newChild;
	}

	/// Add node as my last child
	void addLastChild(T * newChild){
		newChild->removeFromParent();
		newChild->parent = self();
		if(!child){	// No children, so make first child
			child = newChild;
		}
		else{		// Have children, so add to end of children
			lastChild()->sibling = newChild;
		}
	}
	
	/// Remove self from parent leaving my own descendent tree intact
	void removeFromParent(){

		if(parent && parent->child){

			// re-patch parent's child?
			if(parent->child == self()){
				// I'm my parent's first child 
				// - remove my reference, but keep the sibling list healthy
				parent->child = sibling;
			}
			
			// re-patch the sibling chain?
			else{
				// I must be one of parent->child's siblings
				T * temp = parent->child;
				while(temp->sibling){
					if(temp->sibling == self()) {
						// I'm temp's sibling
						// - remove my reference, keep the sibling list healthy
						temp->sibling = this->sibling; 
						break; 
					}
					temp = temp->sibling;
				}
			}
			
			parent=0; sibling=0; // no more parent or sibling, but child is still valid
		}
	}

	T * lastChild(){
		T * n = child;
		while(n->sibling) n = n->sibling;
		return n;
	}

	/// Returns next node using depth-first traversal
	
	/// Returns 0 when the next node equals the terminal node.
	///
	T * next(const T * const terminal){		
		T * n = self();
		
		if(n->child){
			n = n->child;
			return (n != terminal) ? n : 0;
		}
		else{
			return nextBreadth(terminal);
		}
	}

	const T * next(const T * const terminal) const {		
		return const_cast<const T*>(next(terminal));
	}
	
	/// Returns next node using breadth-first traversal
	T * nextBreadth(const T * const terminal){
		T * n = self();
		if(n->sibling){
			n = n->sibling;
		}
		else{
			while(n != terminal && n->sibling == 0){
				n = n->parent;
			}
			if(n != terminal && n->sibling){
				n = n->sibling;
			}
		}
		return (n != terminal) ? n : 0;
	}

private:
	T * self(){ return static_cast<T*>(this); }
	const T * self() const { return static_cast<T*>(this); }
};

// Implementation_______________________________________________________________

// Node2

template <class T>
Node2<T>::Node2()
:	nodeL(0), nodeR(0)
{
}

template <class T>
Node2<T>::Node2(bool zeroLinks){
	if(zeroLinks){
		nodeL = 0;
		nodeR = 0;
	}
}

template <class T>
Node2<T>::~Node2(){
	//printf("~Node2\n");
	nodeRemove();
}

template <class T>
void Node2<T>::nodeInsertL(T & node){
	nodeL = node.nodeL;
	nodeR = &node;
	if(nodeL) nodeL->nodeR = (T *)this;
	node.nodeL = (T *)this;
}

template <class T>
void Node2<T>::nodeInsertR(T & node){
	nodeR = node.nodeR;
	nodeL = &node;
	if(nodeR) nodeR->nodeL = (T *)this;
	node.nodeR = (T *)this;
}

template <class T>
void Node2<T>::nodeRemove(){
	// connect neighbors
	if(nodeL){ nodeL->nodeR = nodeR; }
	if(nodeR){ nodeR->nodeL = nodeL; }
	nodeL = 0;
	nodeR = 0;
}

template <class T>
bool Node2<T>::linked() const {
	return nodeL || nodeR;
}

template <class T>
void Node2<T>::print(const char * append, FILE * fp) const {
	fprintf(fp, "%p %p %p%s", nodeL, this, nodeR, append);
}

template <class T>
void Node2<T>::printAll(const char * append, FILE * fp) const {

	Node2<T> const * n = &leftmost();
	
	while(n){
		fprintf(fp, n == this ? "(%p) " : " %p  ", n);
		n = n->nodeR;
	}
	fprintf(fp, "%s", append);
}

} // gam::

#endif

