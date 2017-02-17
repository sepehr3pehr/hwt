#include "root_node.h"
#include "bitops.h"
#include <stdlib.h>

RootNode::RootNode(int _B, int _capacity) {
	B = _B;
	capacity = _capacity;
	B_over_8 = B/8;
	children = (Node**) calloc(B,sizeof(Node*));
	for(int i=0;i<B;i++)
		children[i] = NULL;
}

bool RootNode::insert(UINT8* code,Node* &curr_node) {

	UINT8 norm[1];
	norm[0] = l1norm(code,B_over_8);
	//printf("norm = %u\n",norm[0]);
	if(children[norm[0]]==NULL){
		children[norm[0]] = new Node(0,norm);
		curr_node = children[norm[0]];
		return true;
	}
	curr_node = children[norm[0]];
	return false;



}

RootNode::~RootNode(){
	for(int i=0;i<B;i++)
		if(children[i] != NULL)
			delete children[i];
	delete[] children;
	printf("Done deleting ROOT node\n");
}
