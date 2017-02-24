#include <stdlib.h>
#include <iostream>
#include "types.h"
#include <unordered_map>
#include <math.h>
#include <stdio.h>
#include "bitops.h"
#include "sparsepp.h"

#ifndef NODEOPS_H_
#define NODEOPS_H_

UINT64 find_offset_of_child(UINT8* code,Node* curr_node, UINT8* subnorms, int B_over_8) {
	int depth = curr_node->depth;
	int numsubs = pow(2,depth);
	UINT64  total_offset=0;
	UINT64 offsets[(UINT64)pow(2,depth)];
	UINT8* subsubnorms = (UINT8*) calloc(pow(2,depth+1),sizeof(UINT8));
	norm_chunks(subsubnorms,depth+1,code,B_over_8);
	offsets[numsubs-1] = 1;
	for(int i=numsubs-2;i>=0;i--)
		offsets[i] = offsets[i+1]*subnorms[i];
	norm_chunks(subsubnorms,depth+1,code,B_over_8);
	for(int j=0;j<pow(2,depth);j++)
		total_offset += subsubnorms[2*j]*offsets[j];

	return total_offset;
}

Node* find_the_child(Node* curr_node,UINT8* code,UINT8* subnorms, int B_over_8,int capacity)
{
	Node* target;
	UINT8* subsubnorms = (UINT8*) calloc(pow(2,curr_node->depth+1),sizeof(UINT8));
	norm_chunks(subsubnorms,curr_node->depth+1,code,B_over_8);

	UINT64 total_offset = find_offset_of_child(code, curr_node, subnorms, B_over_8);
	sparse_hash_map<UINT64,Node*>::const_iterator got = curr_node->children.find(total_offset);
	if(got==curr_node->children.end()){
		target = new Node(curr_node->depth+1,subsubnorms);
	}
	else
		target = got->second;
	return target;


}



void push_to_node(UINT8* _code, Node* curr_node, int B_over_8, UINT32 index) {
	if(!curr_node->isleaf){
		printf("error: calling push on an internal node\n");
		return;
	}
	else{
		UINT8* codet = _code;
		curr_node->tail_list_codes->next = new LL_node;
		curr_node->tail_list_codes = curr_node->tail_list_codes->next;
		curr_node->tail_list_codes->index = index;
		curr_node->size++;
	}
}
#endif
