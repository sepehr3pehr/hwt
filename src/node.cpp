#include "node.h"
#include "types.h"
#include "bitops.h"
#include <stdlib.h>
#include <iostream>
#include <unordered_map>
#include "sparsepp.h"
#include <LL_node.h>

UINT64 find_offset_of_child3(UINT8* subnorms,UINT8* child_pattern, int depth) {
	UINT64 p2d = pow2(depth);
	UINT64 offsets[p2d];
	UINT64 index = 0;
	offsets[0] = 1;
	for(UINT64 i=1;i<p2d;i++)
		offsets[i] = (subnorms[i-1]+1) * offsets[i-1];
	for(UINT64 i=0;i<p2d;i++)
		index += child_pattern[2*i]*offsets[i];
	return index;
}


Node::Node(int _depth,UINT8* chunk_subnorms) {
	isleaf = true;
	depth = _depth;
	size = 0;
	num_children = 0;
	if(depth < 0)
		isleaf = false;
	if(depth >= 0) {
		subnorms = (UINT8*)calloc(pow(2,depth),sizeof(UINT8));
		head_list_codes = tail_list_codes = NULL;
		
		for(int i = 0;i<pow(2,depth);i++) 
			subnorms[i] = chunk_subnorms[i];
 	}
	branchfactor = 0;

}

Node::~Node() {
	if(this==NULL)
		return;
	if(!isleaf)
	{

		for ( auto it = children.begin(); it != children.end(); ++it )
			delete it->second;
	}

	else {
		LL_node* it,*it2 = head_list_codes;
		it=it2;
		for(int i=0;i<size;i++) {
			it2 = it->next;
			delete it;
			it=it2;
		}

	}
}

void Node::expand() {
	if(!isleaf)
	{
		printf("Error: calling expand on a non-leaf node");
		return;
	}
	isleaf = false;



	UINT64  total_offset=0;
	UINT8* subsubnorms = (UINT8*) calloc(pow(2,depth+1),sizeof(UINT8));
	UINT8* code;
	Node* targetchild;
	LL_node* iterat = head_list_codes;
	LL_node* iterat2;
	for(int i=0;i<size;i++){
		code = dbcode+(UINT64)iterat->index*(UINT64)B_over_8;
		
		iterat2 = iterat;
		iterat = iterat->next;

		
		norm_chunks(subsubnorms,depth+1,code,B_over_8);
		total_offset = find_offset_of_child(code, subnorms);
		sparse_hash_map<UINT64,Node*>::const_iterator got = children.find(total_offset);
		if(got==children.end())
		{

			Node* newchild = new Node(this->depth+1,subsubnorms);
			children.insert({total_offset,newchild});
			targetchild = newchild;
			num_children++;

		}
		else{
			
			targetchild = got->second;
			
		}
		targetchild->push_to_node(code,iterat2->index);
		delete iterat2;

	}
	delete[] subsubnorms;
}


Node* Node::find_the_child(UINT8* code,UINT8* subnorms)
{
	Node* target;
	UINT8* subsubnorms = (UINT8*) calloc(pow(2, depth+1),sizeof(UINT8));
	norm_chunks(subsubnorms,depth+1,code,B_over_8);

	UINT64 total_offset = find_offset_of_child(code, subnorms);
	sparse_hash_map<UINT64,Node*>::const_iterator got = children.find(total_offset);
	if(got==children.end()){
		target = new Node(depth+1,subsubnorms);
		this->num_children++;
		children.insert({total_offset,target});
	}
	else {
		target = got->second;
		}
	return target;


}



UINT64 Node::find_offset_of_child(UINT8* code,UINT8* subnorms) {
	if(depth == -1) //root node
		return subnorms[0];
	int numsubs = pow(2,depth);
	UINT64  total_offset=0;
	UINT64 offsets[(UINT64)pow2(depth)];
	UINT8* subsubnorms = (UINT8*) calloc(pow2(depth+1),sizeof(UINT8));
	norm_chunks(subsubnorms,depth+1,code,B_over_8);

	
	UINT64 temp = find_offset_of_child3(this->subnorms,subsubnorms,depth);
	
	temp = find_offset_of_child3(this->subnorms,subsubnorms,depth);
	delete subsubnorms;
	return temp;
}




void Node::push_to_node(UINT8* _code, UINT32 index) {
	if(!isleaf){
		printf("error: calling push on an internal node\n");
		return;
	}
	else{
		if(head_list_codes == NULL){
			head_list_codes = tail_list_codes = new LL_node;
			head_list_codes->next = NULL;
			tail_list_codes->index = index;
			size++;
		}

		else{
			tail_list_codes->next = new LL_node;
			tail_list_codes = tail_list_codes->next;
			tail_list_codes->next = NULL;
			tail_list_codes->index = index;
			size++;
		}
	}
}

//
//UINT64 Node::findchild(UINT8* code) {
//	UINT64  total_offset=0;
//	UINT64 offsets[(UINT64)pow(2,depth)];
//	UINT8* subsubnorms = (UINT8*) calloc(pow(2,depth+1),sizeof(UINT8));
//	int numsubs = pow(2,depth);
//	norm_chunks(subsubnorms,depth+1,code,B_over_8);
//	offsets[numsubs-1] = 1;
//	for(int i=numsubs-2;i>=0;i--)
//		offsets[i] = offsets[i+1]*subnorms[i];
//	norm_chunks(subsubnorms,depth+1,code,B_over_8);
//	for(int j=0;j<pow(2,depth);j++)
//		total_offset += subsubnorms[2*j]*offsets[j];
//
//	return total_offset;
//}

