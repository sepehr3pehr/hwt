/*
 * node.h
 *
 *  Created on: 2016-03-09
 *      Author: s2eghbal
 */

#ifndef NODE_H_
#define NODE_H_

#include <stdio.h>
#include <iostream>
#include <unordered_map>
#include <string>
#include "types.h"
#include <LL_node.h>
#include <sparsepp.h>

using namespace std;
using spp::sparse_hash_map;
class Node {

private:
	//UINT64 findchild(UINT8* code);

public:

	int level;

	static int B_over_8;
	static UINT8* dbcode;
	UINT64 num_children;
	UINT64* assigned_codes;

	bool isleaf; // When a node is created it is initially a leaf


	UINT64 branchfactor;

	LL_node* head_list_codes;
	LL_node* tail_list_codes;


	int size;



	int depth; // Depth of the node

	UINT8* subnorms;

	sparse_hash_map<UINT64,Node*> children;

	Node(int _depth,UINT8* chunk_subnorms);

	~Node();

	void expand(); // Expand if the number nodes is more than a threshold

 //  void insert(UINT8* _code, int _index, int B_over_8); // It should exapand if the number of nodes is more than the threshold
   Node* find_the_child(UINT8* code,UINT8* subnorms);
   UINT64 find_offset_of_child(UINT8* code,UINT8* subnorms);
   void push_to_node(UINT8* _code, UINT32 index);

};




#endif /* NODE_H_ */
