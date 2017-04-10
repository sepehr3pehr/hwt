#include "types.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <utility>
#include <stack.h>

#define nullptr 0

UINT64 child_index(UINT8* subnorms, UINT8* child_pattern, int depth) {
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

struct stack_entry {
	Node* node;
	stack_entry *next;
};

struct stack_info {
	stack_entry* head;
	stack_entry* tail;
	UINT64 size;
};


class SearchNode {
private:
	UINT32 B_over_8;
	UINT32 num_found;
	UINT8* query;
	UINT8* code;
	int B;
	int hamm_radius;
	UINT8* query_normchunks;
	UINT32* num_res;
	int max_lvl;
	UINT32 K;
	Node* root;
	int NQ;
	UINT32* num_resutls; //all results
	UINT32* r;
	stack* main_stack;
	stack** stacks_hammd;

	// added variabes for r-solution:
	int max_c;
	int max_r;
	std::pair<int,int> pair_results[100][100][100];
	int num_pairs[100][100];

public:
	void HNN_search();
	SearchNode(UINT8*,UINT32, int, int, Node*,int);
	~SearchNode();
	void setQuery(UINT8* query);
	void all_children(Node*);
	int lower_bound_distance(Node*);
	void linear_scan_node(Node*);
	int num_compare;
	UINT32** results;
	UINT32* res;

	// added functions for r-solution:
	void par_solutions();
	int get_size(int r, int c);
	int num_possible_child(int* hamm_ws, int* r_dis, int two_to_d, int* num_partial_result);
	void check_r_solution(Node*);
	void node_search(entry*, bool, int );






};

SearchNode::SearchNode(UINT8* _code,UINT32 _B_over_8, int _max_lvl, int _K, Node* _root, int _NQ) {

	root = _root;
	max_lvl = _max_lvl;
	B_over_8 = _B_over_8;
	code = _code;
	K = _K;	
	B = B_over_8 * 8;
	NQ = _NQ;
	main_stack = new stack();
	max_c=50;
	max_r=100;
	



	results = (UINT32**) malloc(sizeof(UINT32*)*NQ);
	results[0] = (UINT32 *) malloc(sizeof(UINT32)*K*NQ);
	r = results[0];
	res = new UINT32[K*(B/2)];

	for (size_t i=1; i<NQ; i++)
		results[i] = results[i-1] + K;

	int num_chunks = pow2(max_lvl+1)-1;
	
	query_normchunks = (UINT8*) malloc(num_chunks * sizeof(UINT8));


	num_res = (UINT32*) malloc((B+1) * sizeof(UINT32));
	par_solutions();


}

SearchNode::~SearchNode() {
	free(results[0]);
	free(results);
	free(num_res);
	free(query_normchunks);

	delete[] res;

	for(int i=0;i<B+1;i++)
		delete stacks_hammd[i];
	delete[] stacks_hammd;

}
void SearchNode::setQuery(UINT8* _query) {
	query = _query;
	num_compare = 0;
	num_found = 0;

	stacks_hammd = new stack*[B+1];
	for(int i=0;i<B+1;i++)
		stacks_hammd[i] = new stack();



	memset(num_res, 0, (B+1)*sizeof(*num_res));

	UINT8* c_query_normchunks = query_normchunks;	 
	for(int i=0;i<max_lvl;i++) {
		norm_chunks(c_query_normchunks,(UINT32)i, query,(int)B_over_8);
		c_query_normchunks += (int)pow(2,i);
	}
}

void SearchNode::HNN_search() {
	num_found = 0;
	hamm_radius = 0;
	entry* cur = nullptr;
	UINT32 total_offset;

	entry* current_entry;
	sparse_hash_map<UINT64,Node*>::const_iterator got;

	bool all_chilld = true;
	for(int i=0;i<=B;i++) {
		if(query_normchunks[0]+i <= B) {
			total_offset = query_normchunks[0]+i;
			got = root->children.find(total_offset);
			if(got != root->children.end()) {
				Node* h = got->second;
				stacks_hammd[i]->push(h);
				if(lower_bound_distance(h) != i)
					printf("error");
			}
			if(i==0)
				continue;
		}
		if(query_normchunks[0]-i>=0) {
			total_offset = query_normchunks[0]-i;
			got = root->children.find(total_offset);
			if(got != root->children.end()) {
				Node* h = got->second;
				if(lower_bound_distance(h) != i)
					printf("error");
				stacks_hammd[i]->push(h);
			}


		}
	}
	while (true) {
		//printf("Hamming radius = %d\n",hamm_radius);
		

		for(int stack_iter = hamm_radius%2 ;stack_iter <=hamm_radius;stack_iter = stack_iter + 2) {
			cur = stacks_hammd[stack_iter]->head;
			
			while(cur!=nullptr) {
				if(false)
					all_chilld =true;
				else
					all_chilld = false;
				node_search(cur,all_chilld,stack_iter);
				if(all_chilld) { //  The node is no longer needed in the linked list since we have checked all of its children
					entry* ent = cur;
					cur = cur->next;
					stacks_hammd[stack_iter]->size--;
					if(ent->next==nullptr && ent->prev==nullptr)
					{
						stacks_hammd[stack_iter]->head = nullptr;
						stacks_hammd[stack_iter]->tail = nullptr;
						delete ent;

					}
					else if(ent->next==nullptr) {//last element of linked list
						stacks_hammd[stack_iter]->tail = ent->prev;
						ent->prev->next = nullptr;
						delete ent;

					}
					else if(ent->prev == nullptr) {//head element of linked list
						stacks_hammd[stack_iter]->head = ent->next;
						ent->next->prev = nullptr;
						delete ent;

					}
					else {
						ent->next->prev = ent->prev;
						ent->prev->next = ent->next;
						delete ent;
					}
				}

				else {
					if(cur->current->isleaf) { // for the first elements of stacks that are leaf (only these can be leaf)
						stacks_hammd[stack_iter]->size--;
						entry* ent = cur;
						if(ent->next==nullptr && ent->prev==nullptr)
						{
							stacks_hammd[stack_iter]->head = nullptr;
							stacks_hammd[stack_iter]->tail = nullptr;
							delete ent;

						}
						else if(ent->next==nullptr) {//last element of linked list
							stacks_hammd[stack_iter]->tail = ent->prev;
							ent->prev->next = nullptr;
							delete ent;

						}
						else if(ent->prev == nullptr) {//head element of linked list
							stacks_hammd[stack_iter]->head = ent->next;
							ent->next->prev = nullptr;
							delete ent;

						}
						else {
							ent->next->prev = ent->prev;
							ent->prev->next = ent->next;
							delete ent;
						}
					}
					cur = cur->next;
				}



			}
		}
		//break;
		num_found += num_res[hamm_radius];
		//printf("num_found = %d\n",num_found);

		if(num_found + num_res[hamm_radius+1] >= K )
		{
			//printf("Broke at %d\n",hamm_radius);
			//printf("num_found %d\n",num_found);
			break;
		}
		hamm_radius++;
		
		if(hamm_radius>B/2)
		{	//printf("found so-far: %d\n",num_found + num_res[hamm_radius+1]);
			break;
		}


	}
	int n=0;
	for(int s=0; s<=B/2 && n < K; s++)
		for(int c = 0; c < num_res[s] &&  n < K; c++) {
			//printf("Nearest neighbor: %d\n",res[s*K+c]);
			//results[n++] = res[s*K+c];

		}
	
	//r += K;



}

void SearchNode::node_search(entry* ent, bool all_child,int stack_index) {
	if(ent->current->isleaf) {
		linear_scan_node(ent->current);
		return;
	}
	if(all_child) {
		all_children(ent->current);
	}
	else
		check_r_solution(ent->current);

}

int SearchNode::lower_bound_distance(Node* node) {
	int s_index = pow2(node->depth)-1;
	int e_index = pow2(node->depth+1)-1;
	UINT8* c_norm_chunk = query_normchunks + s_index;
	int range = e_index - s_index;
	int sum_diff = 0;
	for(int i=0;i<range;i++)
		sum_diff += abs(node->subnorms[i] - c_norm_chunk[i]);
	return sum_diff;
}

void SearchNode::linear_scan_node(Node* node) {
	int hammd;

	if(lower_bound_distance(node) == hamm_radius) {
		
		LL_node* curr= node->head_list_codes;
		
		while(curr != nullptr) {

			

			hammd = match(query,code + (UINT64)curr->index*(B_over_8), B_over_8);

			num_compare++;
			if(hammd<B/2 && num_res[hammd] < K) {
				res[hammd * K + num_res[hammd]] = curr->index + 1;
			}
			num_res[hammd]++;
			if(hammd == hamm_radius) {
				num_found++;
				if(num_found == K)
					break;
			}
			curr = curr->next;
		}
	}
	else
		printf("error2, hamm radius is%d but lower is %d, and depth is %d\n",hamm_radius,lower_bound_distance(node),node->depth);
	return;

}


void SearchNode::all_children(Node* node) {
	for ( auto it = node->children.begin(); it != node->children.end(); ++it ) {
		Node* temp = it->second;
		int a = lower_bound_distance(temp);
		if(a==hamm_radius && temp->isleaf)
			linear_scan_node(temp);
		else if(a>=hamm_radius)
			stacks_hammd[a]->push(temp);
	}

}

int SearchNode::get_size(int r, int c) {
	return num_pairs[r][c+max_c];
}




void SearchNode::check_r_solution(Node* node) {
	int r = hamm_radius;
	UINT64 total_offset;
	sparse_hash_map<UINT64,Node*>::const_iterator got;
	int two_to_d = pow2(node->depth);
	UINT8* query_pattern = query_normchunks + two_to_d - 1;
	UINT8* query_pattern_next_lvl = query_normchunks + 2*two_to_d - 1;
	UINT8* node_pattern = node->subnorms;
	int hamm_ws[two_to_d];
	int num_partial_result[two_to_d];
	UINT8 child_pattern[2*two_to_d];
	int s = two_to_d - 1;
	int r_dis[two_to_d]; // distribution of radius
	int bit = s-1;
	int power[100];
	int r_index;
	int c_index;

	bool valid_child = true;


	// used for stopping criterion (location of (s+1)th 1)

	int sum_diff = 0;
	for(int i=0; i<two_to_d;i++) {
		hamm_ws[i] = node_pattern[i] - query_pattern[i];
		sum_diff += abs(hamm_ws[i]);
	}
	r = (hamm_radius - sum_diff);

	if(r<0) {
		printf("negative r\n");
		return;
	}

	if(r%2 != 0) {
		printf("odd r\n");
		return;
	}
	r = r/2;
	

	bit = s-1;
	for (int i=0; i<s; i++)
		power[i] = i;			// power[i] stores the location of the i'th 1
	power[s] = s+r+1;
	while (true) {
		// the loop to distribute radius
		if (bit != -1) {
			
			power[bit]++;
			bit--;
		}
		else {
			r_dis[0] = 2 * (power[0]-1)+abs(hamm_ws[0]);
			for (int i=1;i<s;i++) {
				r_dis[i] = 2 * (power[i]-power[i-1]-1)+abs(hamm_ws[i]);
			}
			if(s!=0)
				r_dis[s] = 2 * (power[s] - power[s-1] - 1)+abs(hamm_ws[s]);

			int count_child = num_possible_child(hamm_ws, r_dis, two_to_d, num_partial_result);
			
			if(count_child != 0) {
				
				int child_to_check = count_child;
				int j=0;
				int index = 0;
				while(j<count_child) {
					
					int temp = j;
					valid_child = true;
					for(int k=0; k<2*two_to_d;k = k+2) {
						
						index = temp % num_partial_result[k/2];
						temp = temp / num_partial_result[k/2];
						
						r_index = r_dis[k/2];
						c_index = hamm_ws[k/2]+max_c;
						child_pattern[k] = pair_results[r_index][c_index][index].first + query_pattern_next_lvl[k];
						child_pattern[k+1] = pair_results[r_index][c_index][index].second + query_pattern_next_lvl[k+1];
						if(child_pattern[k] <0 || child_pattern[k] > (B_over_8*8)/pow2(node->depth) || child_pattern[k+1] <0 ||
						 child_pattern[k+1] > (B_over_8*8)/pow2(node->depth)) {
							valid_child = false;
							break;
						}
					}
					j++;
					if(!valid_child)
						continue;
					
					total_offset = child_index(node->subnorms, child_pattern,node->depth);
					got = node->children.find(total_offset);
					if(got != node->children.end()) {
						Node* h = got->second;
						int cq = lower_bound_distance(h);
						if(cq != hamm_radius) {
							
							total_offset = child_index(node->subnorms, child_pattern,node->depth);
						}
						if(h->isleaf)
							linear_scan_node(h);
						else
							stacks_hammd[hamm_radius]->push(h);
					}
				}
			}
			while (++bit < s && power[bit] == power[bit+1]-1) {
				power[bit] = bit;
			}
			if (bit == s)
				break;

		}
	}
	

}

int SearchNode::num_possible_child(int* hamm_ws, int* r_dis, int two_to_d, int* num_partial_result) {
	int num_res = 1;
	for (int i=0;i<two_to_d;i++) {
		
		num_res *= get_size(r_dis[i],hamm_ws[i]);
		num_partial_result[i] = get_size(r_dis[i],hamm_ws[i]);
		if(num_res == 0)
			return 0;

	}
	return num_res;
}

void SearchNode::par_solutions() {
	for (int r=0;r<max_r;r++)
		for(int c=0;c<2*max_c;c++)
			num_pairs[r][c] = 0;
	for (int r=0;r<max_r;r++)
	{
		
		for (int c=-1*max_c;c<max_c;c++) {

			if(r < abs(c) || r%2 != abs(c)%2)
				continue;

			// case 1:
			if( c == r) {
				for( int temp = 0;temp<=c;temp++) {
					pair_results[r][c+max_c][num_pairs[r][c+max_c]] = std::make_pair(temp,c-temp);
					num_pairs[r][c+max_c]++;
				}
				continue;
			}
			//case 4:
			if(c == -1*r) {
				for(int temp = c;temp<=0;temp++) {
					pair_results[r][c+max_c][num_pairs[r][c+max_c]] = std::make_pair(temp,c-temp);
					num_pairs[r][c+max_c]++;
				}
				continue;
			}

			// case 2 and 3:
			pair_results[r][c+max_c][num_pairs[r][c+max_c]] = std::make_pair((c-r)/2,(c+r)/2);
			num_pairs[r][c+max_c]++;
			pair_results[r][c+max_c][num_pairs[r][c+max_c]] = std::make_pair((c+r)/2,(c-r)/2);
			num_pairs[r][c+max_c]++;


		}


	}
}



