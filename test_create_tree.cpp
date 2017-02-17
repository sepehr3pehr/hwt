#include <stdio.h>
#include "root_node.h"
//#include "node.h"
#include "io.h"
#include <time.h>
#include "nodeops.h"
#include <stdlib.h>
#include "bitops.h"
#include "SearchNode.h"
#include "memusage.h"

int Node::B_over_8;
UINT8* Node::dbcode;

int main(int argc, char**argv)
{
	if (argc < 3) {
		printf("Usage:\n\nmih <infile> <outfile> [options]\n\n");
		printf("Options:\n");
		printf(" -N <number>          Set the number of binary codes from the beginning of the dataset file to be used\n");
		printf(" -Q <number>          Set the number of query points to use from <infile>, default all\n");
		printf(" -B <number>          Set the number of bits per code, default autodetect\n");
		printf(" -m <number>          Set the number of chunks to use, default 1\n");
		printf(" -K <number>          Set number of nearest neighbors to be retrieved\n");
		printf(" -R <number>          Set the number of codes (in Millions) to use in computing the optimal bit reordering, default OFF (0)\n");
		printf("\n");
		return 0;
	}

	char *infile = argv[1];
	char *outfile = argv[2];

	UINT32 N = 0;
	UINT32 NQ = 0, Q0 = 0, Q1 = 0;
	int B = 0;
	int m = 1;
	UINT32 K = -1;
	size_t R = 0;

	for (int argnum = 3; argnum < argc; argnum++) {
		if (argv[argnum][0] == '-') {
			switch (argv[argnum][1]) {
			case 'B':
				B = atoi(argv[++argnum]);
				break;
			case 'K':
				K = atoi(argv[++argnum]);
				break;
			case 'N':
				N = atoi(argv[++argnum]);
				break;
			case 'Q':
				Q0 = atoi(argv[++argnum]);
				if (++argnum < argc) {
					if (argv[argnum][0] != '-') {
						Q1 = atof(argv[argnum]);
					} else {
						argnum--;
						Q1 = Q0;
						Q0 = 0;
					}
				}
				NQ = Q1-Q0;
				break;
			case 'm':
				m = atoi(argv[++argnum]);
				break;
			case 'R':
				R = atoi(argv[++argnum])*1000000;
				break;
			default:
				printf("Unrecognized Option or Missing Parameter when parsing: %s\n", argv[argnum]);
				return EXIT_FAILURE;
			}
		} else {
			printf("Invalid Argument: %s\n", argv[argnum]);
			return EXIT_FAILURE;
		}
	}

	if (!NQ) {
		printf("-Q is required.\n");
		return EXIT_FAILURE;
	}

	if (B % 8 != 0) {		// in case of B == 0 this should be fine
		printf("Non-multiple of 8 code lengths are not currently supported.\n");
		return EXIT_FAILURE;
	}

	if (R > N) {
		printf("R was greater than N, R will now be equal to N.\n");
		R = N;
	}

	if (K < 1 || K > N) {
		printf("A valid K is not provided.\n");
		return EXIT_FAILURE;
	}

	int B_over_8 = B/8;

	UINT8 *codes_db;
	int dim1codes;
	UINT8 *codes_query;
	int dim1queries;

	clock_t start0, end0;
	time_t start1, end1;

	printf("Loading codes... ");
	fflush(stdout);

	codes_db = (UINT8*)malloc((size_t)N * (B/8) * sizeof(UINT8));

	load_bin_codes(infile, "B", codes_db, &N, &B_over_8);
	printf("Loaded\n");
	fflush(stdout);

	//RootNode root = RootNode(B);

	start1 = time(NULL);
	start0 = clock();
	//	for(UINT32 i=0;i<N;i++)
	//	{
	//		if(i%500000==0)
	//			printf("%d items are load\n",i);
	//		root.insert(codes_db+i*B_over_8,i);
	//	}
	int capacity = 100;
	Node* curr_node;
	UINT8* ccode;
	int tree_lvl = 0;
	int max_tree_lvl = 0;
	//RootNode *root = new RootNode(B,capacity);
	
	
	Node* root_node = new Node(-1,NULL);
	root_node->isleaf = false;
	Node::B_over_8 = B_over_8;
	Node::dbcode = codes_db;
	int max_lvl = 3;
	int expansions=0;
	//UINT8* subnroms = (UINT8*)malloc((B*2-1)*sizeof(UINT8)); //all possible norms of a binary string with B
	printf("Started creating tree\n");
	fflush(stdout);
	for(UINT32 i=0;i<N;i++)
	{
		
		
		if(i==37805)
			printf("i = %d\n",i);
		tree_lvl = 0;
		//printf("%d items are load\n",i);

		
		ccode = codes_db + i*B_over_8;
		
		//root->insert(ccode,curr_node);
		//tree_lvl++;
		curr_node = root_node;
		while(true){
			if(curr_node->isleaf) {

				curr_node->push_to_node(ccode,i);
				if(curr_node->size == capacity && tree_lvl<max_lvl){
					curr_node->expand();
					expansions++;
				}
				break;
			}
			else{
				UINT8 subnorms[pow2(tree_lvl)];
				norm_chunks(subnorms,tree_lvl,ccode,B_over_8);
				//printf("here\n");
				curr_node = curr_node->find_the_child(ccode,subnorms);
			//	if(tree_lvl>=5)
								//	printf("i = %d\n",i);
			}
			//	curr_node == root;
			tree_lvl++;
		}
		if(tree_lvl>max_tree_lvl)
			max_tree_lvl = tree_lvl;
		if(i%500000==0)
		{
			printf("%d items are loaded\n",i);
			fflush(stdout);
		}


	}
	end0 = clock();
	end1 = time(NULL);
	double cput = (double)(end0-start0) / (CLOCKS_PER_SEC) ;
	double wt = (double)(end1-start1) ;
	printf("cput = %f, wt = %f\n",cput,wt);
	printf("Done Loading\n");
	//delete root_node;
	double vm = -1;
	double rss = -1;
	process_mem_usage(&vm, &rss);
    vm  /= double(1024*1024);
    rss /= double(1024*1024);
    printf("VM %.1fgb | RSS %.1fgb\n",vm,rss);
	printf("max tree level = %d\n",max_tree_lvl);
	printf("expansions = %d\n",expansions);
	printf("Done\n");
	codes_query = (UINT8*)malloc((size_t)NQ * (B/8) * sizeof(UINT8));
	UINT8 *codesq, *ccodeq;
	printf("Loading queries\n");
        load_bin_codes(infile, "Q", codes_query, &NQ, &B_over_8, Q0);
        printf("Queries Loaded\n");
	ccodeq = codesq = codes_query;

	UINT64 total_num_comp=0;

	SearchNode* sn = new SearchNode(codes_db, B_over_8,max_tree_lvl, K, root_node,NQ);
	start0 = clock();
	start1 = time(NULL);

	for (UINT32 i=0; i < NQ; i++) {
		sn->setQuery(ccodeq);
		//printf("query is set\n");
		ccodeq += B_over_8;
		//UINT8* set_norm_chunks = new UINT8[(int)pow(2,max_tree_lvl+1)];
		//printf("i = %d\n",i);
		//hammin_radius = 0;
		//visited_parents_list* = NULL; 
		//printf("i = %d\n",i);
		sn->HNN_search();
		//if(sn->num_compare>N)
			//printf("i = %d compare = %d\n",i,sn->num_compare);
		total_num_comp += sn->num_compare;
		//delete norm_chunks;

	}

	end0 = clock();
	end1 = time(NULL);
	cput = (double)(end0-start0) / (CLOCKS_PER_SEC) / NQ;
	wt = (double)(end1-start1)/NQ;
	delete sn;
	printf("cput = %f, wall =%f \n",cput,wt);
	printf("Average number of comparisons = %ld\n",total_num_comp/NQ);
	printf("Search Finished");

}




