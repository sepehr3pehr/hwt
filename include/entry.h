#include "node.h"
#ifndef entry_H_
#define entry_H_

struct entry {
	Node* current;
	entry* next;
	entry* prev;
};

#endif
