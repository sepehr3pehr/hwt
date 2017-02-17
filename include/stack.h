
#include "entry.h"
#ifndef STACK_H_
#define STACK_H_

#define nullptr 0

class stack{
public:
	entry* head;
	entry* tail;
	UINT64 size;

	stack();
	void push(Node*);
	Node* pop();
	Node* pop_nofree();
	entry* top();
	~stack();
};

stack::stack() {
	head = nullptr;
	tail = nullptr;
	size = 0;
}

void stack::push(Node* node) {
	if(head==nullptr) {
		head = new entry;
		tail = head;
		head->current = node;
		head->next = nullptr;
		head->prev = nullptr;
		size++;
		return;
	}
	tail->next = new entry;
	tail->next->current = node;
	tail->next->prev = tail;
	tail = tail->next;
	tail->next = nullptr;
	size++;
	return;
}

Node* stack::pop() {
	/*Node* temp = head->current;
	delete head;
	head = head->next;
	size--;*/
	Node* temp = tail->current;
	entry* temp2 = tail->prev;
	delete  tail;
	tail = temp2;
	size--;
	return temp;
}

Node* stack::pop_nofree() {
	entry* temp = tail;
	tail = tail->prev;
	size--;
	return temp->current;
}

entry* stack::top() {
	return tail;
}

stack::~stack() {
	if(size==0)
		return;
	entry* current,*n;
	current = head;
	while(current!=nullptr) {
		n = current->next;
		delete current;
		current = n;
	}

}
#endif
