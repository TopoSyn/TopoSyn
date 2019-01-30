
//

#ifndef BUCKETALLOCATOR_H
#define BUCKETALLOCATOR_H

enum bcAllocHint
{
	BC_ALLOC_PERM,		
	BC_ALLOC_TEMP		
};

typedef void* (bcAllocFunc)(int size, bcAllocHint hint);
typedef void (bcFreeFunc)(void* ptr);

void bcAllocSetCustom(bcAllocFunc *allocFunc, bcFreeFunc *freeFunc);

void* bcAlloc(int size, bcAllocHint hint);
void bcFree(void* ptr);

#endif
