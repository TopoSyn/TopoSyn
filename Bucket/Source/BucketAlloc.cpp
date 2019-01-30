

#include <stdlib.h>
#include "BucketAlloc.h"


static void *bcAllocDefault(int size, bcAllocHint)
{
	return malloc(size);
}

static void bcFreeDefault(void *ptr)
{
	free(ptr);
}

static bcAllocFunc* sAllocFunc = bcAllocDefault;
static bcFreeFunc* sFreeFunc =  bcFreeDefault;

void bcAllocSetCustom(bcAllocFunc *allocFunc, bcFreeFunc *freeFunc)
{
	sAllocFunc = allocFunc ? allocFunc : bcAllocDefault;
	sFreeFunc = freeFunc ? freeFunc : bcFreeDefault;
}

void* bcAlloc(int size, bcAllocHint hint)
{
	return sAllocFunc(size, hint);
}

void bcFree(void* ptr)
{
	if (ptr)
		sFreeFunc(ptr);
}
