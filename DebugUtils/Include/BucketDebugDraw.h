//

//

#ifndef BUCKETDEBUGDRAW_H
#define BUCKETDEBUGDRAW_H

#include "DetourNavMesh.h"
#include "VoroDiag.h"
#include "BucketBuilder.h"


void bcDebugDrawIsovist(duDebugDraw* dd,float* verts, int nverts,float h);
void bcDebugDrawBucketGraph(duDebugDraw* dd,bcBucketMesh* bucMesh,BcGraph& g,float h);
void bcDebugDrawBucketGraphColored(duDebugDraw* dd,bcBucketMesh* bucMesh,BcGraph& g,float h,double max);
void bcDebugDrawBucketMesh(struct duDebugDraw* dd, const dtNavMesh& mesh, const dtNavMeshQuery& query, unsigned char flags);


#endif // BUCKETDEBUGDRAW_H