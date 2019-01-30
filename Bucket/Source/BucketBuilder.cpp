

#include "BucketBuilder.h"
#include "DetourCommon.h"
#include "DetourAssert.h"

#include "Recast.h"
#include "Sample.h"

#include "Cluster.h"


#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/betweenness_centrality.hpp>
#include <boost/foreach.hpp> 

// Create constant for PI.
#define PI 3.14159265
// Uncomment this to dump all the requests in stdout.
#define DUMP_REQS


static float distPtLine2d(const float* pt, const float* p, const float* q)
{
	float pqx = q[0] - p[0];
	float pqz = q[2] - p[2];
	float dx = pt[0] - p[0];
	float dz = pt[2] - p[2];
	float d = pqx*pqx + pqz*pqz;
	float t = pqx*dx + pqz*dz;
	if (d != 0) t /= d;
	dx = p[0] + t*pqx - pt[0];
	dz = p[2] + t*pqz - pt[2];
	return dx*dx + dz*dz;
}

static int calcAreaOfIsovist(const int* verts, const int nverts)
{
	int area = 0;
	for (int i = 0, j = nverts-1; i < nverts; j=i++)
	{
		const int* vi = &verts[i*4];
		const int* vj = &verts[j*4];
		area += vi[0] * vj[2] - vj[0] * vi[2];
	}
	return (area+1) / 2;
}
 
static float CalcAngle2(float* p1, float* p2, float* p3, float* p4)
{

float a = p1[0] - p2[0];
float b = p1[2] - p2[2];
float c = p3[0] - p4[0];
float d = p3[2] - p4[2];
//
float cos_angle , angle;
float mag_v1 = sqrt(a*a + b*b);
float mag_v2 = sqrt(c*c + d*d);
//
cos_angle = (a*c + b*d) / (mag_v1 * mag_v2);
angle = acos(cos_angle);
angle = angle * 180.0 / 3.14; 
//
return angle;
}
static float* segmentMidpoint(float* a,float* b)
{
	float c[3];
	c[0]=(a[0]+b[0])/2;
	c[1]=a[1];
	c[2]=(a[2]+b[2])/2;
	return c;
}

bcBucketMesh* bcAllocBucketMesh(BcGraph& bcg,BcGraph& vbcg)
{
	void* mem = bcAlloc(sizeof(bcBucketMesh), BC_ALLOC_PERM);
	if (!mem) return 0;
	return new(mem) bcBucketMesh(bcg,vbcg);
}

void bcFreeBucketMesh(bcBucketMesh* bucmesh)
{
	if (!bucmesh) return;
	bucmesh->~bcBucketMesh();
	dtFree(bucmesh);
}

void bcFreeBucketSet(bcBucketSet* bset)
{
	if (!bset) return;
	
	bcFree(bset);
}

bcBucketMesh::bcBucketMesh(BcGraph& bcg,BcGraph& vbcg):m_bucketGraph(bcg),m_visibilityBucketGraph(vbcg)
{}

bcBucketMesh::~bcBucketMesh()
{
	(m_bucketGraph).clear();
	(m_visibilityBucketGraph).clear();
}

 char bcBucketMesh::bcSegSegIntersect2d(const float* a,const float* b,const float* c,const float* d,double& s, double& t)
 {
	 char code='?';
	 double num ;
	 double denom ;

	 //******************************************************************
	 denom=a[0]*(d[2]-c[2])+
	       b[0]*(c[2]-d[2])+
	       d[0]*(b[2]-a[2])+
	       c[0]*(a[2]-b[2]);
	
	 if (denom==0.0)
	 {
		 code='p';
		 return code;
	 }

	 num=a[0]*(d[2]-c[2])+
	     c[0]*(a[2]-d[2])+
	     d[0]*(c[2]-a[2]);

	 if ((num==0)||(num==denom))
	 {
	   code='e';
	 }

	 s=num/denom;

	 num=-(a[0]*(c[2]-b[2])+
	       b[0]*(a[2]-c[2])+
	       c[0]*(b[2]-a[2]));

	 if ((num==0)||(num==denom))
	 {
	   code='e';
	 }

	 t=num/denom;

	 if ((0.0<=t)&&(t<=1.0))
	 {
		 code='1';
	 }
	 else 
	 {
		 code='0';
	 }
#ifdef DUMP_REQS
			
#endif
	 return code;
 }

 char bcBucketMesh::bcVerSegSegIntersect2d(const float* a,const float* b,const float* c,const float* d,float& s, float& t)
 {
	 char code='?';

	 //**************************AB vertical*******************************
	 if (a[0]-b[0]==0.0)
		 if(c[0]!=d[0])
		{
			s=(c[2]-d[2])/(c[0]-d[0]);
			t=c[2]-s*c[0];
			code='v';
			return code;
		 }
		 else
		 {
			 code='p';
			 return code;
		 }

	 //**************************AB horizontal*******************************
	 if (a[2]-b[2]==0.0)
		if ((c[2]!=d[2])&&(c[0]!=d[0]))
		{
			s=(c[2]-d[2])/(c[0]-d[0]);
			t=c[2]-s*c[0];
			code='h';
			return code;
	    }
		else
		{
			code='p';
			return code;
	    }

 }

 bool bcBucketMesh::bcIntersectLinePoly2D(const float* p0, const float* p1,
							  const float* verts, int nverts,
							  double& smin, double& smax,
							  int& segMin, int& segMax)
{
	segMin = -1;
	segMax = -1;
	double t,s;

	for (int i = 0, j = nverts-1; i < nverts; j=i++)
	{
		char code =bcSegSegIntersect2d(p0,p1,&verts[i*3], &verts[j*3],s,t);
#ifdef DUMP_REQS
			
#endif
		
		if (code != '1')
		{
			continue;
		}
		else 
		{
			smin=s;
		    smax=s;
			segMin = j;
			segMax = j;
			for (int k = i, l = j; k < nverts; l=k++)
	        {
				char code =bcSegSegIntersect2d(p0,p1,&verts[k*3], &verts[l*3],s,t);
				if (code != '1')
				{
					continue;
				}
				else 
				{
					if (s < smin)
					{
						smin = s;
						segMin = l;
				
						if (smin > smax)
							return false;
					}
					else if (s > smax)
					{
						smax = s;
						segMax = l;
				
				        if (smax < smin)
							return false;
					}
				}
			}
			
			return true;
		}
	}

	return false;
}

 bool bcBucketMesh::bcIntersectVerLinePoly2D(const float* p0, const float* p1,
							  const float* verts, int nverts,
							  float& ymin, float& ymax,
							  int& segMin, int& segMax)
{
	segMin = -1;
	segMax = -1;
	float a,b,t,y;
	float min,max;
	
	for (int i = 0, j = nverts-1; i < nverts; j=i++)
	{
		char code =bcVerSegSegIntersect2d(p0,p1,&verts[i*3],&verts[j*3],a,b);
		const float* vp1=&verts[i*3];
		const float* vp2=&verts[j*3];

		if (code != 'v')
		{
			continue;
		}
		else 
		{
			ymin=a*p0[0]+b;
			min=vp1[2];
			max=vp2[2];
			if (min>max)
			{
				t=min;
				min>max;
				max=t;
			}
			if ((ymin<min)||(ymin>max))
				continue;

		    ymax=ymin;
			segMin = j;
			segMax = j;
			for (int k = i, l = j; k < nverts; l=k++)
	        {
				char code =bcVerSegSegIntersect2d(p0,p1,&verts[k*3], &verts[l*3],a,b);
				const float* vp1=&verts[k*3];
				const float* vp2=&verts[l*3];
				if (code != 'v')
				{
					continue;
				}
				else 
				{
					y=a*p0[0]+b;
					min=vp1[2];
					max=vp2[2];
					if (min>max)
					{
						t=min;
						min>max;
						max=t;
					}
					if ((y<min)||(y>max))
						continue;

					if (y < ymin)
					{
						ymin = y;
						segMin = l;
				
						if (ymin > ymax)
							return false;
					}
					else if (y > ymax)
					{
						ymax = y;
						segMax = l;
				
				        if (ymax < ymin)
							return false;
					}
				}
			}
#ifdef DUMP_REQS
			//printf("Isovist ridge2  %f %f  \n", ymax,ymin); 
#endif
			return true;
		}
	}
#ifdef DUMP_REQS  
			//printf("Isovist ridge3  %f %f  \n", ymax,ymin); 
#endif
	return false;
}

void bcBucketMesh::addEdge(MaGraph& graph,const float* v0,const float* v1)
{
	ma_vertex_t vg,vg0,vg1;
	std::pair<ma_vertex_iterator_t, ma_vertex_iterator_t> vp;
	bool exist= false;
	for (vp = vertices(graph); vp.first != vp.second; ++vp.first)
		{
			if ((graph[*vp.first].x==v0[0])&&(graph[*vp.first].y==v0[2]))
			{
				vg=*vp.first;
				exist= true;
				break;
			}
		}
	if (exist)
	{
		vg0=vg;
	}
	else
	{
		vg0 = boost::add_vertex(VertexProperties(v0[0], v0[2]),graph);	
	}
	
	exist= false;
	
	for (vp = vertices(graph); vp.first != vp.second; ++vp.first)
	{
		if ((graph[*vp.first].x==v1[0])&&(graph[*vp.first].y==v1[2]))
		{
			vg=*vp.first;
			exist= true;
			break;
		}
	}
	if (exist)
	{
		vg1=vg;
	}
	else
	{
		vg1 = boost::add_vertex(VertexProperties(v1[0], v1[2]),graph);
	}

	boost::add_edge(vg0, vg1, graph);
	boost::add_edge(vg1, vg0, graph);
}

int bcBucketMesh::getNeighbours(dtNavMesh* mesh,dtPolyRef queryPointRef,dtPolyRef* neighRef)
{
    
	int polyIndex=0;

	const dtMeshTile* queryPolyTile = 0;
	const dtPoly* queryPoly = 0;
	if (dtStatusFailed(mesh->getTileAndPolyByRef(queryPointRef, &queryPolyTile, &queryPoly)))
		return 0;
	
	neighRef[0]=queryPointRef;
	int j=0;
	
	for (unsigned int i = queryPoly->firstLink; i != DT_NULL_LINK; i = queryPolyTile->links[i].next)
	{
		dtPolyRef neighbourRef = queryPolyTile->links[i].ref;
		
		if (neighbourRef)
		{
			neighRef[j]=neighbourRef;
			j++;
		}
	}
	return j;
}

int bcBucketMesh::getWalls(dtNavMesh* mesh,dtPolyRef queryPointRef,int* e)
{
	int polyIndex=0;

	const dtMeshTile* queryPolyTile = 0;
	const dtPoly* queryPoly = 0;
	if (dtStatusFailed(mesh->getTileAndPolyByRef(queryPointRef, &queryPolyTile, &queryPoly)))
		return 0;

	int k=0;
	for (int i = 0, j = (int)queryPoly->vertCount-1; i < (int)queryPoly->vertCount; j = i++)
		{
			// Skip non-solid edges.
			if (queryPoly->neis[j] & DT_EXT_LINK)
			{
				// Tile border.
				bool solid = true;
				for (unsigned int k = queryPoly->firstLink; k != DT_NULL_LINK; k = queryPolyTile->links[k].next)
				{
					const dtLink* link = &queryPolyTile->links[k];
					if (link->edge == j)
					{
						if (link->ref != 0)
							solid = false;
						break;
					}
				}
				if (!solid) continue;
			}
			else if (queryPoly->neis[j])
			{
				// Internal edge
				continue;
			}
			e[k] =j;
		    k++;
		}
	return k;
}

bool bcBucketMesh::getABPointsRef(dtNavMesh* mesh,dtPolyRef from, dtPolyRef to, float* left, float* right,
								unsigned char& fromType, unsigned char& toType) 
{
	const dtMeshTile* fromTile = 0;
	const dtPoly* fromPoly = 0;
	if (dtStatusFailed(mesh->getTileAndPolyByRef(from, &fromTile, &fromPoly)))
		return DT_FAILURE | DT_INVALID_PARAM;
	fromType = fromPoly->getType();

	const dtMeshTile* toTile = 0;
	const dtPoly* toPoly = 0;
	if (dtStatusFailed(mesh->getTileAndPolyByRef(to, &toTile, &toPoly)))
		return DT_FAILURE | DT_INVALID_PARAM;
	toType = toPoly->getType();

	return getABPoints(from, fromPoly, fromTile,
						   to, toPoly, toTile,
						   left, right);
}

bool bcBucketMesh::getABPoints(dtPolyRef from, const dtPoly* fromPoly, const dtMeshTile* fromTile,
								dtPolyRef to, const dtPoly* toPoly, const dtMeshTile* toTile,
								float* left, float* right)
{
	// Find the link that points to the 'to' polygon.
	const dtLink* link = 0;
	for (unsigned int i = fromPoly->firstLink; i != DT_NULL_LINK; i = fromTile->links[i].next)
	{
		if (fromTile->links[i].ref == to)
		{
			link = &fromTile->links[i];
			break;
		}
	}
	if (!link)
		return false;
	
	// Handle off-mesh connections.
	if (fromPoly->type == DT_POLYTYPE_OFFMESH_CONNECTION)
	{
		// Find link that points to first vertex.
		for (unsigned int i = fromPoly->firstLink; i != DT_NULL_LINK; i = fromTile->links[i].next)
		{
			if (fromTile->links[i].ref == to)
			{
				const int v = fromTile->links[i].edge;
				dtVcopy(left, &fromTile->verts[fromPoly->verts[v]*3]);
				dtVcopy(right, &fromTile->verts[fromPoly->verts[v]*3]);
				return true;
			}
		}
		return false;
	}

	if (toPoly->type == DT_POLYTYPE_OFFMESH_CONNECTION)
	{
		for (unsigned int i = toPoly->firstLink; i != DT_NULL_LINK; i = toTile->links[i].next)
		{
			if (toTile->links[i].ref == from)
			{
				const int v = toTile->links[i].edge;
				dtVcopy(left, &toTile->verts[toPoly->verts[v]*3]);
				dtVcopy(right, &toTile->verts[toPoly->verts[v]*3]);
				return true;
			}
		}
		return false;
	}
		
	// Find portal vertices.
	const int v0 = fromPoly->verts[link->edge];
	const int v1 = fromPoly->verts[(link->edge+1) % (int)fromPoly->vertCount];
	dtVcopy(left, &fromTile->verts[v0*3]);
	dtVcopy(right, &fromTile->verts[v1*3]);
	
	// If the link is at tile boundary, dtClamp the vertices to
	// the link width.
	if (link->side == 0 || link->side == 4)
	{
		// Unpack portal limits.
		const float smin = dtMin(left[2],right[2]);
		const float smax = dtMax(left[2],right[2]);
		const float s = (smax-smin) / 255.0f;
		const float lmin = smin + link->bmin*s;
		const float lmax = smin + link->bmax*s;
		left[2] = dtMax(left[2],lmin);
		left[2] = dtMin(left[2],lmax);
		right[2] = dtMax(right[2],lmin);
		right[2] = dtMin(right[2],lmax);
	}
	else if (link->side == 2 || link->side == 6)
	{
		// Unpack portal limits.
		const float smin = dtMin(left[0],right[0]);
		const float smax = dtMax(left[0],right[0]);
		const float s = (smax-smin) / 255.0f;
		const float lmin = smin + link->bmin*s;
		const float lmax = smin + link->bmax*s;
		left[0] = dtMax(left[0],lmin);
		left[0] = dtMin(left[0],lmax);
		right[0] = dtMax(right[0],lmin);
		right[0] = dtMin(right[0],lmax);
	}
	
	return true;
}

void bcBucketMesh::bcView(dtNavMesh* mesh,MaGraph& vg,const float* portalApex,dtPolyRef poly,dtPolyRef parentpoly,float* A,float* B,float* Right,float* Left)
{
	dtPolyRef neighpoly[6]; 
	float portalRight[3], portalLeft[3];
	float LeftPrime[3];
	float RightPrime[3];
	float newLeft[3];
	float newRight[3];
	unsigned char queryType, neighType;

	
	float polyVerts[12];
	const float* v=0;
	int numBorVerts;
	
	double s,t;
	int polyIdx=0;

	const dtMeshTile* currentPolyTile = 0;
	const dtPoly* currentPoly = 0;
	dtStatusFailed(mesh->getTileAndPolyByRef(poly, &currentPolyTile, &currentPoly));
		
	
//******************************************border edges***************************************
	int walls[6];
	int ne=getWalls(mesh,poly,walls);
	for (int j = 0;j < ne; ++j)
	{
		int w=walls[j];
		const int v0 = currentPoly->verts[w];
		dtVcopy(portalLeft,&currentPolyTile->verts[v0*3]);
		const int v1 = currentPoly->verts[(w+1) % (int)currentPoly->vertCount];
		dtVcopy(portalRight, &currentPolyTile->verts[v1*3]);
		
		if ((dtTriArea2D(portalApex, Left, portalRight) > 0.0f)||(dtTriArea2D(portalApex, Right, portalLeft) < 0.0f))
		{
			continue;
		}
		if (dtTriArea2D(portalApex, Left, portalLeft) <= 0.0f)
		{
			dtVcopy(LeftPrime,portalLeft);
		}
		else
		{
			 char c=bcSegSegIntersect2d(portalApex,Left,portalLeft,portalRight,s,t);
			 if (c=='1')
			 {
				LeftPrime[0] = portalApex[0] + (Left[0] - portalApex[0]) * s;
				LeftPrime[1] = portalApex[1] + (Left[1] - portalApex[1]) * s;
				LeftPrime[2] = portalApex[2] + (Left[2] - portalApex[2]) * s;
			 }
			 if (c=='v')
			 {
				 dtVcopy(LeftPrime,portalLeft);
			 }
		}

		if (dtTriArea2D(portalApex,Right,portalRight) >= 0.0f)
		{
			dtVcopy(RightPrime,portalRight);
		}
		else
		{		
			char c=bcSegSegIntersect2d(portalApex,Right,portalLeft,portalRight,s,t);
			if (c=='1')
			{
				RightPrime[0] = portalApex[0] + (Right[0] - portalApex[0]) * s;
				RightPrime[1] = portalApex[1] + (Right[1] - portalApex[1]) * s;
				RightPrime[2] = portalApex[2] + (Right[2] - portalApex[2]) * s;
			}
			if (c=='v')
			{
				 dtVcopy(RightPrime,portalRight);
			}
		}

		addEdge(vg,LeftPrime,portalRight);  
		addEdge(vg,RightPrime,portalLeft); 
	
	}

//***********************************Internal Poly********************************************

	int neighSize=getNeighbours(mesh,poly,neighpoly);
	int c=0;
	for (int i = 0; i < neighSize; ++i)
	{
			if (neighpoly[i]==parentpoly)
			{
				c=c+1;
				continue;
			}

			getABPointsRef(mesh,poly,neighpoly[i],portalRight,portalLeft,queryType,neighType);
			if ((dtTriArea2D(portalApex, Left, portalRight) > 0.0f)||(dtTriArea2D(portalApex, Right, portalLeft) < 0.0f))
			{
				c=c+1;
				continue;
			}

			dtVcopy(newLeft,Left);
			dtVcopy(newRight,Right);

			// portal Left vertex.
			if (dtTriArea2D(portalApex,Left,portalLeft) <= 0.0f)
			{
				dtVcopy(newLeft,portalLeft);
				addEdge(vg,newLeft,Left);
			}

			// Right vertex.
			if (dtTriArea2D(portalApex,Right,portalRight) >= 0.0f)
			{
				dtVcopy(newRight, portalRight);
				addEdge(vg,newRight,Right);
						
			}

			bcView(mesh,vg,portalApex,neighpoly[i],poly,portalLeft,portalRight,newRight,newLeft); 
			
	}
	if (c==neighSize)
	{
		addEdge(vg,Left,Right);
	}
}

int bcBucketMesh::bcCalcIsovist(dtNavMeshQuery* navQuery,dtNavMesh* mesh,dtPolyRef queryPointRef,const float* queryPointPos,MaGraph& vg)
{
	float closestportalApex[3];
	navQuery->closestPointOnPoly(queryPointRef, queryPointPos, closestportalApex,0);
	
	float reflexLeft[3], reflexRight[3];
	float portalA[3], portalB[3];
	int neighSize;
	dtPolyRef poly[6];
	unsigned char queryType, neighType;
	float* borderVerts[12];
	int numBorVerts;
	bool exist;

	int polyIdx=0;

	const dtMeshTile* queryPolyTile = 0;
	const dtPoly* queryPoly = 0;
	if (dtStatusFailed(mesh->getTileAndPolyByRef(queryPointRef, &queryPolyTile, &queryPoly)))
		return 0;

	for (int j = 0, nj = (int)queryPoly->vertCount; j < nj; ++j)
		{
			const float* v = &queryPolyTile->verts[queryPoly->verts[j]*3];
			addEdge(vg,queryPointPos,v);
	    }
	
	int walls[6];
	int ne=getWalls(mesh,queryPointRef,walls);
	for (int j = 0;j < ne; ++j)
	{
		int w=walls[j];
		const int v0 = queryPoly->verts[w];
		const int v1 = queryPoly->verts[(w+1) % (int)queryPoly->vertCount];
		addEdge(vg,&queryPolyTile->verts[v0*3],&queryPolyTile->verts[v1*3]);          
	}


	neighSize=getNeighbours(mesh,queryPointRef,poly);
#ifdef DUMP_REQS

#endif 

	for (int i = 0; i < neighSize; ++i)
		{
			getABPointsRef(mesh,queryPointRef,poly[i], portalA,portalB, queryType,neighType);
			
			dtVcopy(reflexRight,portalA);
			dtVcopy(reflexLeft,portalB);
			
			bcView(mesh,vg,closestportalApex,poly[i],queryPointRef,portalA,portalB,reflexRight,reflexLeft);
	    }			
	int numV=boost::num_vertices(vg);		
	
	return numV;

}

int bcBucketMesh::bcRaycast1(dtNavMesh* mesh,dtPolyRef centerRef, const float* startPos, const float* endPos, const dtQueryFilter* filter,double& sa) 
{
	sa = 0.0;
	
	if (!centerRef || !mesh->isValidPolyRef(centerRef))
		return 0;
	
	dtPolyRef curRef = centerRef;
	float verts[DT_VERTS_PER_POLYGON*3];	
	int n = 0;
	double smin, smax;
	int segMin, segMax;

	const dtMeshTile* curTile = 0;
	const dtPoly* curPoly = 0;
	if (dtStatusFailed(mesh->getTileAndPolyByRef(curRef, &curTile, &curPoly)))
		return DT_FAILURE | DT_INVALID_PARAM;

	//**********************Max hit*************************************************
	// Collect vertices.
	int nv = 0;
	for (int i = 0; i < (int)curPoly->vertCount; ++i)
	{
		dtVcopy(&verts[nv*3], &curTile->verts[curPoly->verts[i]*3]);
		nv++;
	}		

	if (!bcIntersectLinePoly2D(startPos, endPos,verts, nv,smin, smax, segMin, segMax))
		return 0;

	do
    {
		dtPolyRef nextRef = 0;
		for (unsigned int i = curPoly->firstLink; i != DT_NULL_LINK; i = curTile->links[i].next)
		{
			const dtLink* link = &curTile->links[i];
			
			// Find link which contains this edge.
			if ((int)link->edge != segMax)
				continue;
			
			const dtMeshTile* nextTile = 0;
			const dtPoly* nextPoly = 0;
			if (dtStatusFailed(mesh->getTileAndPolyByRef(link->ref, &nextTile, &nextPoly)))
				return DT_FAILURE | DT_INVALID_PARAM;

			// Skip off-mesh connections.
			if (nextPoly->type == DT_POLYTYPE_OFFMESH_CONNECTION)
				continue;
				
			// Skip links based on filter.
			/*if (!passFilter(filter, nextPoly->flags))
				continue;*/
		
			// If the link is internal, just return the ref.
			if (link->side == 0xff)
			{
				nextRef = link->ref;
				break;
			}

			// If the link is at tile boundary,
			const int v0 = curPoly->verts[link->edge];
			const int v1 = curPoly->verts[(link->edge+1) % curPoly->vertCount];
			const float* left = &curTile->verts[v0*3];
			const float* right = &curTile->verts[v1*3];
			
			// Check that the intersection lies inside the link portal.
			if (link->side == 0 || link->side == 4)
			{
				// Calculate link size.
				const float ssmin = dtMin(left[2],right[2]);
				const float ssmax = dtMax(left[2],right[2]);
				const float ss = (ssmax-ssmin) / 255.0f;
				const float lmin = ssmin + link->bmin*ss;
				const float lmax = ssmin + link->bmax*ss;
				// Find Z intersection.
				float z = startPos[2] + (endPos[2]-startPos[2])*smax;
				if (z >= lmin && z <= lmax)
				{
					nextRef = link->ref;
					break;
				}
			}
			else if (link->side == 2 || link->side == 6)
			{
				// Calculate link size.
				const float ssmin = dtMin(left[0],right[0]);
				const float ssmax = dtMax(left[0],right[0]);
				const float ss = (ssmax-ssmin) / 255.0f;
				const float lmin = ssmin + link->bmin*ss;
				const float lmax = ssmin + link->bmax*ss;
				// Find X intersection.
				float x = startPos[0] + (endPos[0]-startPos[0])*smax;
				if (x >= lmin && x <= lmax)
				{
					nextRef = link->ref;
					break;
				}
			}
		
		}
		n++;

		if (!nextRef)
		{
			sa=smax;
			return n;
		}
		
		const dtMeshTile* nextTile = 0;
		const dtPoly* nextPoly = 0;
		if (dtStatusFailed(mesh->getTileAndPolyByRef(nextRef, &nextTile, &nextPoly)))
			return DT_FAILURE | DT_INVALID_PARAM;

		// Collect vertices.
		nv = 0;
		for (int i = 0; i < (int)nextPoly->vertCount; ++i)
		{
			dtVcopy(&verts[nv*3], &nextTile->verts[nextPoly->verts[i]*3]);
			nv++;
		}		

		bcIntersectLinePoly2D(startPos, endPos,verts, nv,smin, smax, segMin, segMax);
		
		curRef=nextRef;
		curTile = nextTile;
		curPoly = nextPoly;
		
	} while (curRef);

	sa=smax;
	return n;
}

int bcBucketMesh::bcRaycast2(dtNavMesh* mesh,dtPolyRef centerRef, const float* startPos, const float* endPos, const dtQueryFilter* filter, double& si) 
{
	si = 0.0;
	
	if (!centerRef || !mesh->isValidPolyRef(centerRef))
		return 0;
	
	dtPolyRef curRef = centerRef;
	float verts[DT_VERTS_PER_POLYGON*3];	
	int n = 0;
	double smin, smax;
	int segMin, segMax;

	const dtMeshTile* curTile = 0;
	const dtPoly* curPoly = 0;
	if (dtStatusFailed(mesh->getTileAndPolyByRef(curRef, &curTile, &curPoly)))
		return DT_FAILURE | DT_INVALID_PARAM;

	//**********************Min hit*************************************************
	// Collect vertices.
	int nv = 0;
	for (int i = 0; i < (int)curPoly->vertCount; ++i)
	{
		dtVcopy(&verts[nv*3], &curTile->verts[curPoly->verts[i]*3]);
		nv++;
	}		

	if (!bcIntersectLinePoly2D(startPos, endPos,verts, nv,smin, smax, segMin, segMax))
		return 0;

	do
    {

		dtPolyRef nextRef = 0;
		for (unsigned int i = curPoly->firstLink; i != DT_NULL_LINK; i = curTile->links[i].next)
		{
			const dtLink* link = &curTile->links[i];
			
			// Find link which contains this edge.
			if ((int)link->edge != segMin)
				continue;
				
			
			const dtMeshTile* nextTile = 0;
			const dtPoly* nextPoly = 0;
			if (dtStatusFailed(mesh->getTileAndPolyByRef(link->ref, &nextTile, &nextPoly)))
				return DT_FAILURE | DT_INVALID_PARAM;
			
			// Skip off-mesh connections.
			if (nextPoly->type == DT_POLYTYPE_OFFMESH_CONNECTION)
				continue;
		
			// If the link is internal, just return the ref.
			if (link->side == 0xff)
			{
				nextRef = link->ref;
				break;
			}

			// If the link is at tile boundary,
			const int v0 = curPoly->verts[link->edge];
			const int v1 = curPoly->verts[(link->edge+1) % curPoly->vertCount];
			const float* left = &curTile->verts[v0*3];
			const float* right = &curTile->verts[v1*3];
			
			// Check that the intersection lies inside the link portal.
			if (link->side == 0 || link->side == 4)
			{
				// Calculate link size.
				const float ssmin = dtMin(left[2],right[2]);
				const float ssmax = dtMax(left[2],right[2]);
				const float ss = (ssmax-ssmin) / 255.0f;
				const float lmin = ssmin + link->bmin*ss;
				const float lmax = ssmin + link->bmax*ss;
				// Find Z intersection.
				float z = startPos[2] + (endPos[2]-startPos[2])*smin;
				if (z >= lmin && z <= lmax)
				{
					nextRef = link->ref;
					break;
				}
			}
			else if (link->side == 2 || link->side == 6)
			{
				// Calculate link size.
				const float ssmin = dtMin(left[0],right[0]);
				const float ssmax = dtMax(left[0],right[0]);
				const float ss = (ssmax-ssmin) / 255.0f;
				const float lmin = ssmin + link->bmin*ss;
				const float lmax = ssmin + link->bmax*ss;
				// Find X intersection.
				float x = startPos[0] + (endPos[0]-startPos[0])*smin;
				if (x >= lmin && x <= lmax)
				{
					nextRef = link->ref;
					break;
				}
			}
		
		}
		n++;

		if (!nextRef)
		{
			si=smin;
			return n;
		}

		const dtMeshTile* nextTile = 0;
		const dtPoly* nextPoly = 0;
		if (dtStatusFailed(mesh->getTileAndPolyByRef(nextRef, &nextTile, &nextPoly)))
			return DT_FAILURE | DT_INVALID_PARAM;

		// Collect vertices.
		nv = 0;
		for (int i = 0; i < (int)nextPoly->vertCount; ++i)
		{
			dtVcopy(&verts[nv*3], &nextTile->verts[nextPoly->verts[i]*3]);
			nv++;
		}		

		bcIntersectLinePoly2D(startPos, endPos,verts, nv,smin, smax, segMin, segMax);
		
		curRef=nextRef;
		curTile = nextTile;
		curPoly = nextPoly;
		
	} while (curRef);
	
	si=smin;
	return n;
}

void bcBucketMesh::bcHitPos(const float* spos,float* epos,float* hitPosMax, float* hitPosMin,double smin, double smax)
{
	hitPosMax[0] = spos[0] + (epos[0] - spos[0]) * (float)smax;
	hitPosMax[1] = spos[1] + (epos[1] - spos[1]) * (float)smax;
	hitPosMax[2] = spos[2] + (epos[2] - spos[2]) * (float)smax;

	hitPosMin[0] = spos[0] + (epos[0] - spos[0]) * (float)smin;
	hitPosMin[1] = spos[1] + (epos[1] - spos[1]) * (float)smin;
	hitPosMin[2] = spos[2] + (epos[2] - spos[2]) * (float)smin;
}

float bcBucketMesh::bcDistancePtPtSqr2D(float* A, float* B)
{
	float dx = (B[0]-A[0]);
	float dz = (B[2]-A[2]);
	return dx*dx + dz*dz;
}

void bcBucketMesh::bcBuildIsovistPoly(MaGraph isoG,std::vector<polygon_t>& resultVect)
{
	std::map<vertex_t,size_t>::iterator itc;
	std::pair<ma_vertex_iterator_t, ma_vertex_iterator_t> vp;

	// extraction of connected components*****************************************
	compMap_t compMap;
	boost::associative_property_map<compMap_t> componentMap(compMap); 
	int num = boost::connected_components(isoG, componentMap); 
	
	// creation of polys***********************************************************
	float fx[2];
	pointsVect_t polyPoints;
	polygon_t poly;
	std::vector<polygon_t> polysVect; 
	int k=0;
	for (int i=0;i<num;i++)
	{
	    for(itc=compMap.begin(); itc!=compMap.end(); ++itc)
		{
			if (itc->first!=i)
				continue;
			for (vp = vertices(isoG); vp.first != vp.second; ++vp.first)
				{
					if (itc->first==*vp.first)
					{
						fx[0] = isoG[*vp.first].x;
						fx[1] = isoG[*vp.first].y;
						point_t point(fx[0],fx[1]);
						polyPoints.push_back(point);
						k++;
					}
				}
		}
		boost::geometry::append(poly,polyPoints);
        boost::geometry::correct(poly);
		polysVect.push_back(poly);
		polyPoints.clear();
	}
}

static void getPolyCenter(dtNavMesh* navMesh, dtPolyRef ref, float* center)
{
	center[0] = 0;
	center[1] = 0;
	center[2] = 0;
	
	const dtMeshTile* tile = 0;
	const dtPoly* poly = 0;
	dtStatus status = navMesh->getTileAndPolyByRef(ref, &tile, &poly);
	if (dtStatusFailed(status))
		return;
		
	for (int i = 0; i < (int)poly->vertCount; ++i)
	{
		const float* v = &tile->verts[poly->verts[i]*3];
		center[0] += v[0];
		center[1] += v[1];
		center[2] += v[2];
	}
	const float s = 1.0f / poly->vertCount;
	center[0] *= s;
	center[1] *= s;
	center[2] *= s;
}

static void calcCentroid(dtNavMesh* navMesh, dtPolyRef ref, float* center)
{
	const dtMeshTile* polyTile = 0;
	const dtPoly* poly = 0;

	if (dtStatusFailed(navMesh->getTileAndPolyByRef(ref, &polyTile, &poly)))
			return;

	polygon_t p;
	std::vector<point_t> vectorPoints;
    int i=0;
	
	for (int j = 0, nj = (int)poly->vertCount; j < nj; ++j)
	{
		float* v0 = &polyTile->verts[poly->verts[j]*3];
		double t[2];
		t[0]=double(v0[0]); t[1]=double(v0[2]);
		center[1] = v0[1];
		point_t polyPoint(t[0],t[1] );
		vectorPoints.push_back(polyPoint);
	}
	
	boost::geometry::append(p, vectorPoints);
    boost::geometry::correct(p);
	vectorPoints.clear();

	point_t c;
    boost::geometry::centroid(p, c);
	center[0] = float(c.x());
	center[2] = float(c.y());
}

bool bcBucketMesh::bcBuildBucket(dtNavMeshQuery* navQuery,const dtMeshTile* tile,dtNavMesh* mesh,bcBucketSet& bset)
{
	
	dtQueryFilter filter;
	dtPolyRef centreRef;

	dtPolyRef neighbour[6];

	float lvert[6*3];
	float rvert[6*3];

	int m=0;

	filter.setIncludeFlags(SAMPLE_POLYFLAGS_ALL ^ SAMPLE_POLYFLAGS_DISABLED);
	filter.setExcludeFlags(0);

	std::vector<dtPoly*> poly_vector;
	std::vector<dtPolyRef> polyRef_vector;
	const dtMeshTile* polyTile = tile;
	const dtPolyRef base = mesh->getPolyRefBase(tile);

	for (int i = 0; i < tile->header->polyCount; ++i)
	{
		dtPoly* poly = &tile->polys[i];

		int nb=0;
		const int nv = poly->vertCount;
		for (int j = 0; j < nv; ++j)
		{
			// Skip non-portal edges.
			if ((poly->neis[j]) == 0)
				continue;
			nb++;
		}

		if ((nb>2)||(nb==1))
		{
			const dtPolyRef polyRef = base | (dtPolyRef)i;
			polyRef_vector.push_back(polyRef);
			poly_vector.push_back(poly);
		}
	}

	int maxBuckets=polyRef_vector.size();
	bset.buckets = (bcBucket*)bcAlloc(sizeof(bcBucket)*maxBuckets, BC_ALLOC_PERM);
	if (!bset.buckets)
		return false;
	bset.nbuckets = 0;
	unsigned short int Id=0;

	for (std::size_t i = 0; i < maxBuckets; ++i) 
	    {
			dtPoly* poly=poly_vector[i];
			dtPolyRef polyRef=polyRef_vector[i];
			
			//*** Bucket Portal Neighbours******//
			dtStatus status = navQuery->findBucketPortalNeighbours(polyRef,neighbour,lvert,rvert,&filter,&m);
			
			if (dtStatusSucceed(status))
			{
				bcBucket* bucket = &bset.buckets[bset.nbuckets++];
			
				//*** Bucket index******//
				bucket->bcBucRef=Id;
				Id++;
			
				//*** Bucket center point******//
				calcCentroid(mesh,polyRef, bucket->centre);
				//getPolyCenter(mesh,polyRef, bucket->centre);
			
				//*** Bucket Polygons******//
				bucket->npolys=1;
				bucket->bucketPoly = (dtPolyRef*)bcAlloc(sizeof(dtPolyRef)*bucket->npolys, BC_ALLOC_PERM);
				memcpy(bucket->bucketPoly,&polyRef, sizeof(dtPolyRef)*bucket->npolys);

				bucket->portals = (dtPolyRef*)bcAlloc(sizeof(dtPolyRef)*bucket->npolys, BC_ALLOC_PERM);
				memcpy(bucket->portals, &polyRef, sizeof(dtPolyRef));
			
				//*** Bucket Portal Neighbours******//
				bucket->nportals=m;
				bucket->portalsNeigh = (dtPolyRef*)bcAlloc(sizeof(dtPolyRef)*bucket->npolys, BC_ALLOC_PERM);
				memcpy(bucket->portalsNeigh, &neighbour, sizeof(dtPolyRef)*bucket->nportals);

				//*** Bucket Portal neiLines******//
				bucket->neiLine = (bcNeiLine*)bcAlloc(sizeof(bcNeiLine)*bucket->nportals, BC_ALLOC_PERM);
				

				//*** Bucket Portal leftVert******//
				bucket->leftVert = (float*)bcAlloc(sizeof(float)*bucket->nportals*3, BC_ALLOC_PERM);
				memcpy(bucket->leftVert, &lvert, sizeof(float)*bucket->nportals*3);

				//*** Bucket Portal rightVert******//
				bucket->rightVert = (float*)bcAlloc(sizeof(float)*bucket->nportals*3, BC_ALLOC_PERM);
				memcpy(bucket->rightVert, &rvert, sizeof(float)*bucket->nportals*3);
			}
			else
			{
				return false;
			}
		}

	return true;
}

vertex_t bcBucketMesh::bcvertex(vertex_t vparent,vertex_t v,const MaGraph& mg)
{
	int i=out_degree(v, mg);
	if (i<3)
		{
			if (i==1)
			{
				return vparent;
			}
			if (i==2)
			{
				std::pair<ma_adjacency_iterator_t, ma_adjacency_iterator_t> neighbors = boost::adjacent_vertices(v,mg);
				for(; neighbors.first != neighbors.second; ++neighbors.first)
				{
					if (*neighbors.first!=vparent)
					{
						return bcvertex(v,*neighbors.first,mg);
					}
				}
			}
		}
	else
		{
		return v;
		}
}

int bcBucketMesh::bcBuildPartialIso(dtNavMesh* mesh,const dtMeshTile* tile,dtPolyRef bcPolyRef,float* verts,int nv,std::list<polygon_t>& output)
{
	const dtMeshTile* bcPolyTile = 0;
	const dtPoly* bcPoly = 0;

	if (dtStatusFailed(mesh->getTileAndPolyByRef(bcPolyRef, &bcPolyTile, &bcPoly)))
			return 0;
	
	polygon_t isoPolygon;
	polygon_t bcPolygon;
	std::vector<point_t> IsoVectorPoints;
	std::vector<point_t> bcVectorPoints;

	int j=0;
	while (j < nv*4)
	{
		point_t IsoPolyPoint(verts[j],verts[j+1] );
		IsoVectorPoints.push_back(IsoPolyPoint);
		j=j+2;
	}
	
	for (int j = 0, nj = (int)bcPoly->vertCount; j < nj; ++j)
	{
		float* v0 = &tile->verts[bcPoly->verts[j]*3];
		double t[2];
		t[0]=double(v0[0]); t[1]=double(v0[2]);
		point_t bcPolyPoint(t[0],t[1] );
		bcVectorPoints.push_back(bcPolyPoint);
	}
	
	boost::geometry::append(isoPolygon, IsoVectorPoints);
	boost::geometry::correct(isoPolygon);
	boost::geometry::append(bcPolygon, bcVectorPoints);
    boost::geometry::correct(bcPolygon);

	float ia=boost::geometry::area(isoPolygon);
	float bca=boost::geometry::area(bcPolygon);

	IsoVectorPoints.clear();
	bcVectorPoints.clear();

	boost::geometry::difference(isoPolygon, bcPolygon, output);

	return output.size();

}

static double distanceFar(polygon_t const& poly,point_t point)
{
   
	double max_d = 0.0;
    for(auto it = boost::begin(boost::geometry::exterior_ring(poly)); it != boost::end(boost::geometry::exterior_ring(poly)); ++it)
    {
		double x = boost::geometry::get<0>(*it);
		double y = boost::geometry::get<1>(*it);
		point_t bcPolyPoint(x,y);
        double d = boost::geometry::comparable_distance(point, bcPolyPoint);
        if (d > max_d)
        {
            max_d = d;
        }
    }
	return max_d;
}

bool bcBucketMesh::bcAffectIsoToNei(dtNavMesh* mesh,const dtMeshTile* tile,float* queryPoint,dtPolyRef bcPolyRef,std::list<polygon_t> polygons,float* a,double* d)
{
	const dtMeshTile* bcPolyTile = 0;
	const dtPoly* bcPoly = 0;
	std::vector<point_t> vectorPoints;
	std::list<polygon_t> output;
	polygon_t neiPolygon;

	if (dtStatusFailed(mesh->getTileAndPolyByRef(bcPolyRef, &bcPolyTile, &bcPoly)))
			return 0;

	for (int j = 0, nj = (int)bcPoly->vertCount; j < nj; ++j)
	{
		float* v0 = &tile->verts[bcPoly->verts[j]*3];
		double t[2];
		t[0]=double(v0[0]); t[1]=double(v0[2]);
		point_t bcPolyPoint(t[0],t[1]);
		vectorPoints.push_back(bcPolyPoint);
	}
	
	boost::geometry::append(neiPolygon, vectorPoints);
	boost::geometry::correct(neiPolygon);
	float neia=boost::geometry::area(neiPolygon);

	vectorPoints.clear();
	float maxarea=0.0;
	polygon_t maxit;
	double t[2];
	t[0]=double(queryPoint[0]); t[1]=double(queryPoint[2]);
	point_t isoCenterPoint(t[0],t[1]);
	
	for (std::list<polygon_t>::iterator it=polygons.begin(); it != polygons.end(); ++it)
    {
		boost::geometry::intersection(neiPolygon,*it, output);

		for (std::list<polygon_t>::iterator iit=output.begin(); iit != output.end(); ++iit)
		{
				*a=boost::geometry::area(*it);
				*d = distanceFar(*it,isoCenterPoint);
				return true;
		}
	}

	return true;

}

void bcBucketMesh::bcBuildBucketGraph(dtNavMesh* mesh,dtNavMeshQuery* navQuery,bcBucketSet& bcset,BcGraph& cg)
{
	 boost::associative_property_map<bucket_property_t> bc_map(m_bct_prop);
	

	 std::pair<b_vertex_iterator_t, b_vertex_iterator_t> vpg;
	 b_vertex_t vgn,vg,vg0,vg1;
	 bool exist;
	 bcBucket* bc=NULL;
	 bcBucket* bucket=NULL;
	 bcBucket* bucketN=NULL;
	 const int maxPath=256;
	 int pathCount;

	 for (int i = 0; i < bcset.nbuckets; ++i)
	 {
		bucket=&bcset.buckets[i];	
		dtPolyRef* portpoly=bucket->bucketPoly;

		exist=false;
		for (vpg = vertices(cg); vpg.first != vpg.second; ++vpg.first)
		{
			bc=bc_map[*vpg.first];
			if (bc==bucket)
			{
				vg=*vpg.first;
				exist= true;
				break;
			}

		}
		if (exist)
		{
			vg0=vg;
		}
		else
		{
			vg0 = boost::add_vertex(cg);
			put(bc_map,vg0,bucket);
		}

		for (int j = 0; j < bucket->nportals; ++j)
		{
			dtPolyRef portNeigh=bucket->portalsNeigh[j];
			unsigned short int areaId;
			
			pathCount=-1;
			dtPolyRef path[maxPath];
			dtPolyRef startRef;
			dtPolyRef endRef;

			dtStatus status=bcfindRoad(mesh,*portpoly,portNeigh,&areaId,path,&pathCount,maxPath,&startRef,&endRef);
			if (status!=DT_SUCCESS)
				continue;

			if (pathCount<0)
			{
				continue;
			}

			exist= false;
			for (int t = 0; t < bcset.nbuckets; ++t)
			{
				bucketN=&bcset.buckets[t];	
		        unsigned short ref=bucketN->bcBucRef;
				if (areaId==ref)
				{
					exist=true;
					break;
				}

			}

			if (exist==false)
				continue;

			exist=false;
			for (vpg = vertices(cg); vpg.first != vpg.second; ++vpg.first)
			{
				bc=bc_map[*vpg.first];
				if (bc==bucketN)
				{
					vg=*vpg.first;
					exist= true;
					break;
				}

			}
			if (exist)
			{
				vg1=vg;
			}
			else
			{
				vg1 = boost::add_vertex(cg);
				put(bc_map,vg1,bucketN);
			}
			if (boost::edge(vg0,vg1,cg).second) continue;
			boost::add_edge(vg0,vg1,(1),cg);
		}
	}
}

dtStatus bcBucketMesh::bcSetBucketType(dtNavMesh* mesh,bcBucketSet& bcset)
{
	int n=0;
	for (int i = 0; i < bcset.nbuckets; ++i)
	{
		for (int j = 0; j < bcset.buckets[i].npolys; ++j)
		{
			dtPolyRef polyRef=bcset.buckets[i].bucketPoly[j];
			
			if (!polyRef)
				return DT_FAILURE | DT_INVALID_PARAM;

			if (!mesh->isValidPolyRef(polyRef))
				return DT_FAILURE | DT_INVALID_PARAM;

			mesh->setPolyArea(polyRef,SAMPLE_POLYAREA_CROSSROAD); 
			mesh->setPolyAreaId(polyRef,bcset.buckets[i].bcBucRef);
			n=n+1;
		}
	}

	if (bcset.nbuckets==n)
		return DT_SUCCESS;
	else 
		return DT_FAILURE;
}

dtStatus bcBucketMesh::bcBuildRoads(dtNavMesh* mesh,const dtMeshTile* tile,BcGraph& bcg)
{
	boost::associative_property_map<bucket_property_t> bc_map(m_bct_prop);
	boost::associative_property_map<line_property_t> line_map(m_ln_prop);
	dtStatus status=DT_FAILURE;

	std::pair<b_vertex_iterator_t, b_vertex_iterator_t> vg;
	b_vertex_t vg0,vg1;
	bcBucket* bc=NULL;
	
	const int maxPath=256;
	int pathCount;
	int n=0;
	

	for (vg = vertices(bcg); vg.first != vg.second; ++vg.first)
		 {
			 bc=bc_map[*vg.first];
			 bc->nverts=1024;
			 std::list<polygon_t> output;
			 bc->verts = (float*)bcAlloc(sizeof(float)*bc->nverts*4, BC_ALLOC_PERM);
	
			 dtPolyRef portpoly=bc->portals[0];
			 n++;
			 for (int j = 0; j < bc->nportals; ++j)
				{
					
					dtPolyRef portNeigh=bc->portalsNeigh[j];
					
					unsigned short int areaId;

					pathCount=-1;
					dtPolyRef path[maxPath];
					dtPolyRef startRef;
					dtPolyRef endRef;
					status=bcfindRoad(mesh,portpoly,portNeigh,&areaId,path,&pathCount,maxPath,&startRef,&endRef);
			
					if (status==DT_SUCCESS)
					{
						if (pathCount<0)
						{
							continue;
						}
						else
						{
							std::pair<b_out_edge_iterator_t, b_out_edge_iterator_t> it = boost::out_edges(*vg.first, bcg);

							for(; it.first != it.second; ++it.first)
							{
								struct bcLine* bcl=line_map[*it.first];
								bcBucket* bcn=NULL;

								bool exist=false;

								if (bcl->bucketL==bc)
								{
									bcn=bcl->bucketR;
									bcl->polyR=endRef;
									bcl->polyL=startRef;
									exist=true;
									
								}
									
								if (bcl->bucketR==bc)
								{
									bcn=bcl->bucketL;
									bcl->polyR=startRef;
									bcl->polyL=endRef;
									exist=true;
									
								}
								
								if (exist==false) continue;

								if (bcn->bcBucRef==areaId)
								{
									bcl->linePoly=path;
									bcl->npolys=pathCount;
									
									if (pathCount==0)
										continue;

									bcl->area=bcCalcRoadArea(mesh,bcl->linePoly,bcl->npolys);
									bc->neiLine[j].neiLineRef=bcl->bcLineRef;

									for (int i = 0; i < pathCount; ++i)
									{
										mesh->setPolyArea(path[i],SAMPLE_POLYAREA_ROAD); 
										mesh->setPolyAreaId(path[i],bcl->bcLineRef);
									}
								}
							}
						}

					}
					else
						continue;
			 }
		 }
	return status;
}

dtStatus bcBucketMesh::bcfindRoad(dtNavMesh* mesh,dtPolyRef polyRef,dtPolyRef neiRef,unsigned short int* resultAreaId,
					dtPolyRef* path, int* pathCount, const int maxPath,dtPolyRef* startRef,dtPolyRef* endRef)
{
	if (!maxPath)
		return DT_FAILURE | DT_INVALID_PARAM;
	
	if (!mesh->isValidPolyRef(polyRef) || !mesh->isValidPolyRef(neiRef))
		return DT_FAILURE | DT_INVALID_PARAM;

	unsigned char resultArea;
	mesh->getPolyArea(neiRef,&resultArea);
	mesh->getPolyAreaId(neiRef,resultAreaId);
	int j = 0;
	if (resultArea == SAMPLE_POLYAREA_CROSSROAD)
	{
		path[0] = neiRef;
		*pathCount = 0;
		return DT_SUCCESS;
	}
	
	unsigned short int areaId=10000;
	*startRef=neiRef;
	bcneiiter(mesh,neiRef,polyRef,path,&j,&areaId);
	*endRef=polyRef;

	if (areaId==10000)
	{
		return DT_FAILURE;
	}
	else
	{
		*pathCount = j;	
		*resultAreaId = areaId;
		return DT_SUCCESS;
	}
}

static bool exist(dtPolyRef polyRef,dtPolyRef* path,int* j)
{
	for (int i=0;i<*j;++i)
	{
		if (path[i] == polyRef)
			return true;
	}
	
	return false;
}

void bcBucketMesh::bcneiiter(dtNavMesh* mesh,dtPolyRef polyRef,dtPolyRef parentRef,dtPolyRef* path,int* j,unsigned short int* areaId)
{
	const dtMeshTile* polyTile = 0;
	const dtPoly* poly = 0;
	mesh->getTileAndPolyByRefUnsafe(polyRef, &polyTile, &poly);
	unsigned char resultArea;
	
	mesh->getPolyArea(polyRef,&resultArea);
	
	// Expand to neighbour
	if (resultArea == SAMPLE_POLYAREA_CROSSROAD)
	{
		mesh->getPolyAreaId(polyRef,areaId);
	}
	else
	{
		path[*j] = polyRef;
		*j=*j+1;
		unsigned int i = poly->firstLink;
		while ( i != DT_NULL_LINK)
		{
			const dtLink* link = &polyTile->links[i];
			dtPolyRef neighbourRef = link->ref;

			// Skip invalid neighbours and do not follow back to parent.
			if (!neighbourRef || neighbourRef == parentRef)
			{
				i = polyTile->links[i].next;
				continue;
			}

			// Skip visited neighbours
			if (exist(neighbourRef,path,j))
				continue;

			bcneiiter(mesh,neighbourRef,polyRef,path,j,areaId);
			
			i = polyTile->links[i].next;
		}
	}

}

int bcBucketMesh::bcBuildConnectivityGraph(BcGraph& bcg,ConGraph& cg)
{
	b_edge_iterator_t ei, ei_end;
	std::pair<c_vertex_iterator_t, c_vertex_iterator_t> vi;
	c_vertex_t vg,vg0,vg1;

	bool exist;

	for (boost::tie(ei, ei_end) = edges(bcg); ei != ei_end; ++ei)
	 {
		 b_vertex_t v1 = boost::source(*ei, bcg);
		 b_vertex_t v2 = boost::target(*ei, bcg);
		 exist=false;

		 for (vi = vertices(cg); vi.first != vi.second; ++vi.first)
		 {
			 ConVertexProperties& vProp = cg[*vi.first];
			 if (((vProp.source==v1)&&(vProp.target==v2))||((vProp.source==v2)&&(vProp.target==v1)))
			 {
				vg=*vi.first;
				exist=true;
				break;
			 }
		 }

		 if (exist)
		 {
			vg0=vg;
		 }
		 else
		 {
			vg0 = boost::add_vertex(ConVertexProperties(v1,v2),cg);
		 }
		 
		 std::pair<b_adjacency_iterator_t, b_adjacency_iterator_t> sneighbors = boost::adjacent_vertices(v1,bcg);
		 for (; sneighbors.first != sneighbors.second; ++sneighbors.first)
		 {
			exist=false;

			for (vi = vertices(cg); vi.first != vi.second; ++vi.first)
			{
				ConVertexProperties& vProp = cg[*vi.first];
				if (((vProp.source==v1)&&(vProp.target==*sneighbors.first))||((vProp.source==*sneighbors.first)&&(vProp.target==v1)))
				{
					vg=*vi.first;
					exist=true;
					break;
				}
			}

			if (exist)
			{
				vg1=vg;
			}
			else
			{
				vg1 = boost::add_vertex(ConVertexProperties(v1,*sneighbors.first),cg);
			}

			if (boost::edge(vg0,vg1,cg).second) continue;
			boost::add_edge(vg0,vg1,(1),cg);

		 }

		 std::pair<b_adjacency_iterator_t, b_adjacency_iterator_t> tneighbors = boost::adjacent_vertices(v2,bcg);
		 for (; tneighbors.first != tneighbors.second; ++tneighbors.first)
		 {
			exist=false;

			for (vi = vertices(cg); vi.first != vi.second; ++vi.first)
			{
				ConVertexProperties& vProp = cg[*vi.first];
				if (((vProp.source==v2)&&(vProp.target==*tneighbors.first))||((vProp.source==*tneighbors.first)&&(vProp.target==v2)))
				{
					vg=*vi.first;
					exist=true;
					break;
				}
			}

			if (exist)
			{
				vg1=vg;
			}
			else
			{
				vg1 = boost::add_vertex(ConVertexProperties(v2,*tneighbors.first),cg);
			}

			if (boost::edge(vg0,vg1,cg).second) continue;
			boost::add_edge(vg0,vg1,(1),cg);

		 }
	}
	return boost::num_vertices(cg);
}

int bcBucketMesh::bcBuildJusGraph(BcGraph& bcg,ConGraph& cg)
{
	b_edge_iterator_t ei, ei_end;
	std::pair<c_vertex_iterator_t, c_vertex_iterator_t> vi;
	c_vertex_t vg,vg0,vg1;

	bool exist;

	for (boost::tie(ei, ei_end) = edges(bcg); ei != ei_end; ++ei)
	 {
		 b_vertex_t v1 = boost::source(*ei, bcg);
		 b_vertex_t v2 = boost::target(*ei, bcg);
		 exist=false;

		 for (vi = vertices(cg); vi.first != vi.second; ++vi.first)
		 {
			 ConVertexProperties& vProp = cg[*vi.first];
			 if (((vProp.source==v1)&&(vProp.target==v2))||((vProp.source==v2)&&(vProp.target==v1)))
			 {
				vg=*vi.first;
				exist=true;
				break;
			 }
		 }

		 if (exist)
		 {
			vg0=vg;
		 }
		 else
		 {
			vg0 = boost::add_vertex(ConVertexProperties(v1,v2),cg);
		 }
		 
		 std::pair<b_adjacency_iterator_t, b_adjacency_iterator_t> sneighbors = boost::adjacent_vertices(v1,bcg);
		 for (; sneighbors.first != sneighbors.second; ++sneighbors.first)
		 {
			exist=false;

			for (vi = vertices(cg); vi.first != vi.second; ++vi.first)
			{
				ConVertexProperties& vProp = cg[*vi.first];
				if (((vProp.source==v1)&&(vProp.target==*sneighbors.first))||((vProp.source==*sneighbors.first)&&(vProp.target==v1)))
				{
					vg=*vi.first;
					exist=true;
					break;
				}
			}

			if (exist)
			{
				vg1=vg;
			}
			else
			{
				vg1 = boost::add_vertex(ConVertexProperties(v1,*sneighbors.first),cg);
			}

			if (boost::edge(vg0,vg1,cg).second) continue;

			struct bcBucket* LV1=m_bct_prop[v1]; float P1[3];dtVcopy(P1,LV1->centre);
		    struct bcBucket* RV2=m_bct_prop[v2]; float P2[3];dtVcopy(P2,RV2->centre);
			
		    struct bcBucket* NV1=m_bct_prop[*sneighbors.first];float P4[3];dtVcopy(P4,NV1->centre);
			float a = CalcAngle2(P1, P2, P1, P4);
			boost::add_edge(vg0,vg1,(1),cg);

		 }

		 std::pair<b_adjacency_iterator_t, b_adjacency_iterator_t> tneighbors = boost::adjacent_vertices(v2,bcg);
		 for (; tneighbors.first != tneighbors.second; ++tneighbors.first)
		 {
			exist=false;

			for (vi = vertices(cg); vi.first != vi.second; ++vi.first)
			{
				ConVertexProperties& vProp = cg[*vi.first];
				if (((vProp.source==v2)&&(vProp.target==*tneighbors.first))||((vProp.source==*tneighbors.first)&&(vProp.target==v2)))
				{
					vg=*vi.first;
					exist=true;
					break;
				}
			}

			if (exist)
			{
				vg1=vg;
			}
			else
			{
				vg1 = boost::add_vertex(ConVertexProperties(v2,*tneighbors.first),cg);
			}

			if (boost::edge(vg0,vg1,cg).second) continue;

			struct bcBucket* LV1=m_bct_prop[v1]; float P1[3];dtVcopy(P1,LV1->centre);
		    struct bcBucket* RV2=m_bct_prop[v2]; float P2[3];dtVcopy(P2,RV2->centre);
			
		    struct bcBucket* NV1=m_bct_prop[*tneighbors.first];float P4[3];dtVcopy(P4,NV1->centre);
			float a = CalcAngle2(P1, P2, P2, P4);
			boost::add_edge(vg0,vg1,(1),cg);

		 }
	}
	return boost::num_vertices(cg);
}

bcBucketSet* bcBucketMesh::bcAllocBucketSet()
{
	bcBucketSet* bset = (bcBucketSet*)bcAlloc(sizeof(bcBucketSet), BC_ALLOC_PERM);
	memset(bset, 0, sizeof(bcBucketSet));
	return bset;
}

int bcBucketMesh::buildDistMat(ConGraph& cg,std::vector< std::vector<int> >& distMat)
{
	std::pair<c_vertex_iterator_t, c_vertex_iterator_t> vi;
	std::pair<c_vertex_iterator_t, c_vertex_iterator_t> vj;
	
	//*******************dijkstra********************************************************
	boost::property_map<ConGraph,boost::edge_weight_t>::type weightmap = get(boost::edge_weight, cg);

    for (vi = vertices(cg); vi.first != vi.second; ++vi.first)
	{
		std::vector<c_vertex_t> p(num_vertices(cg));
		std::vector<int> d(num_vertices(cg));
		boost::dijkstra_shortest_paths(cg, *vi.first,
               predecessor_map(boost::make_iterator_property_map(p.begin(), get(boost::vertex_index, cg))).
               distance_map(boost::make_iterator_property_map(d.begin(), get(boost::vertex_index, cg))));
		for (vj = vertices(cg); vj.first != vj.second; ++vj.first)
		{
			distMat[*vi.first][*vj.first]= d[*vj.first];
		}
		p.clear();d.clear();
	}
    //*******************dijkstra********************************************************
	
	return boost::num_vertices(cg);
}

int bcBucketMesh::TotalDepth(int node,std::vector< std::vector<int> >& distMat,int k,bool t) 
{
	if (t==true)
	{
		int TD=0;
		for (int i=0;i<k;i++)
		{
			TD=TD+distMat[node][i];
		}
		return TD;
	}
	else
	{
		int LD=0;
		for (int i=0;i<k;i++)
		{
			if (distMat[node][i]<=3) 
				LD=LD+distMat[node][i];
		}
		return LD;
	}
}

//t==0 Local Depth
//t==1 global Depth

float bcBucketMesh::MeanDepth(int node,std::vector< std::vector<int> >& distMat, int k,bool t) 
{
	float MD=0.0;
	if (t==true)
	{
		MD=TotalDepth(node,distMat,k,t)/(k-1);
		return MD;
	}
	else
	{
		MD=TotalDepth(node,distMat,k,t)/(2);
		return MD;
	}

}

//t==0 Local Integration
//t==1 global Integration
float bcBucketMesh::Integration(int node,std::vector< std::vector<int> >& distMat, int k,bool t )
{ 
	float MD=MeanDepth(node,distMat,k,t);
	float RA=1.0;
	if (t==true)
	{
		RA=2*(MD-1)/(float(k-2));
	}
	else
	{
		RA=2*(MD-1)/(float(k-2));
	}
	float RRA=1/RA;
	return RRA;
}

int bcBucketMesh::Connectivity(int node,ConGraph& cg)
{ 
	int i=out_degree(node, cg);
	return i;
}

int bcBucketMesh::VConnectivity(int node,BcGraph& bcg)
{ 
	int i=out_degree(node, bcg);
	return i;
}

bool bcBucketMesh::bcLinePolyIntersection(struct bcBucket* L,struct bcBucket* R,polygon_t& poly)
{
	bool visible=false;
	line_t ls,output;
	segment_t s;
	point_t p0(L->centre[0],L->centre[2]);
	point_t p1(R->centre[0],R->centre[2]);

	boost::geometry::append(ls, p0); 
	boost::geometry::append(ls, p1); 

	boost::geometry::intersection(ls, poly, output);

	if (output.size()==0)
		visible=true;
	else
		visible=false;

	return visible;
}

void bcBucketMesh::bcCalcBetweeness(ConGraph& cg)
{
	boost::shared_array_property_map<double, boost::property_map<ConGraph, boost::vertex_index_t>::const_type>
		centrality_map(boost::num_vertices(cg), get(boost::vertex_index, cg));

	boost::relative_betweenness_centrality(cg,centrality_map);
}

int bcBucketMesh::bcCalcSSParam(ConGraph& cg,BcGraph& bcg,std::vector< std::vector<int> >& distMat, int k,polygon_t& poly) 
{
	boost::associative_property_map<line_property_t> line_map(m_ln_prop);

	boost::shared_array_property_map<double, boost::property_map<ConGraph, boost::vertex_index_t>::const_type>
		centrality_map(boost::num_vertices(cg), get(boost::vertex_index, cg));

	boost::brandes_betweenness_centrality(cg,centrality_map);
	
	b_edge_iterator_t ei, ei_end;
	std::pair<c_vertex_iterator_t, c_vertex_iterator_t> vi;
	c_vertex_t vg;
	struct bcLine* bcl=NULL;
	bool exist;
	int n=0;

	for (boost::tie(ei, ei_end) = edges(bcg); ei != ei_end; ++ei)
	{
		 b_vertex_t v1 = boost::source(*ei, bcg);
		 b_vertex_t v2 = boost::target(*ei, bcg);

		 struct bcBucket* L=m_bct_prop[v1]; 
		 struct bcBucket* R=m_bct_prop[v2];

		 exist=false;
		 for (vi = vertices(cg); vi.first != vi.second; ++vi.first)
	     {
			 ConVertexProperties& vProp = cg[*vi.first];
			 if (((vProp.source==v1)&&(vProp.target==v2))||((vProp.source==v2)&&(vProp.target==v1)))
			 {
				vg=*vi.first;
				exist=true;
				break;
			 }
	     }
		 if (exist)
		 {
			bcl = (bcLine*)bcAlloc(sizeof(bcLine), BC_ALLOC_PERM);
		
			bcl->bcLineRef=vg;
			bcl->bucketL=L;
			bcl->bucketR=R;
			bcl->connectivity=Connectivity(vg,cg);
			bcl->betweeness=(centrality_map[vg]);
			bcl->visibility=bcLinePolyIntersection(L,R,poly);
			const float* a=segmentMidpoint(L->centre,R->centre);
			dtVcopy(bcl->midpoint, a);
			n++;
			put(line_map,*ei,bcl);
		 }
	}

	return n;
}

int bcBucketMesh::bcCalcSSParamm(ConGraph& cg,BcGraph& bcg,polygon_t& poly) 
{
	boost::associative_property_map<line_property_t> line_map(m_ln_prop);

	boost::shared_array_property_map<double, boost::property_map<ConGraph, boost::vertex_index_t>::const_type>
		centrality_map(boost::num_vertices(cg), get(boost::vertex_index, cg));

	boost::brandes_betweenness_centrality(cg,centrality_map);
	
	b_edge_iterator_t ei, ei_end;
	std::pair<c_vertex_iterator_t, c_vertex_iterator_t> vi;
	c_vertex_t vg;
	struct bcLine* bcl=NULL;
	bool exist;
	int n=0;

	for (boost::tie(ei, ei_end) = edges(bcg); ei != ei_end; ++ei)
	{
		 b_vertex_t v1 = boost::source(*ei, bcg);
		 b_vertex_t v2 = boost::target(*ei, bcg);

		 struct bcBucket* L=m_bct_prop[v1]; 
		 struct bcBucket* R=m_bct_prop[v2];

		 exist=false;
		 for (vi = vertices(cg); vi.first != vi.second; ++vi.first)
	     {
			 ConVertexProperties& vProp = cg[*vi.first];
			 if (((vProp.source==v1)&&(vProp.target==v2))||((vProp.source==v2)&&(vProp.target==v1)))
			 {
				vg=*vi.first;
				exist=true;
				break;
			 }
	     }
		 if (exist)
		 {
			bcl = (bcLine*)bcAlloc(sizeof(bcLine), BC_ALLOC_PERM);
		
			bcl->bcLineRef=vg;
			bcl->bucketL=L;
			bcl->bucketR=R;
			bcl->connectivity=Connectivity(vg,cg);
			bcl->betweeness=(centrality_map[vg]);
			bcl->visibility=bcLinePolyIntersection(L,R,poly);
			const float* a=segmentMidpoint(L->centre,R->centre);
			dtVcopy(bcl->midpoint, a);
			n++;
			put(line_map,*ei,bcl);
			
		 }
	}

	return n;
}

int bcBucketMesh::bcBuildVisibilityGraph(BcGraph& cgv,BcGraph& bcg) 
{
	boost::associative_property_map<bucket_property_t> bc_map(m_bct_prop);
	boost::associative_property_map<bucket_property_t> bctv_map(m_bctv_prop);
	std::pair<b_vertex_iterator_t, b_vertex_iterator_t> vvg;
	b_edge_iterator_t ei, ei_end;
	bcBucket* bc=NULL;
	bcBucket* bcv=NULL;
	bool exist;
	b_vertex_t vg,vg0,vg1;

	for (boost::tie(ei, ei_end) = edges(bcg); ei != ei_end; ++ei)
	 {
		struct bcLine* L=m_ln_prop[*ei]; 
		if (L->visibility)
		{
			b_vertex_t v1 = boost::source(*ei, bcg);
			b_vertex_t v2 = boost::target(*ei, bcg);
			exist=false;
			bc=bc_map[v1];
			for (vvg = vertices(cgv); vvg.first != vvg.second; ++vvg.first)
			{
				bcv=bctv_map[*vvg.first];
				if (bc==bcv)
				{
					vg=*vvg.first;
					exist =true;
					break;
				}
			}

			if (exist)
			{
				vg0=vg;
			}
			else
			{
				vg0 = boost::add_vertex(cgv);
				put(bctv_map,vg0,bc);
			}
			exist= false;
			bc=bc_map[v2];
			for (vvg = vertices(cgv); vvg.first != vvg.second; ++vvg.first)
			{
				bcv=bctv_map[*vvg.first];
				if (bc==bcv)
				{
					vg=*vvg.first;
					exist =true;
					break;
				}
			}
			if (exist)
			{
				vg1=vg;
			}
			else
			{
				vg1 =boost::add_vertex(cgv);
				put(bctv_map,vg1,bc);
			}

			if (boost::edge(vg0,vg1,cgv).second) continue;
			boost::add_edge(vg0,vg1,(1),cgv);
			
		}
	}

	return boost::num_vertices(cgv);
}

int bcBucketMesh::buildDistMatV(BcGraph& cgv,std::vector< std::vector<int> >& distMat)
{
	std::pair<b_vertex_iterator_t, b_vertex_iterator_t> vi;
	std::pair<b_vertex_iterator_t, b_vertex_iterator_t> vj;
	
	//*******************dijkstra********************************************************
	boost::property_map<ConGraph,boost::edge_weight_t>::type weightmap = get(boost::edge_weight, cgv);
    std::vector<c_vertex_t> p(num_vertices(cgv));
	std::vector<int> d(num_vertices(cgv));
    for (vi = vertices(cgv); vi.first != vi.second; ++vi.first)
	{
		boost::dijkstra_shortest_paths(cgv, *vi.first,
               predecessor_map(boost::make_iterator_property_map(p.begin(), get(boost::vertex_index, cgv))).
               distance_map(boost::make_iterator_property_map(d.begin(), get(boost::vertex_index, cgv))));
		for (vj = vertices(cgv); vj.first != vj.second; ++vj.first)
		{
			distMat[*vi.first][*vj.first]= d[*vj.first];
		}
	}
    //*******************dijkstra********************************************************
	
	return boost::num_vertices(cgv);
}

void bcBucketMesh::bcCalcSSParamVG(dtNavMesh* mesh,BcGraph& cgv,std::vector< std::vector<int> >& distMat, int k)
{
	boost::associative_property_map<bucket_property_t> bcv_map(m_bctv_prop);
	std::pair<b_vertex_iterator_t, b_vertex_iterator_t> vi;
	
	bcBucket* bc=NULL;
	
	for (vi = vertices(cgv); vi.first != vi.second; ++vi.first)
	{
		bc=bcv_map[*vi.first];
		bc->connectivity=VConnectivity(*vi.first,cgv);
		bc->globalIntegration=Integration(*vi.first,distMat,k,true);
		bc->localIntegration=Integration(*vi.first,distMat,k,false);
		bc->area=bcCalcRoadArea(mesh,bc->bucketPoly,1);
	}
}

void bcBucketMesh::bcCalcIntersectionBucParam(BcGraph& cgv)
{
	boost::associative_property_map<bucket_property_t> bc_map(m_bct_prop);
	std::pair<b_vertex_iterator_t, b_vertex_iterator_t> vi;
	
	bcBucket* bc=NULL;
	
	for (vi = vertices(cgv); vi.first != vi.second; ++vi.first)
	{
		bc=bc_map[*vi.first];
		bc->connectivity=VConnectivity(*vi.first,cgv);
	}
}

void bcBucketMesh::bcCalcSSParammVG(dtNavMesh* mesh,BcGraph& cgv)
{
	boost::associative_property_map<bucket_property_t> bcv_map(m_bctv_prop);
	std::pair<b_vertex_iterator_t, b_vertex_iterator_t> vi;
	
	bcBucket* bc=NULL;
	
	for (vi = vertices(cgv); vi.first != vi.second; ++vi.first)
	{
		bc=bcv_map[*vi.first];
		bc->connectivity=VConnectivity(*vi.first,cgv);
		bc->area=bcCalcRoadArea(mesh,bc->bucketPoly,1);
	}
}

float bcBucketMesh::bcCalcArea(const dtMeshTile* tile,const dtPoly* poly)
{
	polygon_t p;
	std::vector<point_t> vectorPoints;
    int i=0;
	
	for (int j = 0, nj = (int)poly->vertCount; j < nj; ++j)
	{
		float* v0 = &tile->verts[poly->verts[j]*3];
		double t[2];
		t[0]=double(v0[0]); t[1]=double(v0[2]);
		point_t polyPoint(t[0],t[1] );
		vectorPoints.push_back(polyPoint);
	}
	
	boost::geometry::append(p, vectorPoints);
    boost::geometry::correct(p);
	vectorPoints.clear();
	float a=boost::geometry::area(p);
	return a;
}

float bcBucketMesh::bcCalcRoadArea(dtNavMesh* mesh,dtPolyRef* linePoly,int n)
{
	float area=0.0;
	const dtMeshTile* polyTile = 0;
	const dtPoly* poly = 0;

	for (int i = 0; i < n; ++i)
	{
		dtPolyRef polyRef=linePoly[i];
		
		if (dtStatusFailed(mesh->getTileAndPolyByRef(polyRef, &polyTile, &poly)))
			return 0;
		area=area+bcCalcArea(polyTile,poly);
	}

	return area;
}

float bcBucketMesh::bcCalcAngle(float* p1,float* p2,float* p3)
{
	float dx21 = p2[0]-p1[0];
	float dy21 = p2[0]-p1[2];

	float dx31 = p3[0]-p1[0];
	float dy31 = p3[0]-p1[2];

    float dot = (dx21 * dx31 + dy21 * dy31);

    float abSqr = dx21 * dx21 + dy21 * dy21;
    float cbSqr = dx31 * dx31 + dy31 * dy31;
  
    float cosSqr = dot * dot / abSqr / cbSqr;

    float cos2 = 2 * cosSqr - 1;

    const float pi = 3.141592f;

    float alpha2 =
        (cos2 <= -1) ? pi :
        (cos2 >= 1) ? 0 :
        acosf(cos2);

    float rslt = alpha2 / 2;

    float rs = rslt * 180. / pi;

    if (dot < 0)
        rs = 180 - rs;

	rs = (float)floor(rs);
    return (180 - rs);
}

static void normalizeVector(std::vector<float>& vect)
{
	int n=vect.size();
	float minPen = FLT_MAX;
	float maxPen = -FLT_MAX;
	for (int i = 0; i < n; ++i)
	{
		minPen = dtMin(minPen, vect[i]);
		maxPen = dtMax(maxPen, vect[i]);
	}
	const float penRange = maxPen-minPen;
	const float s = penRange > 0.001f ? (1.0f / penRange) : 1;
	for (int i = 0; i < n; ++i)
		vect[i] = dtClamp((vect[i]-minPen)*s, 0.0f, 1.0f);
}

static void normalizeIntVector(std::vector<int> vect,std::vector<float>& nvect)
{
	int n=vect.size();
	int maxPen = 0;
	for (int i = 0; i < n; ++i)
	{
		maxPen = dtMax(maxPen,vect[i]);
	}
	
	for (int i = 0; i < n; ++i)
		nvect[i] = vect[i]*(100/maxPen);
}

static void normalizeParams(std::vector<float>& angles,std::vector<float>& dists)
{
	normalizeVector(angles);
	normalizeVector(dists);
}

static void normalizePPParams(std::vector<float>& area,std::vector<float>& dist,std::vector<int> con,std::vector<float>& ncon)
{
	normalizeVector(area);
	normalizeVector(dist);
	normalizeIntVector(con,ncon);
}

float bcBucketMesh::bcCalcdist(float* p1,float* p2,float* p3)
{
	float d12=bcDistancePtPtSqr2D(p1,p2);
	float d23=bcDistancePtPtSqr2D(p2,p3);
	return d12+d23;
}

static bool calcVectProd(float* A,float* B,float* C)
{
	float v1[3];float v2[3];float v3[3]; 
	
	v1[0]=A[0]-B[0]; v2[0]=C[0]-B[0];
	v1[1]=A[1]-B[1]; v2[1]=C[1]-B[1];
	v1[2]=A[2]-B[2]; v2[2]=C[2]-B[2];

	v3[0]=v1[1]*v2[2]-v1[2]*v2[1];
	v3[1]=v1[2]*v2[0]-v1[0]*v2[2];
	v3[2]=v1[0]*v2[1]-v1[1]*v2[0];

	if (v3[1]>0)
	{
		return true;
	}
	else 
	{
		return false;
	}
}

int bcBucketMesh::bcCalcCost(std::vector<float>& costVect,std::vector<float> angles,std::vector<float> dists,
	                         unsigned int fam, unsigned int expl)
{
	float p=0.1;
	for (int i=0;i<angles.size();++i)
		{
			float t=(p*angles[i]+dists[i]);
			costVect.push_back(t);
		}

	if (costVect.size()==angles.size())
		return costVect.size();
	else 
		return 0;
}

int bcBucketMesh::bcCalcBenefit(std::vector<float>& costVect,std::vector<float> angles,std::vector<float> dists)
{
	float p=0.1;

	for (int i=0;i<angles.size();++i)
		{
			float t=(p*angles[i]+dists[i]);
			costVect.push_back(t);
		}

	if (costVect.size()==angles.size())
		return costVect.size();
	else 
		return 0;
}

static int rand_oi() 
{
	float r = ((double) rand() / (RAND_MAX));
	return r;
}

static int bcCalcChoiceBenefit(std::vector<float>& costVect,std::vector<float> area,std::vector<float> dist,std::vector<float> con)
{
	float p[3];
	for (int i=0;i<3;++i)
	{
		float t=rand_oi();
		if (t==0)
			p[i]=t+1;
		p[i]=t;
	}

	for (int i=0;i<area.size();++i)
		{
			float t=(p[0]*area[i]+p[1]*dist[i]+p[2]*con[i]);
			costVect.push_back(t);
		}

	if (costVect.size()==area.size())
		return costVect.size();
	else 
		return 0;
}

void bcBucketMesh::bcCalcPath(dtNavMesh* mesh,bcBucket* bc,std::map<int, int> mapDensity,dtPolyRef parentRef,unsigned short parentLineRef,float* startPos,float* endPos,
	                          dtPolyRef endRef,unsigned int fam, unsigned int expl,struct bcLineStep* path,int* j,unsigned short p,unsigned short* t)    //p:line count
{
	boost::associative_property_map<bucket_property_t> bcv_map(m_bct_prop);                                                        
	boost::associative_property_map<line_property_t> line_map(m_ln_prop);

	std::pair<b_vertex_iterator_t, b_vertex_iterator_t> vg;
	std::vector<float> angles;
	std::vector<float> densities;
	std::vector<float> dists;
	std::vector<float> costV;
	std::vector<struct bcLine*> lines;
	
	unsigned char endArea;
	mesh->getPolyArea(endRef,&endArea);
	unsigned short int endAreaId;
	mesh->getPolyAreaId(endRef,&endAreaId);

	float parentPos[3];
	getPolyCenter(mesh,parentRef, parentPos);

	if (endArea==2)
		{
			if (bc->bcBucRef==endAreaId)
				return;
		}

	if (p!=9999)
		*j=*j+1;

	bcBucket* bcn=NULL;
	for (vg = vertices(m_bucketGraph); vg.first != vg.second; ++vg.first)
	{
		bcBucket* bct= bcv_map[*vg.first];
		
		if (bc!=bct)
			continue;

		unsigned short minRef=-1;
		float minAng=0.0;
		struct bcLine* bclmin=NULL;
		struct bcLine* bcl=NULL;
		float angle=0.0;
		float dist=0.0;

		std::pair<b_out_edge_iterator_t, b_out_edge_iterator_t> it = boost::out_edges(*vg.first,m_bucketGraph);
		std::pair<b_out_edge_iterator_t, b_out_edge_iterator_t> itt = boost::out_edges(*vg.first,m_bucketGraph);
		int c=0;
		for(; itt.first != itt.second; ++itt.first)
		{
			bcl=line_map[*itt.first];
			if (bcl->bcLineRef==parentLineRef)
				continue;
			float density= mapDensity[bcl->bcLineRef]/bcl->area;
			if ((fam==bcl->interval)&&(density<SeuilMax))
				c=c+1;;
		}
		for(; it.first != it.second; ++it.first)
		{
			bcl=line_map[*it.first];
			if (bcl->bcLineRef==parentLineRef)
				continue;
			if (bcl->bcLineRef==endAreaId)
			{
				path[*j].bcLineRef = bcl->bcLineRef;
				path[*j].interBuc = endRef;

				if (p==9999)
		            *j=*j+1;
				return;
			}
			if (fam!=1)
				if ((c>1)&&(expl==0)&&(fam>bcl->interval))
					continue;

			if ((bcl->bucketL==bc)&&(bcl->bucketR->bucketPoly[0]!=parentRef))
				bcn=bcl->bucketR;
			else if ((bcl->bucketR==bc)&&(bcl->bucketL->bucketPoly[0]!=parentRef))
				bcn=bcl->bucketL;
			else
				continue;

			dist=bcCalcdist(bc->centre,bcn->centre,endPos);
			dists.push_back(dist);
			angle=bcCalcAngle(bc->centre,bcn->centre,parentPos);
			angles.push_back(angle);
			lines.push_back(bcl);
		}

		normalizeParams(angles,dists);

		int s = bcCalcCost(costV,angles,dists,fam,expl);

		if (s==0)	
			return;

		minAng=costV[0];
		for (int i=0;i<costV.size();++i)
		{
			if (minAng>=costV[i])
			{
				minAng=costV[i];
				bclmin=lines[i];
				if (angles[i]>=0.25) 
					*t=*t+1;

			}
		}

		path[*j].bcLineRef=bclmin->bcLineRef;
		
		if (bc==bclmin->bucketL)
		{
			bcn=bclmin->bucketR;
			path[*j].interBuc = bcn->bucketPoly[0];
		}
		if (bc==bclmin->bucketR)
		{
			bcn=bclmin->bucketL;
			path[*j].interBuc = bcn->bucketPoly[0];
		}

		costV.clear();
		angles.clear();
		dists.clear();
		lines.clear();
		break;
	}
}

static int rand_o(int n) 
{ 
	srand(time(NULL));
    int nbgen=rand()%n;    
    
    return nbgen;
}

void bcBucketMesh::bcCalcChoice(dtNavMesh* mesh,bcBucket* bc,std::map<int, int> mapDensity,dtPolyRef parentRef,float* startPos,float* endPos,
	                          dtPolyRef endRef,struct bcLineStep* path,int* j,unsigned short p,unsigned short* t)    
{
	boost::associative_property_map<bucket_property_t> bcv_map(m_bct_prop);                                                        
	boost::associative_property_map<line_property_t> line_map(m_ln_prop);

	std::pair<b_vertex_iterator_t, b_vertex_iterator_t> vg;
	std::vector<int> connectivity;
	std::vector<float> nconnectivity;
	std::vector<float> densities;
	std::vector<float> dists;
	std::vector<float> isoArea;
	std::vector<float> lineSight;
	std::vector<float> costV;
	std::vector<struct bcLine*> lines;
	
	unsigned char endArea;
	mesh->getPolyArea(endRef,&endArea);
	unsigned short int endAreaId;
	mesh->getPolyAreaId(endRef,&endAreaId);

	float parentPos[3];
	getPolyCenter(mesh,parentRef, parentPos);

	if (endArea==2)
		{
			if (bc->bcBucRef==endAreaId)
				return;
		}

	if (p!=9999)
		*j=*j+1;

	bcBucket* bcn=NULL;
	for (vg = vertices(m_bucketGraph); vg.first != vg.second; ++vg.first)
	{
		bcBucket* bct= bcv_map[*vg.first];
		
		if (bc!=bct)
			continue;

		unsigned short minRef=-1;
		int maxCost=0;
		struct bcLine* bclmin=NULL;
		struct bcLine* bcl=NULL;
		float dist=0.0;

		std::pair<b_out_edge_iterator_t, b_out_edge_iterator_t> it = boost::out_edges(*vg.first,m_bucketGraph);
		std::pair<b_out_edge_iterator_t, b_out_edge_iterator_t> itt = boost::out_edges(*vg.first,m_bucketGraph);
		
		int c=0;
		for(; itt.first != itt.second; ++itt.first)
		{
			bcl=line_map[*itt.first];
			float density= mapDensity[bcl->bcLineRef]/bcl->area;
			if ((bcl->visibility)&&(density>SeuilMin))
				c=c+1;;
		}
		
		for(; it.first != it.second; ++it.first)
		{
			bcl=line_map[*it.first];

			if (bcl->bcLineRef==endAreaId)
			{
				path[*j].bcLineRef = bcl->bcLineRef;
				path[*j].interBuc = endRef;

				if (p==9999)
		            *j=*j+1;
				return;
			}

			
			if ((c>1)&&(bcl->visibility)==false)
					continue;

			if ((bcl->bucketL==bc)&&(bcl->bucketR->bucketPoly[0]!=parentRef))
				bcn=bcl->bucketR;
			else if ((bcl->bucketR==bc)&&(bcl->bucketL->bucketPoly[0]!=parentRef))
				bcn=bcl->bucketL;
			else
				continue;

			float a=0.0;
			float d=0.0;

			for (int i=0;i<bct->nportals;++i)
			{
				if (bct->neiLine[i].neiLineRef==bcl->bcLineRef)
				{
					a = bct->neiLine[i].areaIso;
					d = float(bct->neiLine[i].farDistIso);
					break;
				}
			
			}

			dist=bcCalcdist(bc->centre,bcn->centre,endPos);
			dists.push_back(dist);
			isoArea.push_back(a);
			lineSight.push_back(d);
			connectivity.push_back(bcn->connectivity);
			lines.push_back(bcl);
		}

		int s=0;
		unsigned int e=rand_o(1);

		if (e==0)
		{
			normalizePPParams(isoArea,lineSight,connectivity,nconnectivity);
			s = bcCalcChoiceBenefit(costV,isoArea,lineSight,nconnectivity);
		}
		else
		{
			normalizePPParams(isoArea,dists,connectivity,nconnectivity);
			s = bcCalcChoiceBenefit(costV,isoArea,dists,nconnectivity);
		}

		if (s==0)	
			return;


		maxCost=costV[0];
		for (int i=0;i<costV.size();++i)
		{
			if (maxCost<=costV[i])
			{
				maxCost=costV[i];
				bclmin=lines[i];
			}
		}

		path[*j].bcLineRef=bclmin->bcLineRef;
		if (bc==bclmin->bucketL)
		{
			bcn=bclmin->bucketR;
			path[*j].interBuc = bcn->bucketPoly[0];
		}
		if (bc==bclmin->bucketR)
		{
			bcn=bclmin->bucketL;
			path[*j].interBuc = bcn->bucketPoly[0];
		}

		nconnectivity.clear();
		connectivity.clear();
		isoArea.clear();
		lineSight.clear();
		lines.clear();
		break;
	}
}

void bcBucketMesh::bcCalcPathLight(dtNavMesh* mesh,bcBucket* bc,dtPolyRef parentRef,unsigned short parentLineRef,float* startPos,float* endPos,
	                          dtPolyRef endRef,unsigned int fam, unsigned int expl,struct bcLineStep* path,int* j,unsigned short p,unsigned short* t)    //p:line count
{

	boost::associative_property_map<bucket_property_t> bcv_map(m_bct_prop);                                                        
	boost::associative_property_map<line_property_t> line_map(m_ln_prop);

	std::pair<b_vertex_iterator_t, b_vertex_iterator_t> vg;
	std::vector<float> angles;
	std::vector<float> dists;
	std::vector<float> costV;
	std::vector<struct bcLine*> lines;
	
	unsigned char endArea;
	mesh->getPolyArea(endRef,&endArea);
	unsigned short int endAreaId;
	mesh->getPolyAreaId(endRef,&endAreaId);

	float parentPos[3];
	getPolyCenter(mesh,parentRef, parentPos);

	if (endArea==2)
		{
			if (bc->bcBucRef==endAreaId)
				return;
		}

	if (p!=9999)
		*j=*j+1;
	
	bcBucket* bcn=NULL;
	for (vg = vertices(m_bucketGraph); vg.first != vg.second; ++vg.first)
	{
		bcBucket* bct= bcv_map[*vg.first];
		
		if (bc!=bct)
			continue;

		unsigned short minRef=-1;
		float minAng=0.0;
		struct bcLine* bclmin=NULL;
		struct bcLine* bcl=NULL;
		float angle=0.0;
		float dist=0.0;

		std::pair<b_out_edge_iterator_t, b_out_edge_iterator_t> it = boost::out_edges(*vg.first,m_bucketGraph);
		std::pair<b_out_edge_iterator_t, b_out_edge_iterator_t> itt = boost::out_edges(*vg.first,m_bucketGraph);
		int c=0;
		for(; itt.first != itt.second; ++itt.first)
		{
			bcl=line_map[*itt.first];
			if (bcl->bcLineRef==parentLineRef)
				continue;
			if (bcl->bcLineRef==endAreaId)
			{
				path[*j].bcLineRef = bcl->bcLineRef;
				path[*j].interBuc = endRef;

				if (p==9999)
		            *j=*j+1;
				return;
			}
			if (bcl->interval>=fam)
				c=c+1;;
		}

		for(; it.first != it.second; ++it.first)
		{
			bcl=line_map[*it.first];
			if (bcl->bcLineRef==parentLineRef)
				continue;
			if ((c>1)&&(expl==0)&&(fam>bcl->interval))
				continue;

			if ((bcl->bucketL==bc)&&(bcl->bucketR->bucketPoly[0]!=parentRef))
				bcn=bcl->bucketR;
			else if ((bcl->bucketR==bc)&&(bcl->bucketL->bucketPoly[0]!=parentRef))
				bcn=bcl->bucketL;
			else
				continue;

			dist=bcCalcdist(bc->centre,bcn->centre,endPos);
			dists.push_back(dist);
			angle=bcCalcAngle(bc->centre,bcn->centre,parentPos);
			angles.push_back(angle);
			lines.push_back(bcl);
		}

		normalizeParams(angles,dists);

		int s = bcCalcCost(costV,angles,dists,fam,expl);

		//if (s==0)	
		//	return;


		minAng=costV[0];
		int minidx=0;
		for (int i=0;i<costV.size();++i)
		{
			if (minAng>=costV[i])
			{
				minAng=costV[i];
				bclmin=lines[i];
				minidx=i;
			}
		}

		if (angles[minidx]>=0.25) 
			*t=*t+1;

		dtPolyRef pRef; 
		path[*j].bcLineRef=bclmin->bcLineRef;
		parentLineRef = bclmin->bcLineRef;
		if (bc==bclmin->bucketL)
		{
			pRef=bcl->bucketL->bucketPoly[0];
			bcn=bclmin->bucketR;
			path[*j].interBuc = bcn->bucketPoly[0];
			
		}
		if (bc==bclmin->bucketR)
		{
			pRef=bcl->bucketR->bucketPoly[0];
			bcn=bclmin->bucketL;
			path[*j].interBuc = bcn->bucketPoly[0];
		}

		costV.clear();
		angles.clear();
		dists.clear();
		lines.clear();

		bcCalcPathLight(mesh,bcn,pRef,parentLineRef,startPos,endPos,endRef,fam,expl,path,j,p,t);
		break;
	}
	
}

dtStatus bcBucketMesh::bcFindPathLight(dtNavMesh* mesh,BcGraph& bcg,dtPolyRef parentRef,dtPolyRef startRef, dtPolyRef endRef,float* startPos,
	                              float* endPos,unsigned int fam, unsigned int expl,struct bcLineStep* path, int* pathCount,unsigned short* t)
{
	boost::associative_property_map<bucket_property_t> bcv_map(m_bct_prop);
	boost::associative_property_map<line_property_t> line_map(m_ln_prop);

	std::pair<b_vertex_iterator_t, b_vertex_iterator_t> vg;
	b_edge_iterator_t ei, ei_end;

	bcBucket* bc=NULL;
	bool top=false;

	*pathCount = 0;
	unsigned short p=0;

	if (!startRef || !endRef)
		return DT_FAILURE | DT_INVALID_PARAM;

	if (!mesh->isValidPolyRef(startRef) || !mesh->isValidPolyRef(endRef))
		return DT_FAILURE | DT_INVALID_PARAM;

	unsigned short int startAreaId;
	unsigned char startArea;
	mesh->getPolyAreaId(startRef,&startAreaId);
	mesh->getPolyArea(startRef,&startArea);

	unsigned short int endAreaId;
	unsigned char endArea;
	mesh->getPolyAreaId(endRef,&endAreaId);
	mesh->getPolyArea(endRef,&endArea);

	if (startRef == endRef)
	{
		path[0].bcLineRef = startAreaId;
		path[0].interBuc = endRef;
		*pathCount = 1;
		return DT_SUCCESS;
	}

	int j=0;
	if (startArea==2)  
	{
		bcBucket* bcn=NULL;
		for (vg = vertices(bcg); vg.first != vg.second; ++vg.first)
		{
			bcn= bcv_map[*vg.first];
			if (bcn->bcBucRef!=startAreaId)
				continue;
			p=9999;
			bcCalcPathLight(mesh,bcn,parentRef,0,startPos,endPos,endRef,fam,expl,path,&j,p,t);
			*pathCount = 1;	
			return DT_SUCCESS;
			
		}
	}
	else
	{
		if (endArea==1)
		{
			if (startAreaId == endAreaId)
			{
				path[0].bcLineRef = startAreaId;
				path[0].interBuc = endRef;
				*pathCount = 1;
				return DT_SUCCESS;
			}
		}
	
		bcBucket* bcnL=NULL;
		bcBucket* bcnR=NULL;
		struct bcLine* bcl=NULL;
		bool exist=false;
		for (boost::tie(ei, ei_end) = edges(bcg); ei != ei_end; ++ei)
		{
			bcl=line_map[*ei];
			
			if (bcl->bcLineRef!=startAreaId)
				continue;
			
			exist=true;

			bcnL=bcl->bucketL;
			bcnR=bcl->bucketR;

			break;
		}

		if (exist==false)
			return DT_FAILURE;

		float distL=bcCalcdist(startPos,bcnL->centre,endPos);
		float distR=bcCalcdist(startPos,bcnR->centre,endPos);

		if (distL<distR)
		{
			path[0].bcLineRef = bcl->bcLineRef;
			path[0].interBuc = bcnL->bucketPoly[0];
			bcCalcPathLight(mesh,bcnL,parentRef,bcl->bcLineRef,startPos,endPos,endRef,fam,expl,path,&j,p,t);
		}
		else
		{
			path[0].bcLineRef = bcl->bcLineRef;
			path[0].interBuc = bcnR->bucketPoly[0];
			bcCalcPathLight(mesh,bcnR,parentRef,bcl->bcLineRef,startPos,endPos,endRef,fam,expl,path,&j,p,t);
		}
    }

	if (j==0)
	{
		return DT_FAILURE;
	}
	else
	{
		*pathCount = j+1;	
		return DT_SUCCESS;
	}

}

dtStatus bcBucketMesh::bcFindPath(dtNavMesh* mesh,BcGraph& bcg,std::map<int, int> mapDensity,dtPolyRef parentRef,dtPolyRef startRef, dtPolyRef endRef,float* startPos,
	                              float* endPos,unsigned int fam, unsigned int expl,struct bcLineStep* path, int* pathCount,unsigned short* t)
{
	boost::associative_property_map<bucket_property_t> bcv_map(m_bct_prop);
	boost::associative_property_map<line_property_t> line_map(m_ln_prop);

	std::pair<b_vertex_iterator_t, b_vertex_iterator_t> vg;
	b_edge_iterator_t ei, ei_end;

	bcBucket* bc=NULL;
	bool top=false;

	*pathCount = 0;
	*t=0;
	unsigned short p=0;

	if (!startRef || !endRef)
		return DT_FAILURE | DT_INVALID_PARAM;

	if (!mesh->isValidPolyRef(startRef) || !mesh->isValidPolyRef(endRef))
		return DT_FAILURE | DT_INVALID_PARAM;

	unsigned short int startAreaId;
	unsigned char startArea;
	mesh->getPolyAreaId(startRef,&startAreaId);
	mesh->getPolyArea(startRef,&startArea);

	unsigned short int endAreaId;
	unsigned char endArea;
	mesh->getPolyAreaId(endRef,&endAreaId);
	mesh->getPolyArea(endRef,&endArea);

	if (startRef == endRef)
	{
		path[0].bcLineRef = startAreaId;
		path[0].interBuc = endRef;
		*pathCount = 1;
		return DT_SUCCESS;
	}

	int j=0;

	if (startArea==2)  
	{
		bcBucket* bcn=NULL;
		for (vg = vertices(bcg); vg.first != vg.second; ++vg.first)
		{
			bcn= bcv_map[*vg.first];
			if (bcn->bcBucRef!=startAreaId)
				continue;
			p=9999;
			bcCalcPath(mesh,bcn,mapDensity,parentRef,0,startPos,endPos,endRef,fam,expl,path,&j,p,t);
			*pathCount = 1;	
			return DT_SUCCESS;
			
		}
	}
	else
	{
		if (endArea==1)
		{
			if (startAreaId == endAreaId)
			{
				path[0].bcLineRef = startAreaId;
				path[0].interBuc = endRef;
				*pathCount = 1;
				return DT_SUCCESS;
			}
		}
	
		bcBucket* bcnL=NULL;
		bcBucket* bcnR=NULL;
		struct bcLine* bcl=NULL;
		bool exist=false;
		for (boost::tie(ei, ei_end) = edges(bcg); ei != ei_end; ++ei)
		{
			bcl=line_map[*ei];
			
			if (bcl->bcLineRef!=startAreaId)
				continue;
			
			exist=true;

			bcnL=bcl->bucketL;
			bcnR=bcl->bucketR;

			break;
		}

		if (exist==false)
			return DT_FAILURE;

		float distL=bcCalcdist(startPos,bcnL->centre,endPos);
		float distR=bcCalcdist(startPos,bcnR->centre,endPos);

		if (distL<distR)
		{
			path[0].bcLineRef = bcl->bcLineRef;
			path[0].interBuc = bcnL->bucketPoly[0];
			bcCalcPath(mesh,bcnL,mapDensity,parentRef,bcl->bcLineRef,startPos,endPos,endRef,fam,expl,path,&j,p,t);
		}
		else
		{
			path[0].bcLineRef = bcl->bcLineRef;
			path[0].interBuc = bcnR->bucketPoly[0];
			bcCalcPath(mesh,bcnR,mapDensity,parentRef,bcl->bcLineRef,startPos,endPos,endRef,fam,expl,path,&j,p,t);
		}
    }

	if (j==0)
	{
		return DT_FAILURE;
	}
	else
	{
		*pathCount = j+1;
		return DT_SUCCESS;
	}
}

dtStatus bcBucketMesh::bcChoosePath(dtNavMesh* mesh,BcGraph& bcg,std::map<int, int> mapDensity,dtPolyRef parentRef,dtPolyRef startRef, dtPolyRef endRef,float* startPos,
	                              float* endPos,struct bcLineStep* path, int* pathCount,unsigned short* t)
{
	boost::associative_property_map<bucket_property_t> bcv_map(m_bct_prop);
	boost::associative_property_map<line_property_t> line_map(m_ln_prop);

	std::pair<b_vertex_iterator_t, b_vertex_iterator_t> vg;
	b_edge_iterator_t ei, ei_end;

	bcBucket* bc=NULL;
	bool top=false;

	*pathCount = 0;
	*t=0;
	unsigned short p=0;

	if (!startRef || !endRef)
		return DT_FAILURE | DT_INVALID_PARAM;

	if (!mesh->isValidPolyRef(startRef) || !mesh->isValidPolyRef(endRef))
		return DT_FAILURE | DT_INVALID_PARAM;

	unsigned short int startAreaId;
	unsigned char startArea;
	mesh->getPolyAreaId(startRef,&startAreaId);
	mesh->getPolyArea(startRef,&startArea);

	unsigned short int endAreaId;
	unsigned char endArea;
	mesh->getPolyAreaId(endRef,&endAreaId);
	mesh->getPolyArea(endRef,&endArea);

	if (startRef == endRef)
	{
		path[0].bcLineRef = startAreaId;
		path[0].interBuc = endRef;
		*pathCount = 1;
		return DT_SUCCESS;
	}

	int j=0;
	if (startArea==2)  
	{
		bcBucket* bcn=NULL;
		for (vg = vertices(bcg); vg.first != vg.second; ++vg.first)
		{
			bcn= bcv_map[*vg.first];
			if (bcn->bcBucRef!=startAreaId)
				continue;
			p=9999;
			bcCalcChoice(mesh,bcn,mapDensity,parentRef,startPos,endPos,endRef,path,&j,p,t);
			*pathCount = 1;	
			return DT_SUCCESS;
			
		}
	}
	else
	{
		if (endArea==1)
		{
			if (startAreaId == endAreaId)
			{
				path[0].bcLineRef = startAreaId;
				path[0].interBuc = endRef;
				*pathCount = 1;
				return DT_SUCCESS;
			}
		}
	
		bcBucket* bcnL=NULL;
		bcBucket* bcnR=NULL;
		struct bcLine* bcl=NULL;
		bool exist=false;
		for (boost::tie(ei, ei_end) = edges(bcg); ei != ei_end; ++ei)
		{
			bcl=line_map[*ei];
			
			if (bcl->bcLineRef!=startAreaId)
				continue;
			
			exist=true;

			bcnL=bcl->bucketL;
			bcnR=bcl->bucketR;

			break;
		}

		if (exist==false)
			return DT_FAILURE;

		float distL=bcCalcdist(startPos,bcnL->centre,endPos);
		float distR=bcCalcdist(startPos,bcnR->centre,endPos);

		if (distL<distR)
		{
			path[0].bcLineRef = bcl->bcLineRef;
			path[0].interBuc = bcnL->bucketPoly[0];
			bcCalcChoice(mesh,bcnL,mapDensity,parentRef,startPos,endPos,endRef,path,&j,p,t);
		}
		else
		{
			path[0].bcLineRef = bcl->bcLineRef;
			path[0].interBuc = bcnR->bucketPoly[0];
			bcCalcChoice(mesh,bcnR,mapDensity,parentRef,startPos,endPos,endRef,path,&j,p,t);
		}
    }

	if (j==0)
	{
		return DT_FAILURE;
	}
	else
	{
		*pathCount = j+1;
		return DT_SUCCESS;
	}
}

BcGraph& bcBucketMesh::getBucGraph()
{
	return m_bucketGraph;
}

bucket_property_t& bcBucketMesh::getBc_prop ()
{
	return m_bct_prop;
}

bucket_property_t& bcBucketMesh::getVisiBc_prop ()
{
	return m_bctv_prop;
}

line_property_t& bcBucketMesh::getLine_prop ()
{
	return m_ln_prop;
}

bool bcBucketMesh::bcBuildRoadClasses(BcGraph& bcg, double pts[], int n)
{
	boost::associative_property_map<line_property_t> line_map(m_ln_prop);
	b_edge_iterator_t ei, ei_end;
	struct bcLine* bcl=NULL;

	
	int total_points=boost::num_edges(bcg);
	vector<double> values;
	vector<Point> points;
	for (boost::tie(ei, ei_end) = edges(bcg); ei != ei_end; ++ei)
	{
		double value; int i;
		bcl=line_map[*ei];
		i=bcl->bcLineRef;
		value= bcl->betweeness;
		values.push_back(value);

		Point p(i, values);
		points.push_back(p);
		values.clear();
		
	}

    int max_iterations=100;
	KMeans kmeans(n, total_points, 1, max_iterations);
	int iter= kmeans.run(points,pts,3);

	if (iter==0)
		return false;

	for (boost::tie(ei, ei_end) = edges(bcg); ei != ei_end; ++ei)
	{
		bcl=line_map[*ei];
		
		if (bcl->betweeness >= pts[1])
		{
			bcl->interval=3;
		}
		if ((bcl->betweeness >= pts[0])&&(bcl->betweeness < pts[1]))
		{
			bcl->interval=2;
		}
		if (bcl->betweeness < pts[0])
		{
			bcl->interval=1;
		}
	}

	return true;
}

void bcBucketMesh::setBucGraph(BcGraph& bcg)
{
	m_bucketGraph=bcg;
}

void bcBucketMesh::setVisiBucGraph(BcGraph& bcg)
{
	m_visibilityBucketGraph=bcg;
}