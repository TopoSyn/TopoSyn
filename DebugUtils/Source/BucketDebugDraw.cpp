//

//

//#include <math.h>
#include "DebugDraw.h"
#include "BucketDebugDraw.h"
#include "DetourNavMesh.h"
#include "DetourCommon.h"
#include "DetourDebugDraw.h"
#include <vector>




static float bcdistancePtLine2d(const float* pt, const float* p, const float* q)
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

static void drawPolyBoundary(duDebugDraw* dd, const dtMeshTile* tile,const dtNavMesh& mesh,
							   const unsigned int col, const float linew,int i,int j)
{
dd->begin(DU_DRAW_LINES, linew);
unsigned int c = col;
const dtPoly* p = &tile->polys[i];
const dtPolyDetail* pd = &tile->detailMeshes[i];
int nj = (int)p->vertCount;
const float* v0 = &tile->verts[p->verts[j]*3];
const float* v1 = &tile->verts[p->verts[(j+1) % nj]*3];
dd->vertex(v0, c);
dd->vertex(v1, c);
dd->end();
}

void bcDebugDrawIsovist(duDebugDraw* dd,float* verts, int nverts,float h)
{
	if (!dd) return;
	const float linew=5.0f;
	unsigned int c = duRGBA(128,25,0,192);
	int i=0;

	dd->begin(DU_DRAW_LINES, linew);
	while ( i<nverts*4)
		{
			PointType vert0(verts[i], verts[i+1]);
			dd->vertex(vert0.x(),h,vert0.y(), c);
			i=i+2;
		}
	dd->end();

}

void bcDebugDrawBucketGraph(duDebugDraw* dd,bcBucketMesh* bucMesh,BcGraph& g,float h)
{
	if (!dd) return;

	const float linew=2.5f;
	unsigned int c = duRGBA(255,255,255,64);
	
	bucket_property_t& bc_prop=bucMesh->getBc_prop();

	std::vector<PointType> samples;

	b_edge_iterator_t ei, ei_end;
	
	for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
	{
		 struct bcBucket* v0=(bc_prop)[source(*ei, g)]; 
         struct bcBucket* v1=(bc_prop)[target(*ei, g)];
		
		 PointType vertex0(v0->centre[0], v0->centre[2]);
	     samples.push_back(vertex0);
	     PointType vertex1(v1->centre[0],v1->centre[2]);
	     samples.push_back(vertex1);
		 
	}

	dd->begin(DU_DRAW_LINES, linew);
	for (std::size_t i = 0; i < samples.size(); ++i) 
	{
         PointType vert = deconvolve(samples[i], shift_);
         dd->vertex(vert.x(),h,vert.y(), c);
	}
	dd->end();
	
}

static unsigned int ColorBetweenness(float f, double max)
{
    unsigned int col;
	float r,g,b,a;
	double lmax= max+10;
	int lmin=(int)floor(lmax/4);
	
	float fa=(lmax-f)/lmin;	
	int X=(int)floor(fa);	
	float Y=floor(255*(fa-X)); 

	switch(X)
	{
		case 0: r=255;g=Y;b=0;break;
		case 1: r=255-Y;g=255;b=0;break;
		case 2: r=0;g=255;b=Y;break;
		case 3: r=0;g=255-Y;b=255;break;
		case 4: r=0;g=0;b=255;break;
	}

	a=192;
	col = duRGBA(r,g,b,a);
	return col;
}

void bcDebugDrawBucketGraphColored(duDebugDraw* dd,bcBucketMesh* bucMesh,BcGraph& g,float h, double max)
{
	BcGraph& mcg = bucMesh->getBucGraph();
	bucket_property_t& bc_prop=bucMesh->getBc_prop();
	line_property_t& ln_prop=bucMesh->getLine_prop();
	b_edge_iterator_t ei, ei_end;

	dd->begin(DU_DRAW_LINES, 4.0);
	for (boost::tie(ei, ei_end) = edges(mcg); ei != ei_end; ++ei)
	{
		struct bcLine* L=(ln_prop)[*ei];
		float f = (float)L->betweeness;
		if (f>17)
		{
			unsigned int col = ColorBetweenness(f,max);
			float v0[3],v1[3];
			dtVcopy(v0,L->bucketL->centre);
			dtVcopy(v1,L->bucketR->centre);
			dd->vertex(v0, col);
			dd->vertex(v1, col);
		}
	}
	dd->end();

}

static float distancePtLine2d(const float* pt, const float* p, const float* q)
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

static float normalizeValue(float min, float max, float b)
{
	// Normalize penaly range.
	float minPen = min;
	float maxPen = max;
	
	const float penRange = maxPen-minPen;
	const float s = penRange > 0.001f ? (1.0f / penRange) : 1;
	
	float f = dtClamp((b-minPen)*s, 0.0f, 1.0f);

	return f;
}

static void drawPolyBoundaries(duDebugDraw* dd, const dtMeshTile* tile,
							   const unsigned int col, const float linew,
							   bool inner)
{
	static const float thr = 0.01f*0.01f;

	dd->begin(DU_DRAW_LINES, linew);

	for (int i = 0; i < tile->header->polyCount; ++i)
	{
		const dtPoly* p = &tile->polys[i];
		
		if (p->getType() == DT_POLYTYPE_OFFMESH_CONNECTION) continue;
		
		const dtPolyDetail* pd = &tile->detailMeshes[i];
		
		for (int j = 0, nj = (int)p->vertCount; j < nj; ++j)
		{
			unsigned int c = col;
			if (inner)
			{
				if (p->neis[j] == 0) continue;
				if (p->neis[j] & DT_EXT_LINK)
				{
					bool con = false;
					for (unsigned int k = p->firstLink; k != DT_NULL_LINK; k = tile->links[k].next)
					{
						if (tile->links[k].edge == j)
						{
							con = true;
							break;
						}
					}
					if (con)
						c = duRGBA(255,255,255,0);
					else
						c = duRGBA(0,0,0,48);
				}
				else
					c = duRGBA(0,48,64,32);
			}
			else
			{
				if (p->neis[j] != 0) continue;
			}
			
			const float* v0 = &tile->verts[p->verts[j]*3];
			const float* v1 = &tile->verts[p->verts[(j+1) % nj]*3];
			
			// Draw detail mesh edges which align with the actual poly edge.
			// This is really slow.
			for (int k = 0; k < pd->triCount; ++k)
			{
				const unsigned char* t = &tile->detailTris[(pd->triBase+k)*4];
				const float* tv[3];
				for (int m = 0; m < 3; ++m)
				{
					if (t[m] < p->vertCount)
						tv[m] = &tile->verts[p->verts[t[m]]*3];
					else
						tv[m] = &tile->detailVerts[(pd->vertBase+(t[m]-p->vertCount))*3];
				}
				for (int m = 0, n = 2; m < 3; n=m++)
				{
					if (((t[3] >> (n*2)) & 0x3) == 0) continue;	// Skip inner detail edges.
					if (distancePtLine2d(tv[n],v0,v1) < thr &&
						distancePtLine2d(tv[m],v0,v1) < thr)
					{
						dd->vertex(tv[n], c);
						dd->vertex(tv[m], c);
					}
				}
			}
		}
	}
	dd->end();
}

static void drawMeshTile(duDebugDraw* dd, const dtNavMesh& mesh, const dtNavMeshQuery* query,
						 const dtMeshTile* tile, unsigned char flags)
{
	dtPolyRef base = mesh.getPolyRefBase(tile);

	int tileNum = mesh.decodePolyIdTile(base);
	
	dd->depthMask(false);

	dd->begin(DU_DRAW_TRIS);
	for (int i = 0; i < tile->header->polyCount; ++i)
	{
		const dtPoly* p = &tile->polys[i];
		if (p->getType() == DT_POLYTYPE_OFFMESH_CONNECTION)	// Skip off-mesh links.
			continue;
			
		const dtPolyDetail* pd = &tile->detailMeshes[i];

		unsigned int col;
		//if (query && query->isInClosedList(base | (dtPolyRef)i))
		//	col = duRGBA(0,192,255,64);//255,196,0,64
		//else
		//{
			if (flags & DU_DRAWNAVMESH_COLOR_TILES)
			{
				col = duIntToCol(tileNum, 128);
			}
			else
			{
				
				if (p->getArea() == 0) // Treat zero area type as default.
					col = duRGBA(0,192,255,64);
				/*else
					col = duIntToCol(p->getArea(), 64);*/
				if (p->getArea() == 1) 
					col = duRGBA(255,196,0,64);
				if (p->getArea() == 2) 
					col = duRGBA(0,0,64,64);
				if (p->getArea() == 3) 
					col = duRGBA(64,64,64,64);
			}
		//}
		
		for (int j = 0; j < pd->triCount; ++j)
		{
			const unsigned char* t = &tile->detailTris[(pd->triBase+j)*4];
			for (int k = 0; k < 3; ++k)
			{
				if (t[k] < p->vertCount)
					dd->vertex(&tile->verts[p->verts[t[k]]*3], col);
				else
					dd->vertex(&tile->detailVerts[(pd->vertBase+t[k]-p->vertCount)*3], col);
			}
		}
	}
	dd->end();
	
	// Draw inter poly boundaries
	drawPolyBoundaries(dd, tile, duRGBA(0,48,64,32), 1.5f, true);
	
	// Draw outer poly boundaries
	drawPolyBoundaries(dd, tile, duRGBA(0,48,64,220), 2.5f, false);

	if (flags & DU_DRAWNAVMESH_OFFMESHCONS)
	{
		dd->begin(DU_DRAW_LINES, 2.0f);
		for (int i = 0; i < tile->header->polyCount; ++i)
		{
			const dtPoly* p = &tile->polys[i];
			if (p->getType() != DT_POLYTYPE_OFFMESH_CONNECTION)	// Skip regular polys.
				continue;
			
			unsigned int col, col2;
			if (query && query->isInClosedList(base | (dtPolyRef)i))
				col = duRGBA(255,196,0,220);
			else
				col = duDarkenCol(duIntToCol(p->getArea(), 220));
			
			const dtOffMeshConnection* con = &tile->offMeshCons[i - tile->header->offMeshBase];
			const float* va = &tile->verts[p->verts[0]*3];
			const float* vb = &tile->verts[p->verts[1]*3];

			// Check to see if start and end end-points have links.
			bool startSet = false;
			bool endSet = false;
			for (unsigned int k = p->firstLink; k != DT_NULL_LINK; k = tile->links[k].next)
			{
				if (tile->links[k].edge == 0)
					startSet = true;
				if (tile->links[k].edge == 1)
					endSet = true;
			}
			
			// End points and their on-mesh locations.
			dd->vertex(va[0],va[1],va[2], col);
			dd->vertex(con->pos[0],con->pos[1],con->pos[2], col);
			col2 = startSet ? col : duRGBA(220,32,16,196);
			duAppendCircle(dd, con->pos[0],con->pos[1]+0.1f,con->pos[2], con->rad, col2);

			dd->vertex(vb[0],vb[1],vb[2], col);
			dd->vertex(con->pos[3],con->pos[4],con->pos[5], col);
			col2 = endSet ? col : duRGBA(220,32,16,196);
			duAppendCircle(dd, con->pos[3],con->pos[4]+0.1f,con->pos[5], con->rad, col2);
			
			// End point vertices.
			dd->vertex(con->pos[0],con->pos[1],con->pos[2], duRGBA(0,48,64,196));
			dd->vertex(con->pos[0],con->pos[1]+0.2f,con->pos[2], duRGBA(0,48,64,196));
			
			dd->vertex(con->pos[3],con->pos[4],con->pos[5], duRGBA(0,48,64,196));
			dd->vertex(con->pos[3],con->pos[4]+0.2f,con->pos[5], duRGBA(0,48,64,196));
			
			// Connection arc.
			duAppendArc(dd, con->pos[0],con->pos[1],con->pos[2], con->pos[3],con->pos[4],con->pos[5], 0.25f,
						(con->flags & 1) ? 0.6f : 0, 0.6f, col);
		}
		dd->end();
	}
	
	const unsigned int vcol = duRGBA(0,0,0,196);
	dd->begin(DU_DRAW_POINTS, 3.0f);
	for (int i = 0; i < tile->header->vertCount; ++i)
	{
		const float* v = &tile->verts[i*3];
		dd->vertex(v[0], v[1], v[2], vcol);
	}
	dd->end();

	dd->depthMask(true);
}

void bcDebugDrawBucketMesh(struct duDebugDraw* dd, const dtNavMesh& mesh, const dtNavMeshQuery& query, unsigned char flags)
{
	if (!dd) return;

	const dtNavMeshQuery* q = (flags & DU_DRAWNAVMESH_CLOSEDLIST) ? &query : 0;
	
	for (int i = 0; i < mesh.getMaxTiles(); ++i)
	{
		const dtMeshTile* tile = mesh.getTile(i);
		if (!tile->header) continue;
		drawMeshTile(dd, mesh, q, tile, flags);
	}
}







