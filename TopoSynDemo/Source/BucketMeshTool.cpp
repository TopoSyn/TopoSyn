//

//

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "SDL.h"
#include "SDL_opengl.h"
#include "imgui.h"
#include "BucketMeshTool.h"
#include "Sample.h"
#include "Recast.h"
#include "RecastDebugDraw.h"
#include "DetourNavMesh.h"
#include "DetourNavMeshBuilder.h"
#include "DetourDebugDraw.h"
#include "DetourCommon.h"
#include "VoroDiag.h"
#include "BucketBuilder.h"
#include "BucketDebugDraw.h"

#ifdef WIN32
#	define snprintf _snprintf
#endif

// Uncomment this to dump all the requests in stdout.
#define DUMP_REQS

inline bool inRange(const float* v1, const float* v2, const float r, const float h)
{
	const float dx = v2[0] - v1[0];
	const float dy = v2[1] - v1[1];
	const float dz = v2[2] - v1[2];
	return (dx*dx + dz*dz) < r*r && fabsf(dy) < h;
}


static int fixupCorridor(dtPolyRef* path, const int npath, const int maxPath,
						 const dtPolyRef* visited, const int nvisited)
{
	int furthestPath = -1;
	int furthestVisited = -1;
	
	// Find furthest common polygon.
	for (int i = npath-1; i >= 0; --i)
	{
		bool found = false;
		for (int j = nvisited-1; j >= 0; --j)
		{
			if (path[i] == visited[j])
			{
				furthestPath = i;
				furthestVisited = j;
				found = true;
			}
		}
		if (found)
			break;
	}

	// If no intersection found just return current path. 
	if (furthestPath == -1 || furthestVisited == -1)
		return npath;
	
	// Concatenate paths.	

	// Adjust beginning of the buffer to include the visited.
	const int req = nvisited - furthestVisited;
	const int orig = rcMin(furthestPath+1, npath);
	int size = rcMax(0, npath-orig);
	if (req+size > maxPath)
		size = maxPath-req;
	if (size)
		memmove(path+req, path+orig, size*sizeof(dtPolyRef));
	
	// Store visited
	for (int i = 0; i < req; ++i)
		path[i] = visited[(nvisited-1)-i];				
	
	return req+size;
}

// This function checks if the path has a small U-turn, that is,
// a polygon further in the path is adjacent to the first polygon
// in the path. If that happens, a shortcut is taken.
// This can happen if the target (T) location is at tile boundary,
// and we're (S) approaching it parallel to the tile edge.
// The choice at the vertex can be arbitrary, 
//  +---+---+
//  |:::|:::|
//  +-S-+-T-+
//  |:::|   | <-- the step can end up in here, resulting U-turn path.
//  +---+---+
static int fixupShortcuts(dtPolyRef* path, int npath, dtNavMeshQuery* navQuery)
{
	if (npath < 3)
		return npath;

	// Get connected polygons
	static const int maxNeis = 16;
	dtPolyRef neis[maxNeis];
	int nneis = 0;

	const dtMeshTile* tile = 0;
	const dtPoly* poly = 0;
	if (dtStatusFailed(navQuery->getAttachedNavMesh()->getTileAndPolyByRef(path[0], &tile, &poly)))
		return npath;
	
	for (unsigned int k = poly->firstLink; k != DT_NULL_LINK; k = tile->links[k].next)
	{
		const dtLink* link = &tile->links[k];
		if (link->ref != 0)
		{
			if (nneis < maxNeis)
				neis[nneis++] = link->ref;
		}
	}

	// If any of the neighbour polygons is within the next few polygons
	// in the path, short cut to that polygon directly.
	static const int maxLookAhead = 6;
	int cut = 0;
	for (int i = dtMin(maxLookAhead, npath) - 1; i > 1 && cut == 0; i--) {
		for (int j = 0; j < nneis; j++)
		{
			if (path[i] == neis[j]) {
				cut = i;
				break;
			}
		}
	}
	if (cut > 1)
	{
		int offset = cut-1;
		npath -= offset;
		for (int i = 1; i < npath; i++)
			path[i] = path[i+offset];
	}

	return npath;
}

static bool getSteerTarget(dtNavMeshQuery* navQuery, const float* startPos, const float* endPos,
						   const float minTargetDist,
						   const dtPolyRef* path, const int pathSize,
						   float* steerPos, unsigned char& steerPosFlag, dtPolyRef& steerPosRef,
						   float* outPoints = 0, int* outPointCount = 0)							 
{
	// Find steer target.
	static const int MAX_STEER_POINTS = 3;
	float steerPath[MAX_STEER_POINTS*3];
	unsigned char steerPathFlags[MAX_STEER_POINTS];
	dtPolyRef steerPathPolys[MAX_STEER_POINTS];
	int nsteerPath = 0;
	navQuery->findStraightPath(startPos, endPos, path, pathSize,
							   steerPath, steerPathFlags, steerPathPolys, &nsteerPath, MAX_STEER_POINTS);
	if (!nsteerPath)
		return false;
		
	if (outPoints && outPointCount)
	{
		*outPointCount = nsteerPath;
		for (int i = 0; i < nsteerPath; ++i)
			dtVcopy(&outPoints[i*3], &steerPath[i*3]);
	}

	
	// Find vertex far enough to steer to.
	int ns = 0;
	while (ns < nsteerPath)
	{
		// Stop at Off-Mesh link or when point is further than slop away.
		if ((steerPathFlags[ns] & DT_STRAIGHTPATH_OFFMESH_CONNECTION) ||
			!inRange(&steerPath[ns*3], startPos, minTargetDist, 1000.0f))
			break;
		ns++;
	}
	// Failed to find good point to steer to.
	if (ns >= nsteerPath)
		return false;
	
	dtVcopy(steerPos, &steerPath[ns*3]);
	steerPos[1] = startPos[1];
	steerPosFlag = steerPathFlags[ns];
	steerPosRef = steerPathPolys[ns];
	
	return true;
}

BucketMeshTool::BucketMeshTool() :
	m_sample(0),
	m_navMesh(0),
	m_bucMesh(0),
	m_navQuery(0),
	m_toolMode(TOOLMODE_PATHFIND),
	m_PathFindStrategy(0),
	m_queryPointRef(0),
	m_dirRef(0),
	m_nverts(-1),
	m_hitResult(false),
	m_qposSet(false),
	m_dposSet(false),
	m_Agent_FamSet(false)
{
	m_filter.setIncludeFlags(SAMPLE_POLYFLAGS_ALL ^ SAMPLE_POLYFLAGS_DISABLED);
	m_filter.setExcludeFlags(0);

	m_polyPickExt[0] = 2;
	m_polyPickExt[1] = 4;
	m_polyPickExt[2] = 2;

	m_toolParams.m_expandDebugDraw = true;
	m_toolParams.m_showConnectivity = false;
	m_toolParams.m_showIntegration = false;
	m_toolParams.m_showDetailAll = false;
	m_toolParams.m_showLabels = false;
	m_toolParams.m_showNodes = false;
	m_toolParams.m_showPath = false;
	m_toolParams.m_showGraph = false;
	
}

BucketMeshTool::~BucketMeshTool()
{
	if (m_sample)
	{
		unsigned char flags = DU_DRAWNAVMESH_CLOSEDLIST;
		if (m_navMesh)
			flags |= DU_DRAWNAVMESH_OFFMESHCONS;
		m_sample->setNavMeshDrawFlags(flags);
	}
	m_vg.clear();
}

void BucketMeshTool::init(Sample* sample)
{
	m_sample = sample;
	m_navMesh = sample->getNavMesh();
	m_bucMesh = sample->getBucketMesh();
	m_navQuery = sample->getNavMeshQuery();
	
	recalc();
	
}

void BucketMeshTool::handleMenu()
{
	imguiIndent();

	if (imguiCheck("Pathfinding Algo", m_toolMode == TOOLMODE_PATHFIND))
	{
		m_toolMode = TOOLMODE_PATHFIND;
		recalc();
	}
	if (m_toolMode == TOOLMODE_PATHFIND)
	{
		imguiIndent();
			imguiLabel("Familiarity Level");
			if (imguiCheck("High", m_Agent_Fam == 1))
			{
				m_Agent_Fam = 1;
				m_Agent_FamSet = true;
				recalc();
			}
			if (imguiCheck("Medium", m_Agent_Fam == 2))
			{
				m_Agent_Fam = 2;
				m_Agent_FamSet = true;
				recalc();
			}
			if (imguiCheck("Low", m_Agent_Fam == 3))
			{
				m_Agent_Fam = 3;
				m_Agent_FamSet = true;
				recalc();
			}

			imguiLabel("Pathfinding strategy");
			if (imguiCheck("Min distance", m_PathFindStrategy == 0))
			{
				m_PathFindStrategy = DT_FINDPATH_MINIMAL_DISTANCE;
				recalc();
			}
			if (imguiCheck("Proposed Algo", m_PathFindStrategy == DT_FINDPATH_MINIMAL_DEVIATION))
			{
				m_PathFindStrategy = DT_FINDPATH_MINIMAL_DEVIATION;
				recalc();
			}
		imguiUnindent();
	}
}

void BucketMeshTool::handleClick(const float* /*s*/, const float* p, bool shift)
{
	if (shift)
	{
		m_dposSet = true;
		rcVcopy(m_dpos, p);
	}
	else
	{
		m_qposSet = true;
		rcVcopy(m_qpos, p);
	}

	recalc();
}

void BucketMeshTool::handleStep()
{
}

void BucketMeshTool::handleToggle()
{

}

void BucketMeshTool::handleUpdate(const float /*dt*/)
{
	
}

void BucketMeshTool::reset()
{
	m_queryPointRef = 0;
	m_dirRef = 0;
	m_npolys = 0;
	m_Agent_FamSet = false;
	memset(m_hitPosMax, 0, sizeof(m_hitPosMax));
	memset(m_hitPosMin, 0, sizeof(m_hitPosMin));
	memset(m_hitNormal, 0, sizeof(m_hitNormal));
	m_vg.clear();
}

static float bcDistancePtPtSqr2D(float* A, float* B)
{
	float dx = (B[0]-A[0]);
	float dz = (B[2]-A[2]);
	return dx*dx + dz*dz;
}

static float calcDist(float* path, int pathCount)
{
	float dist=0.0;
	float prec_pos[3];
	dtVcopy(prec_pos,&path[0]);
	for (int i=0 ; i < pathCount ; i++)
	{
		dist=dist+bcDistancePtPtSqr2D(prec_pos,&path[i*3]);
		dtVcopy(prec_pos, &path[i*3]);
	}
	return dist;
}

void BucketMeshTool::recalc()
{
	if (!m_navMesh)
		return;

	if (!m_bucMesh)
		return;

	BcGraph& bucg = m_bucMesh->getBucGraph();

	if (m_qposSet)
		m_navQuery->findNearestPoly(m_qpos, m_polyPickExt, &m_filter, &m_queryPointRef, 0);
	else
		m_queryPointRef = 0;

	if (m_dposSet)
		m_navQuery->findNearestPoly(m_dpos, m_polyPickExt, &m_filter, &m_dirRef, 0);
	else
		m_dirRef = 0;

	m_pathFindStatus = DT_FAILURE;
	
	
	if (m_toolMode == TOOLMODE_ISOVIST_POLY)
	{
		if (m_qposSet && m_queryPointRef)
		{
			m_nverts=1024;
			m_verts=(float*)bcAlloc(sizeof(float)*m_nverts*4, BC_ALLOC_PERM);
			//m_bucMesh->bcQueryPointVisibility(m_qpos,m_verts,&m_nverts);
		}
	}
	else if(m_toolMode == TOOLMODE_PATHFIND)
	{
#ifdef DUMP_REQS
			printf("TOOLMODE_PATHFIND  \n"); 
#endif

		if (m_qposSet && m_dposSet && m_queryPointRef && m_dirRef)
		{
			m_startRef = m_queryPointRef;
			m_endRef = m_dirRef;
			dtVcopy(m_spos,m_qpos);
            dtVcopy(m_epos,m_dpos);

			if (m_PathFindStrategy == DT_FINDPATH_MINIMAL_DISTANCE)
			{
				m_navQuery->findPath(m_startRef, m_endRef, m_spos, m_epos, &m_filter, m_polys, &m_npolys, MAX_POLYS);
			}
			else if ((m_PathFindStrategy == DT_FINDPATH_MINIMAL_DEVIATION)&&(m_Agent_FamSet)) 
			{
				struct bcLineStep path[128];
				int pathCount=0;
				unsigned short nturns=0;
	
				dtStatus status = m_bucMesh->bcFindPathLight(m_navMesh,bucg,m_startRef,m_startRef,m_endRef,m_spos,m_epos,m_Agent_Fam,0,path,&pathCount,&nturns);
			
				if (status==DT_SUCCESS)
				{
					int j=0;
					for (int i=0;i<pathCount;++i)
					{
						m_endRef=path[i].interBuc;
						m_navQuery->getPolyCenter(m_navMesh,m_endRef,m_epos);
					
						m_navQuery->findPath(m_startRef,m_endRef,m_spos,m_epos,&m_filter,&m_polys[j],&m_npolys,MAX_POLYS);
						j=j+m_npolys;
						m_startRef =m_endRef;
						dtVcopy(m_spos,m_epos);
					}
					m_npolys=j;
				}
			}
			
			m_nsmoothPath = 0;

			if (m_npolys)
			{
				// Iterate over the path to find smooth path on the detail mesh surface.
				dtPolyRef polys[MAX_POLYS];
				memcpy(polys, m_polys, sizeof(dtPolyRef)*m_npolys); 
				int npolys = m_npolys;
				
				float iterPos[3], targetPos[3];
				m_navQuery->closestPointOnPoly(m_queryPointRef, m_qpos, iterPos, 0);
				m_navQuery->closestPointOnPoly(m_dirRef, m_dpos, targetPos, 0);
				
				static const float STEP_SIZE = 0.5f;
				static const float SLOP = 0.01f;
				
				m_nsmoothPath = 0;
				
				dtVcopy(&m_smoothPath[m_nsmoothPath*3], iterPos);
				m_nsmoothPath++;
				
				// Move towards target a small advancement at a time until target reached or
				// when ran out of memory to store the path.
				while (npolys && m_nsmoothPath < MAX_SMOOTH)
				{
					// Find location to steer towards.
					float steerPos[3];
					unsigned char steerPosFlag;
					dtPolyRef steerPosRef;
					
					if (!getSteerTarget(m_navQuery, iterPos, targetPos, SLOP,
										polys, npolys, steerPos, steerPosFlag, steerPosRef))
						break;
					
					bool endOfPath = (steerPosFlag & DT_STRAIGHTPATH_END) ? true : false;
					bool offMeshConnection = (steerPosFlag & DT_STRAIGHTPATH_OFFMESH_CONNECTION) ? true : false;
					
					// Find movement delta.
					float delta[3], len;
					dtVsub(delta, steerPos, iterPos);
					len = dtMathSqrtf(dtVdot(delta, delta));
					// If the steer target is end of path or off-mesh link, do not move past the location.
					if ((endOfPath || offMeshConnection) && len < STEP_SIZE)
						len = 1;
					else
						len = STEP_SIZE / len;
					float moveTgt[3];
					dtVmad(moveTgt, iterPos, delta, len);
					
					// Move
					float result[3];
					dtPolyRef visited[16];
					int nvisited = 0;
					m_navQuery->moveAlongSurface(polys[0], iterPos, moveTgt, &m_filter,
												 result, visited, &nvisited, 16);

					npolys = fixupCorridor(polys, npolys, MAX_POLYS, visited, nvisited);
					npolys = fixupShortcuts(polys, npolys, m_navQuery);

					float h = 0;
					m_navQuery->getPolyHeight(polys[0], result, &h);
					result[1] = h;
					dtVcopy(iterPos, result);

					// Handle end of path and off-mesh links when close enough.
					if (endOfPath && inRange(iterPos, steerPos, SLOP, 1.0f))
					{
						// Reached end of path.
						dtVcopy(iterPos, targetPos);
						if (m_nsmoothPath < MAX_SMOOTH)
						{
							dtVcopy(&m_smoothPath[m_nsmoothPath*3], iterPos);
							m_nsmoothPath++;
						}
						break;
					}
					else if (offMeshConnection && inRange(iterPos, steerPos, SLOP, 1.0f))
					{
						// Reached off-mesh connection.
						float startPos[3], endPos[3];
						
						// Advance the path up to and over the off-mesh connection.
						dtPolyRef prevRef = 0, polyRef = polys[0];
						int npos = 0;
						while (npos < npolys && polyRef != steerPosRef)
						{
							prevRef = polyRef;
							polyRef = polys[npos];
							npos++;
						}
						for (int i = npos; i < npolys; ++i)
							polys[i-npos] = polys[i];
						npolys -= npos;
						
						// Handle the connection.
						dtStatus status = m_navMesh->getOffMeshConnectionPolyEndPoints(prevRef, polyRef, startPos, endPos);
						if (dtStatusSucceed(status))
						{
							if (m_nsmoothPath < MAX_SMOOTH)
							{
								dtVcopy(&m_smoothPath[m_nsmoothPath*3], startPos);
								m_nsmoothPath++;
								// Hack to make the dotted path not visible during off-mesh connection.
								if (m_nsmoothPath & 1)
								{
									dtVcopy(&m_smoothPath[m_nsmoothPath*3], startPos);
									m_nsmoothPath++;
								}
							}
							// Move position at the other side of the off-mesh link.
							dtVcopy(iterPos, endPos);
							float eh = 0.0f;
							m_navQuery->getPolyHeight(polys[0], iterPos, &eh);
							iterPos[1] = eh;
						}
					}
					
					// Store results.
					if (m_nsmoothPath < MAX_SMOOTH)
					{
						dtVcopy(&m_smoothPath[m_nsmoothPath*3], iterPos);
						m_nsmoothPath++;
						
					}

				}

				float d=calcDist(m_smoothPath, m_nsmoothPath);
			}

		}
		else
		{
			m_npolys = 0;
			m_nsmoothPath = 0;
		}
	}
		
}

void BucketMeshTool::handleRender()
{
	DebugDrawGL dd;
	
	static const unsigned int queryCol = duRGBA(128,25,0,192);
	static const unsigned int dirCol = duRGBA(51,102,0,129);
	static const unsigned int isoCol = duRGBA(0,0,0,64);

	
	static const unsigned int pathCol = duRGBA(0,0,0,64);
	
	const float agentRadius = m_sample->getAgentRadius();
	const float agentHeight = m_sample->getAgentHeight();
	const float agentClimb = m_sample->getAgentClimb();
	
	dd.depthMask(false);
	if (m_qposSet)
		drawQueryPoint(m_qpos,agentRadius, agentHeight, agentClimb,queryCol);
	if (m_dposSet)
		drawQueryPoint(m_dpos,agentRadius, agentHeight, agentClimb,dirCol);
	dd.depthMask(true);
	
	if (!m_navMesh)
	{
		return;
	}

	if (m_toolMode == TOOLMODE_ISOVIST_POLY)
	{
		if (m_nverts)
		{
			bcDebugDrawIsovist(&dd,m_verts,m_nverts,m_qpos[1]);
		}
	}
	else if (m_toolMode == TOOLMODE_PATHFIND)
	{
		duDebugDrawNavMeshPoly(&dd, *m_navMesh, m_queryPointRef, queryCol);
		duDebugDrawNavMeshPoly(&dd, *m_navMesh, m_dirRef, dirCol);
		
		if (m_npolys)
		{
			for (int i = 0; i < m_npolys; ++i)
			{
				if (m_polys[i] == m_queryPointRef || m_polys[i] == m_dirRef)
					continue;
				duDebugDrawNavMeshPoly(&dd, *m_navMesh, m_polys[i], pathCol);
			}
		}

		if (m_nsmoothPath)
		{
			dd.depthMask(false);
			const unsigned int spathCol = duRGBA(0,0,0,220);
			dd.begin(DU_DRAW_LINES, 3.0f);
			for (int i = 0; i < m_nsmoothPath; ++i)
				dd.vertex(m_smoothPath[i*3], m_smoothPath[i*3+1]+0.1f, m_smoothPath[i*3+2], spathCol);
			dd.end();
			dd.depthMask(true);
		}
	}
			
}

void BucketMeshTool::handleRenderOverlay(double* proj, double* model, int* view)
{
	GLdouble x, y, z;

	// Draw query point label
	if (m_qposSet && gluProject((GLdouble)m_qpos[0], (GLdouble)m_qpos[1], (GLdouble)m_qpos[2],
								model, proj, view, &x, &y, &z))
	{
		imguiDrawText((int)x, (int)(y-25), IMGUI_ALIGN_CENTER, "Start", imguiRGBA(0,0,0,220));
	}

	if (m_dposSet && gluProject((GLdouble)m_dpos[0], (GLdouble)m_dpos[1], (GLdouble)m_dpos[2],
								model, proj, view, &x, &y, &z))
	{
		imguiDrawText((int)x, (int)(y-25), IMGUI_ALIGN_CENTER, "End", imguiRGBA(0,0,0,220));
	}

	// Tool help
	const int h = view[3];
	int ty = h-40;

	if (m_toolMode == TOOLMODE_PATHFIND)
	{
		imguiDrawText(280, ty, IMGUI_ALIGN_LEFT, "LMB+SHIFT: Set start location  LMB: Set end location", imguiRGBA(255,255,255,192));	
	}
	
}

void BucketMeshTool::drawQueryPoint(const float* pos, float r, float h, float c,const unsigned int col)
{
	DebugDrawGL dd;
	
	glDepthMask(GL_FALSE);
	
	// Agent dimensions.	
	duDebugDrawCylinderWire(&dd, pos[0]-r, pos[1]+0.02f, pos[2]-r, pos[0]+r, pos[1]+h, pos[2]+r, col, 2.0f);

	duDebugDrawCircle(&dd, pos[0],pos[1]+c,pos[2],r,duRGBA(0,0,0,64),1.0f);

	unsigned int colb = duRGBA(0,0,0,196);
	dd.begin(DU_DRAW_LINES);
	dd.vertex(pos[0], pos[1]-c, pos[2], colb);
	dd.vertex(pos[0], pos[1]+c, pos[2], colb);
	dd.vertex(pos[0]-r/2, pos[1]+0.02f, pos[2], colb);
	dd.vertex(pos[0]+r/2, pos[1]+0.02f, pos[2], colb);
	dd.vertex(pos[0], pos[1]+0.02f, pos[2]-r/2, colb);
	dd.vertex(pos[0], pos[1]+0.02f, pos[2]+r/2, colb);
	dd.end();
	
	glDepthMask(GL_TRUE);
}



