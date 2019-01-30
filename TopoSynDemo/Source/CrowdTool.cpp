//
// Copyright (c) 2009-2010 Mikko Mononen memon@inside.org
//
// This software is provided 'as-is', without any express or implied
// warranty.  In no event will the authors be held liable for any damages
// arising from the use of this software.
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
// 1. The origin of this software must not be misrepresented; you must not
//    claim that you wrote the original software. If you use this software
//    in a product, an acknowledgment in the product documentation would be
//    appreciated but is not required.
// 2. Altered source versions must be plainly marked as such, and must not be
//    misrepresented as being the original software.
// 3. This notice may not be removed or altered from any source distribution.
//

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#include <float.h>
#include "SDL.h"
#include "SDL_opengl.h"
#include "imgui.h"
#include "CrowdTool.h"
#include "InputGeom.h"
#include "Sample.h"
#include "DetourCrowd.h"
#include "DetourDebugDraw.h"
#include "DetourObstacleAvoidance.h"
#include "DetourCommon.h"
#include "DetourNode.h"
#include "SampleInterfaces.h"

#ifdef WIN32
#	define snprintf _snprintf
#endif

// Returns a random number [0..1)
static float frand()
{
//	return ((float)(rand() & 0xffff)/(float)0xffff);
	return (float)rand()/(float)RAND_MAX;
}

static bool isectSegAABB(const float* sp, const float* sq,
						 const float* amin, const float* amax,
						 float& tmin, float& tmax)
{
	static const float EPS = 1e-6f;
	
	float d[3];
	dtVsub(d, sq, sp);
	tmin = 0;  // set to -FLT_MAX to get first hit on line
	tmax = FLT_MAX;		// set to max distance ray can travel (for segment)
	
	// For all three slabs
	for (int i = 0; i < 3; i++)
	{
		if (fabsf(d[i]) < EPS)
		{
			// Ray is parallel to slab. No hit if origin not within slab
			if (sp[i] < amin[i] || sp[i] > amax[i])
				return false;
		}
		else
		{
			// Compute intersection t value of ray with near and far plane of slab
			const float ood = 1.0f / d[i];
			float t1 = (amin[i] - sp[i]) * ood;
			float t2 = (amax[i] - sp[i]) * ood;
			// Make t1 be intersection with near plane, t2 with far plane
			if (t1 > t2) dtSwap(t1, t2);
			// Compute the intersection of slab intersections intervals
			if (t1 > tmin) tmin = t1;
			if (t2 < tmax) tmax = t2;
			// Exit with no collision as soon as slab intersection becomes empty
			if (tmin > tmax) return false;
		}
	}
	
	return true;
}

static void getAgentBounds(const dtCrowdAgent* ag, float* bmin, float* bmax)
{
	const float* p = ag->npos;
	const float r = ag->params.radius;
	const float h = ag->params.height;
	bmin[0] = p[0] - r;
	bmin[1] = p[1];
	bmin[2] = p[2] - r;
	bmax[0] = p[0] + r;
	bmax[1] = p[1] + h;
	bmax[2] = p[2] + r;
}

CrowdToolState::CrowdToolState() :
	m_sample(0),
	m_nav(0),
	m_crowd(0),
	m_bucMesh(0),
	m_targetRef(0),
	m_run(true)
{
	m_toolParams.m_expandSelectedDebugDraw = true;
	m_toolParams.m_showCorners = false;
	m_toolParams.m_showCollisionSegments = false;
	m_toolParams.m_showFrequentedRoads = false;
	m_toolParams.m_showBetweenness = false;
	m_toolParams.m_showGlobalIntegration = false;
	m_toolParams.m_showPath = false;
	m_toolParams.m_showVO = false;
	m_toolParams.m_showOpt = false;
	m_toolParams.m_showNeis = false;
	m_toolParams.m_expandDebugDraw = false;
	m_toolParams.m_showLabels = false;
	m_toolParams.m_showGrid = false;
	m_toolParams.m_showNodes = false;
	m_toolParams.m_showPerfGraph = false;
	m_toolParams.m_showDetailAll = false;
	m_toolParams.m_expandOptions = true;
	m_toolParams.m_anticipateTurns = true;
	m_toolParams.m_optimizeVis = true;
	m_toolParams.m_optimizeTopo = true;
	m_toolParams.m_obstacleAvoidance = true;
	m_toolParams.m_obstacleAvoidanceType = 3.0f;
	m_toolParams.m_separation = false;
	m_toolParams.m_separationWeight = 2.0f;
	
	//memset(m_trails, 0, sizeof(m_trails));
	
	m_vod = dtAllocObstacleAvoidanceDebugData();
	m_vod->init(2048);
	
	memset(&m_agentDebug, 0, sizeof(m_agentDebug));
	m_agentDebug.idx = -1;
	m_agentDebug.vod = m_vod;
}

CrowdToolState::~CrowdToolState()
{
	dtFreeObstacleAvoidanceDebugData(m_vod);
}

void CrowdToolState::init(class Sample* sample)
{
	if (m_sample != sample)
	{
		m_sample = sample;
//		m_oldFlags = m_sample->getNavMeshDrawFlags();
//		m_sample->setNavMeshDrawFlags(m_oldFlags & ~DU_DRAWNAVMESH_CLOSEDLIST);
	}
	
	dtNavMesh* nav = m_sample->getNavMesh();
	dtCrowd* crowd = m_sample->getCrowd();
	bcBucketMesh* bucMesh = m_sample->getBucketMesh();
	

	if (nav && crowd && bucMesh && (m_nav != nav || m_crowd != crowd))
	{
		m_nav = nav;
		m_crowd = crowd;
	    m_bucMesh = bucMesh;

		crowd->init(MAX_AGENTS, m_sample->getAgentRadius(), nav,bucMesh);
		
		// Make polygons with 'disabled' flag invalid.
		crowd->getEditableFilter(0)->setExcludeFlags(SAMPLE_POLYFLAGS_DISABLED);
		
		// Setup local avoidance params to different qualities.
		dtObstacleAvoidanceParams params;
		// Use mostly default settings, copy from dtCrowd.
		memcpy(&params, crowd->getObstacleAvoidanceParams(0), sizeof(dtObstacleAvoidanceParams));
		
		// Low (11)
		params.velBias = 0.5f;
		params.adaptiveDivs = 5;
		params.adaptiveRings = 2;
		params.adaptiveDepth = 1;
		crowd->setObstacleAvoidanceParams(0, &params);
		
		// Medium (22)
		params.velBias = 0.5f;
		params.adaptiveDivs = 5; 
		params.adaptiveRings = 2;
		params.adaptiveDepth = 2;
		crowd->setObstacleAvoidanceParams(1, &params);
		
		// Good (45)
		params.velBias = 0.5f;
		params.adaptiveDivs = 7;
		params.adaptiveRings = 2;
		params.adaptiveDepth = 3;
		crowd->setObstacleAvoidanceParams(2, &params);
		
		// High (66)
		params.velBias = 0.5f;
		params.adaptiveDivs = 7;
		params.adaptiveRings = 3;
		params.adaptiveDepth = 3;
		
		crowd->setObstacleAvoidanceParams(3, &params);
	}
}

void CrowdToolState::reset()
{
	m_targetRef=0;
}

static unsigned int ColorPedCount(int f, int max)
{
    unsigned int col;
	float a=128;
	float b=0;
	float g=0;
	float r=0;

	int lmax= max+10;
	int lmin=lmax/4;

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

	col = duRGBA(r,g,b,a);
	return col;
}

static int maxMap (std::map<int,int> frequencyCount)
{
	int currentMax = 0;

	for(auto it = frequencyCount.cbegin(); it != frequencyCount.cend(); ++it ) 
	{
		if (it ->second > currentMax) 
		{
			currentMax = it->second;
		}
	}
	return currentMax;
}

void CrowdToolState::handleRender()
{
	DebugDrawGL dd;
	const float rad = m_sample->getAgentRadius();
	
	dtNavMesh* nav = m_sample->getNavMesh();
	dtCrowd* crowd = m_sample->getCrowd();
	if (!nav || !crowd)
		return;
	
	if (m_toolParams.m_showNodes && crowd->getPathQueue())
	{
		const dtNavMeshQuery* navquery = crowd->getPathQueue()->getNavQuery();
		if (navquery)
			duDebugDrawNavMeshNodes(&dd, *navquery);
	}

	dd.depthMask(false);
	if (m_toolParams.m_showFrequentedRoads)
	{
			BcGraph& mcg = m_bucMesh->getBucGraph();
			bucket_property_t& bc_prop=m_bucMesh->getBc_prop();
			line_property_t& ln_prop=m_bucMesh->getLine_prop();
			b_edge_iterator_t ei, ei_end;
			m_roadCounter = crowd->getRoadCounter();	
			int max=maxMap(m_roadCounter);
			dd.begin(DU_DRAW_LINES, 2.5);
			for (boost::tie(ei, ei_end) = edges(mcg); ei != ei_end; ++ei)
			{
				struct bcLine* L=(ln_prop)[*ei];
				int areaId=L->bcLineRef;
			
				int f = m_roadCounter[areaId];
				unsigned int col = ColorPedCount(f,max);
				float v0[3],v1[3];
				dtVcopy(v0,L->bucketL->centre);
				dtVcopy(v1,L->bucketR->centre);
				dd.vertex(v0, col);
				dd.vertex(v1, col);
			}
			dd.end();
	}

	// Occupancy grid.
	if (m_toolParams.m_showGrid)
	{
		float gridy = -FLT_MAX;
		for (int i = 0; i < crowd->getAgentCount(); ++i)
		{
			const dtCrowdAgent* ag = crowd->getAgent(i);
			if (!ag->active) continue;
			const float* pos = ag->corridor.getPos();
			gridy = dtMax(gridy, pos[1]);
		}
		gridy += 1.0f;
		
		dd.begin(DU_DRAW_QUADS);
		const dtProximityGrid* grid = crowd->getGrid();
		const int* bounds = grid->getBounds();
		const float cs = grid->getCellSize();
		for (int y = bounds[1]; y <= bounds[3]; ++y)
		{
			for (int x = bounds[0]; x <= bounds[2]; ++x)
			{
				const int count = grid->getItemCountAt(x,y); 
				if (!count) continue;
				unsigned int col = duRGBA(128,0,0,dtMin(count*40,255));
				dd.vertex(x*cs, gridy, y*cs, col);
				dd.vertex(x*cs, gridy, y*cs+cs, col);
				dd.vertex(x*cs+cs, gridy, y*cs+cs, col);
				dd.vertex(x*cs+cs, gridy, y*cs, col);
			}
		}
		dd.end();
	}
	
	// Corners & co
	for (int i = 0; i < crowd->getAgentCount(); i++)
	{
		if (m_toolParams.m_showDetailAll == false && i != m_agentDebug.idx)
			continue;
		const dtCrowdAgent* ag =crowd->getAgent(i);
		if (!ag->active)
			continue;
			
		const float radius = ag->params.radius;
		const float* pos = ag->npos;
		
		if (m_toolParams.m_showCorners)
		{
			if (ag->ncorners)
			{
				dd.begin(DU_DRAW_LINES, 2.0f);
				for (int j = 0; j < ag->ncorners; ++j)
				{
					const float* va = j == 0 ? pos : &ag->cornerVerts[(j-1)*3];
					const float* vb = &ag->cornerVerts[j*3];
					dd.vertex(va[0],va[1]+radius,va[2], duRGBA(128,0,0,192));
					dd.vertex(vb[0],vb[1]+radius,vb[2], duRGBA(128,0,0,192));
				}
				if (ag->ncorners && ag->cornerFlags[ag->ncorners-1] & DT_STRAIGHTPATH_OFFMESH_CONNECTION)
				{
					const float* v = &ag->cornerVerts[(ag->ncorners-1)*3];
					dd.vertex(v[0],v[1],v[2], duRGBA(192,0,0,192));
					dd.vertex(v[0],v[1]+radius*2,v[2], duRGBA(192,0,0,192));
				}
				
				dd.end();
				
				
				if (m_toolParams.m_anticipateTurns)
				{
					
				}
			}
		}
		
		if (m_toolParams.m_showNeis)
		{
			duDebugDrawCircle(&dd, pos[0],pos[1]+radius,pos[2], ag->params.collisionQueryRange,
							  duRGBA(0,192,128,128), 2.0f);
			
			dd.begin(DU_DRAW_LINES, 2.0f);
			for (int j = 0; j < ag->nneis; ++j)
			{
				// Get 'n'th active agent.
				// TODO: fix this properly.
				const dtCrowdAgent* nei = crowd->getAgent(ag->neis[j].idx);
				if (nei)
				{
					dd.vertex(pos[0],pos[1]+radius,pos[2], duRGBA(0,192,128,128));
					dd.vertex(nei->npos[0],nei->npos[1]+radius,nei->npos[2], duRGBA(0,192,128,128));
				}
			}
			dd.end();
		}
		
		if (m_toolParams.m_showOpt)
		{
			dd.begin(DU_DRAW_LINES, 2.0f);
			dd.vertex(m_agentDebug.optStart[0],m_agentDebug.optStart[1]+0.3f,m_agentDebug.optStart[2], duRGBA(0,128,0,192));
			dd.vertex(m_agentDebug.optEnd[0],m_agentDebug.optEnd[1]+0.3f,m_agentDebug.optEnd[2], duRGBA(0,128,0,192));
			dd.end();
		}
	}
	
	// Agent cylinders.
	for (int i = 0; i < crowd->getAgentCount(); ++i)
	{
		const dtCrowdAgent* ag = crowd->getAgent(i);
		if (!ag->active) continue;
		
		const float radius = ag->params.radius;
		const float* pos = ag->npos;
		
		unsigned int col = duRGBA(0,0,0,32);
		if (m_agentDebug.idx == i)
			col = duRGBA(255,0,0,128);
			
		duDebugDrawCircle(&dd, pos[0], pos[1], pos[2], radius, col, 2.0f);
	}
	
	for (int i = 0; i < crowd->getAgentCount(); ++i)
	{
		const dtCrowdAgent* ag = crowd->getAgent(i);
		if (!ag->active) continue;
		
		const float height = ag->params.height;
		const float radius = ag->params.radius;
		const float* pos = ag->npos;
		
		unsigned int col = duRGBA(220,220,220,128);
		if (ag->targetState == DT_CROWDAGENT_TARGET_REQUESTING || ag->targetState == DT_CROWDAGENT_TARGET_WAITING_FOR_QUEUE)
			col = duLerpCol(col, duRGBA(128,0,255,128), 32);
		else if (ag->targetState == DT_CROWDAGENT_TARGET_WAITING_FOR_PATH)
			col = duLerpCol(col, duRGBA(128,0,255,128), 128);
		else if (ag->targetState == DT_CROWDAGENT_TARGET_FAILED)
			col = duRGBA(255,32,16,128);
		else if (ag->targetState == DT_CROWDAGENT_TARGET_VELOCITY)
			col = duLerpCol(col, duRGBA(64,255,0,128), 128);
		
		duDebugDrawCylinder(&dd, pos[0]-radius, pos[1]+radius*0.1f, pos[2]-radius,
							pos[0]+radius, pos[1]+height, pos[2]+radius, col);
	}
	
	
	if (m_toolParams.m_showVO)
	{
		for (int i = 0; i < crowd->getAgentCount(); i++)
		{
			if (m_toolParams.m_showDetailAll == false && i != m_agentDebug.idx)
				continue;
			const dtCrowdAgent* ag =crowd->getAgent(i);
			if (!ag->active)
				continue;
		
			// Draw detail about agent sela
			const dtObstacleAvoidanceDebugData* vod = m_agentDebug.vod;
			
			const float dx = ag->npos[0];
			const float dy = ag->npos[1]+ag->params.height;
			const float dz = ag->npos[2];
			
			duDebugDrawCircle(&dd, dx,dy,dz, ag->params.maxSpeed, duRGBA(255,255,255,64), 2.0f);
			
			dd.begin(DU_DRAW_QUADS);
			for (int j = 0; j < vod->getSampleCount(); ++j)
			{
				const float* p = vod->getSampleVelocity(j);
				const float sr = vod->getSampleSize(j);
				const float pen = vod->getSamplePenalty(j);
				const float pen2 = vod->getSamplePreferredSidePenalty(j);
				unsigned int col = duLerpCol(duRGBA(255,255,255,220), duRGBA(128,96,0,220), (int)(pen*255));
				col = duLerpCol(col, duRGBA(128,0,0,220), (int)(pen2*128));
				dd.vertex(dx+p[0]-sr, dy, dz+p[2]-sr, col);
				dd.vertex(dx+p[0]-sr, dy, dz+p[2]+sr, col);
				dd.vertex(dx+p[0]+sr, dy, dz+p[2]+sr, col);
				dd.vertex(dx+p[0]+sr, dy, dz+p[2]-sr, col);
			}
			dd.end();
		}
	}
	
	// Velocity stuff.
	for (int i = 0; i < crowd->getAgentCount(); ++i)
	{
		const dtCrowdAgent* ag = crowd->getAgent(i);
		if (!ag->active) continue;
		
		const float radius = ag->params.radius;
		const float height = ag->params.height;
		const float* pos = ag->npos;
		const float* vel = ag->vel;
		const float* dvel = ag->dvel;
		
		unsigned int col = duRGBA(220,220,220,192);
		if (ag->targetState == DT_CROWDAGENT_TARGET_REQUESTING || ag->targetState == DT_CROWDAGENT_TARGET_WAITING_FOR_QUEUE)
			col = duLerpCol(col, duRGBA(128,0,255,192), 32);
		/*else if (ag->targetState == DT_CROWDAGENT_TARGET_WAITING_FOR_PATH)
			col = duLerpCol(col, duRGBA(128,0,255,192), 128);*/
		else if (ag->targetState == DT_CROWDAGENT_TARGET_FAILED)
			col = duRGBA(255,32,16,192);
		/*else if (ag->targetState == DT_CROWDAGENT_TARGET_VELOCITY)
			col = duLerpCol(col, duRGBA(64,255,0,192), 128);*/
		
		duDebugDrawCircle(&dd, pos[0], pos[1]+height, pos[2], radius, col, 2.0f);
	}
	
	dd.depthMask(true);
}

void CrowdToolState::handleRenderOverlay(double* proj, double* model, int* view)
{
	GLdouble x, y, z;
	
	// Draw start and end point labels
	if (m_targetRef && gluProject((GLdouble)m_targetPos[0], (GLdouble)m_targetPos[1], (GLdouble)m_targetPos[2],
								  model, proj, view, &x, &y, &z))
	{
		//imguiDrawText((int)x, (int)(y+25), IMGUI_ALIGN_CENTER, "TARGET", imguiRGBA(0,0,0,220));
	}
	
	char label[32];
	
	if (m_toolParams.m_showNodes)
	{
		dtCrowd* crowd = m_sample->getCrowd();
		if (crowd && crowd->getPathQueue())
		{
			const dtNavMeshQuery* navquery = crowd->getPathQueue()->getNavQuery();
			const dtNodePool* pool = navquery->getNodePool();
			if (pool)
			{
				const float off = 0.5f;
				for (int i = 0; i < pool->getHashSize(); ++i)
				{
					for (dtNodeIndex j = pool->getFirst(i); j != DT_NULL_IDX; j = pool->getNext(j))
					{
						const dtNode* node = pool->getNodeAtIdx(j+1);
						if (!node) continue;

						if (gluProject((GLdouble)node->pos[0],(GLdouble)node->pos[1]+off,(GLdouble)node->pos[2],
									   model, proj, view, &x, &y, &z))
						{
							const float heuristic = node->total;// - node->cost;
							snprintf(label, 32, "%.2f", heuristic);
							imguiDrawText((int)x, (int)y+15, IMGUI_ALIGN_CENTER, label, imguiRGBA(0,0,0,220));
						}
					}
				}
			}
		}
	}
	
	if (m_toolParams.m_showLabels)
	{
		dtCrowd* crowd = m_sample->getCrowd();
		bcBucketMesh* bucMesh = m_sample->getBucketMesh();
		BcGraph& bucg = bucMesh->getBucGraph();
		bucket_property_t& bc_prop=bucMesh->getBc_prop();
		line_property_t& ln_prop=bucMesh->getLine_prop();
		b_edge_iterator_t ei, ei_end;
		if (crowd)
		{
			for (boost::tie(ei, ei_end) = edges(bucg); ei != ei_end; ++ei)
			{
				struct bcLine* L=(ln_prop)[*ei];
				int areaId=L->bcLineRef;
				
				float pos[3];
				dtVcopy(pos,L->midpoint);
				if (gluProject((GLdouble)pos[0], (GLdouble)pos[1]+0.5, (GLdouble)pos[2],
							   model, proj, view, &x, &y, &z))
				{
					snprintf(label, 32, "%d", areaId);
					imguiDrawText((int)x, (int)y+15, IMGUI_ALIGN_CENTER, label, imguiRGBA(0,0,0,220));
				}
				
			}
		}
	}
	if (m_agentDebug.idx != -1)
	{
		dtCrowd* crowd = m_sample->getCrowd();
		if (crowd) 
		{
			for (int i = 0; i < crowd->getAgentCount(); i++)
			{
				if (m_toolParams.m_showDetailAll == false && i != m_agentDebug.idx)
					continue;
				const dtCrowdAgent* ag =crowd->getAgent(i);
				if (!ag->active)
					continue;
				const float radius = ag->params.radius;
				if (m_toolParams.m_showNeis)
				{
					for (int j = 0; j < ag->nneis; ++j)
					{
						const dtCrowdAgent* nei = crowd->getAgent(ag->neis[j].idx);
						if (!nei->active) continue;
						
						if (gluProject((GLdouble)nei->npos[0], (GLdouble)nei->npos[1]+radius, (GLdouble)nei->npos[2],
									   model, proj, view, &x, &y, &z))
						{
							snprintf(label, 32, "%.3f", ag->neis[j].dist);
							imguiDrawText((int)x, (int)y+15, IMGUI_ALIGN_CENTER, label, imguiRGBA(255,255,255,220));
						}
					}
				}
			}
		}
	}
	
	
	if (m_toolParams.m_showPerfGraph)
	{
		GraphParams gp;
		gp.setRect(300, 10, 500, 200, 8);
		gp.setValueRange(0.0f, 2.0f, 4, "ms");
		
		drawGraphBackground(&gp);
		drawGraph(&gp, &m_crowdTotalTime, 1, "Total", duRGBA(255,128,0,255));
		
		gp.setRect(300, 10, 500, 50, 8);
		gp.setValueRange(0.0f, 2000.0f, 1, "");
		drawGraph(&gp, &m_crowdSampleCount, 0, "Sample Count", duRGBA(96,96,96,128));
	}
	
}

void CrowdToolState::handleUpdate(const float dt)
{
	if (m_run)
		updateTick(dt);
}

static int rand_n(int n) 
{ 
	srand(time(NULL));
    int nbgen=rand()%n+1;    
    
    return nbgen;
}

static int rand_o(int n) 
{ 
	srand(time(NULL));
    int nbgen=rand()%n;    
    
    return nbgen;
}

static float rand_m(int f) 
{ 
    float t[2];
	t[0]=1.0; t[1]=1/f;
	int n=1;
	int partSize   = 1 + (n == RAND_MAX ? 0 : (RAND_MAX - n) / (n + 1)); 
    int maxUsefull = partSize * n + (partSize - 1); 
    int draw; 
  
    do { 
        draw = rand(); 
    } while (draw > maxUsefull); 
  
    return t[draw / partSize]; 
}

void CrowdToolState::addAgent(const float* p,dtCrowdAgentCategory c,unsigned int f)
{
	if (!m_sample) return;
	dtCrowd* crowd = m_sample->getCrowd();
	bcBucketMesh* bucMesh= m_sample->getBucketMesh();
	BcGraph& bg = bucMesh->getBucGraph();

	dtCrowdAgentParams ap;
	memset(&ap, 0, sizeof(ap));
	ap.category=c;

	if (c==DT_CROWDAGENT_NOT_FAMILIAR)
	{
		ap.familiarity=0;
	}
	else
	{ 
		unsigned int e=rand_o(1);   
		ap.familiarity=f;
		ap.exploration=e;
	}

	ap.radius = m_sample->getAgentRadius();
	ap.height = m_sample->getAgentHeight();
	ap.maxAcceleration = 8.0f;
	ap.maxSpeed = 3.5f;
	ap.collisionQueryRange = ap.radius * 12.0f;
	ap.pathOptimizationRange = ap.radius * 30.0f;
	ap.updateFlags = 0; 
	if (m_toolParams.m_anticipateTurns)
		ap.updateFlags |= DT_CROWD_ANTICIPATE_TURNS;
	if (m_toolParams.m_optimizeVis)
		ap.updateFlags |= DT_CROWD_OPTIMIZE_VIS;
	if (m_toolParams.m_optimizeTopo)
		ap.updateFlags |= DT_CROWD_OPTIMIZE_TOPO;
	if (m_toolParams.m_obstacleAvoidance)
		ap.updateFlags |= DT_CROWD_OBSTACLE_AVOIDANCE;
	if (m_toolParams.m_separation)
		ap.updateFlags |= DT_CROWD_SEPARATION;
	ap.obstacleAvoidanceType = (unsigned char)m_toolParams.m_obstacleAvoidanceType;
	ap.separationWeight = m_toolParams.m_separationWeight;
	
	int idx = crowd->addAgent(p, &ap);
	if (idx != -1)
	{
		//m_agentsEnv[idx]=1;

		if (m_targetRef)
			crowd->requestMoveTarget(idx, m_targetRef, m_targetPos,bg);
	}
}

void CrowdToolState::removeAgent(const int idx)
{
	if (!m_sample) return;
	dtCrowd* crowd = m_sample->getCrowd();

	crowd->removeAgent(idx);
	
	if (idx == m_agentDebug.idx)
		m_agentDebug.idx = -1;
}

void CrowdToolState::removeAgents()
{
	if (!m_sample) return;
	dtNavMesh* nav = m_sample->getNavMesh();
	dtCrowd* crowd = m_sample->getCrowd();
	bcBucketMesh* bucMesh = m_sample->getBucketMesh();
	
		crowd->init(MAX_AGENTS, m_sample->getAgentRadius(), nav,bucMesh);
			// Make polygons with 'disabled' flag invalid.
		crowd->getEditableFilter(0)->setExcludeFlags(SAMPLE_POLYFLAGS_DISABLED);
		
		// Setup local avoidance params to different qualities.
		dtObstacleAvoidanceParams params;
		// Use mostly default settings, copy from dtCrowd.
		memcpy(&params, crowd->getObstacleAvoidanceParams(0), sizeof(dtObstacleAvoidanceParams));
		
		// Low (11)
		params.velBias = 0.5f;
		params.adaptiveDivs = 5;
		params.adaptiveRings = 2;
		params.adaptiveDepth = 1;
		crowd->setObstacleAvoidanceParams(0, &params);
		
		// Medium (22)
		params.velBias = 0.5f;
		params.adaptiveDivs = 5; 
		params.adaptiveRings = 2;
		params.adaptiveDepth = 2;
		crowd->setObstacleAvoidanceParams(1, &params);
		
		// Good (45)
		params.velBias = 0.5f;
		params.adaptiveDivs = 7;
		params.adaptiveRings = 2;
		params.adaptiveDepth = 3;
		crowd->setObstacleAvoidanceParams(2, &params);
		
		// High (66)
		params.velBias = 0.5f;
		params.adaptiveDivs = 7;
		params.adaptiveRings = 3;
		params.adaptiveDepth = 3;
		
		crowd->setObstacleAvoidanceParams(3, &params);
		

}

void CrowdToolState::hilightAgent(const int idx)
{
	m_agentDebug.idx = idx;
}

static void calcVel(float* vel, const float* pos, const float* tgt, const float speed)
{
	dtVsub(vel, tgt, pos);
	vel[1] = 0.0;
	dtVnormalize(vel);
	dtVscale(vel, vel, speed);
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

void CrowdToolState::setMoveTarget(const float* p, bool adjust,int n)
{
	if (!m_sample) return;

	// Find nearest point on navmesh and set move request to that location.
	dtNavMeshQuery* navquery = m_sample->getNavMeshQuery();
	dtCrowd* crowd = m_sample->getCrowd();
	bcBucketMesh* bucMesh= m_sample->getBucketMesh();
	BcGraph& bg = bucMesh->getBucGraph();
	const dtQueryFilter* filter = crowd->getFilter(0);
	const float* ext = crowd->getQueryExtents();

	if (adjust)
	{
		float vel[3];
		// Request velocity
		if (m_agentDebug.idx != -1)
		{
			const dtCrowdAgent* ag = crowd->getAgent(m_agentDebug.idx);
			if (ag && ag->active)
			{
				calcVel(vel, ag->npos, p, ag->params.maxSpeed);
				crowd->requestMoveVelocity(m_agentDebug.idx, vel);
			}
		}
		else
		{
			for (int i = 0; i < crowd->getAgentCount(); ++i)
			{
				const dtCrowdAgent* ag = crowd->getAgent(i);
				if (!ag->active) continue;
				calcVel(vel, ag->npos, p, ag->params.maxSpeed);
				crowd->requestMoveVelocity(i, vel);
			}
		}
	}
	else
	{
		if (m_agentDebug.idx != -1)
		{
			navquery->findNearestPoly(p, ext, filter, &m_targetRef, m_targetPos);
			const dtCrowdAgent* ag = crowd->getAgent(m_agentDebug.idx);
			if (ag && ag->active)
				crowd->requestMoveTarget(m_agentDebug.idx, m_targetRef, m_targetPos,bg);
		}
		else
		{
			if (n==1)
			{
				navquery->findNearestPoly(p, ext, filter, &m_targetRef, m_targetPos);
				for (int i = 0; i < crowd->getAgentCount(); ++i)
				{
					const dtCrowdAgent* ag = crowd->getAgent(i);
					if (!ag->active) continue;
					crowd->requestMoveTarget(i, m_targetRef, m_targetPos,bg);
					const dtPolyRef* path = ag->corridor.getPath();
					dtPolyRef startRef=path[0];
					float spos[3];
					float epos[3];
					dtVcopy(spos, ag->npos);
					navquery->findPath(startRef, m_targetRef,spos, m_targetPos, filter, m_polys, &m_npolys, MAX_POLYS);
					m_nstraightPath = 0;
					if (m_npolys)
					{
						// In case of partial path, make sure the end point is clamped to the last polygon.
						dtVcopy(epos, m_targetPos);
						if (m_polys[m_npolys-1] != m_targetRef)
						navquery->closestPointOnPoly(m_polys[m_npolys-1], m_targetPos, epos, 0);
				
						navquery->findStraightPath(spos, epos, m_polys, m_npolys,m_straightPath, m_straightPathFlags,
									m_straightPathPolys, &m_nstraightPath, MAX_POLYS, 0);
						
						float d=calcDist(m_straightPath, m_nstraightPath);
						
					}

				}
			}
			else
			{
					int ncrowd=crowd->getAgentCount();
				
					if (n==ncrowd)
					{
						for (int i = 0; i < ncrowd; ++i)
						{
							const dtCrowdAgent* ag = crowd->getAgent(i);
							navquery->findNearestPoly(&p[i*3], ext, filter, &m_targetRef, m_targetPos);
							if (!ag->active) continue;
							crowd->requestMoveTarget(i, m_targetRef, m_targetPos,bg);
						}
					}
					else
					{
						for (int i = 0; i < ncrowd; ++i)
						{
							const dtCrowdAgent* ag = crowd->getAgent(i);
							int idx=i % n; 
							navquery->findNearestPoly(&p[idx*3], ext, filter, &m_targetRef, m_targetPos);
							if (!ag->active) continue;
							crowd->requestMoveTarget(i, m_targetRef, m_targetPos,bg);
						}
					}
				
			}
		}
	}
}

int CrowdToolState::hitTestAgents(const float* s, const float* p)
{
	if (!m_sample) return -1;
	dtCrowd* crowd = m_sample->getCrowd();
	
	int isel = -1;
	float tsel = FLT_MAX;

	for (int i = 0; i < crowd->getAgentCount(); ++i)
	{
		const dtCrowdAgent* ag = crowd->getAgent(i);
		if (!ag->active) continue;
		float bmin[3], bmax[3];
		getAgentBounds(ag, bmin, bmax);
		float tmin, tmax;
		if (isectSegAABB(s, p, bmin,bmax, tmin, tmax))
		{
			if (tmin > 0 && tmin < tsel)
			{
				isel = i;
				tsel = tmin;
			} 
		}
	}

	return isel;
}

void CrowdToolState::updateAgentParams()
{
	if (!m_sample) return;
	dtCrowd* crowd = m_sample->getCrowd();
	if (!crowd) return;
	
	unsigned char updateFlags = 0;
	unsigned char obstacleAvoidanceType = 0;
	
	if (m_toolParams.m_anticipateTurns)
		updateFlags |= DT_CROWD_ANTICIPATE_TURNS;
	if (m_toolParams.m_optimizeVis)
		updateFlags |= DT_CROWD_OPTIMIZE_VIS;
	if (m_toolParams.m_optimizeTopo)
		updateFlags |= DT_CROWD_OPTIMIZE_TOPO;
	if (m_toolParams.m_obstacleAvoidance)
		updateFlags |= DT_CROWD_OBSTACLE_AVOIDANCE;
	if (m_toolParams.m_obstacleAvoidance)
		updateFlags |= DT_CROWD_OBSTACLE_AVOIDANCE;
	if (m_toolParams.m_separation)
		updateFlags |= DT_CROWD_SEPARATION;
	
	obstacleAvoidanceType = (unsigned char)m_toolParams.m_obstacleAvoidanceType;
	
	dtCrowdAgentParams params;
	
	for (int i = 0; i < crowd->getAgentCount(); ++i)
	{
		const dtCrowdAgent* ag = crowd->getAgent(i);
		if (!ag->active) continue;
		memcpy(&params, &ag->params, sizeof(dtCrowdAgentParams));
		params.updateFlags = updateFlags;
		params.obstacleAvoidanceType = obstacleAvoidanceType;
		params.separationWeight = m_toolParams.m_separationWeight;
		crowd->updateAgentParameters(i, &params);
	}	
}

void CrowdToolState::updateTick(const float dt)
{
	if (!m_sample) return;
	dtNavMesh* nav = m_sample->getNavMesh();
	dtCrowd* crowd = m_sample->getCrowd();
    bcBucketMesh* bucMesh= m_sample->getBucketMesh();
	BcGraph& bg = bucMesh->getBucGraph();

	if (!nav || !crowd) return;
	
	TimeVal startTime = getPerfTime();
	
	crowd->update(dt, &m_agentDebug,bg);
	
	TimeVal endTime = getPerfTime();
	
	m_agentDebug.vod->normalizeSamples();
	
	m_crowdSampleCount.addSample((float)crowd->getVelocitySampleCount());
	m_crowdTotalTime.addSample(getPerfDeltaTimeUsec(startTime, endTime) / 1000.0f);
}


CrowdTool::CrowdTool() :
	m_sample(0),
	m_state(0),
	m_randDestination(0),
	m_mode(TOOLMODE_CREATE)
{
	resetCommonSettings();
}

CrowdTool::~CrowdTool()
{
}

void CrowdTool::resetCommonSettings()
{
	if (m_randDestination)
		free(m_randDestination);
	m_randDestination=0;

	m_Population_Size=0.0;
	
	for(int i=0;i<4;i++)
		m_PopulationTab[i]=0.0;
	
	limH=100.0;limM=100.0;limL=100.0;limU=100.0;
	m_Destinations_Numb=1.0f;

}

void CrowdTool::init(Sample* sample)
{
	if (m_sample != sample)
	{
		m_sample = sample;
	}
	
	if (!sample)
		return;
		
	m_state = (CrowdToolState*)sample->getToolState(type());
	if (!m_state)
	{
		m_state = new CrowdToolState();
		sample->setToolState(type(), m_state);
	}

	m_state->init(sample);
}

void CrowdTool::reset()
{
	if (!m_sample) return;
	if (!m_state) return;
	m_state->setRunning(true);
	m_state->removeAgents();
	m_state->reset();
	resetCommonSettings();
	
}

int CrowdTool::randPoints(dtNavMeshQuery* navQuery,float* randPointsTab,float nrandPoints)
{
	int n = 0;
	int nt=(int)nrandPoints;
	while (n<nt)
	{
		float pt[3];
		dtPolyRef ref;
		dtStatus status = navQuery->findRandomPoint(&m_filter, frand, &ref, pt);
		if (dtStatusSucceed(status))
		{
			const dtMeshTile* t = 0;
			const dtPoly* p = 0;
			dtNavMesh* m_nav=m_sample->getNavMesh();
			m_nav->getTileAndPolyByRefUnsafe(ref, &t, &p);
			if (p->type == DT_POLYTYPE_CITY_GROUND)
			{
				dtVcopy(&randPointsTab[n*3], pt);
				n++;
			}
		}
	}
	return n;
}

int CrowdTool::randPoints2D(dtNavMeshQuery* navQuery,float* randPointsTab,float nrandPoints,float* h)
{
	int n = 0;
	int nt=(int)nrandPoints;
	while (n<nt)
	{
		float pt[3];
		float ptt[2];
		dtPolyRef ref;
		dtStatus status = navQuery->findRandomPoint(&m_filter, frand, &ref, pt);
		*h=pt[1];
		if (dtStatusSucceed(status))
		{
			const dtMeshTile* t = 0;
			const dtPoly* p = 0;
			dtNavMesh* m_nav=m_sample->getNavMesh();
			m_nav->getTileAndPolyByRefUnsafe(ref, &t, &p);
			if (p->type == DT_POLYTYPE_CITY_GROUND)
			{
				ptt[0]=pt[0];
				ptt[1]=pt[2];
				dtVcopy2D(&randPointsTab[n*2], ptt);
				n++;
			}
		}
	}
	return n;
}

float* CrowdTool::randPoint(dtNavMeshQuery* navQuery)
{
		float pt[3];
		dtPolyRef ref;
		dtStatus status = navQuery->findRandomPoint(&m_filter, frand, &ref, pt);
		if (dtStatusSucceed(status))
		{
			return pt;
		}
		return NULL;
}

void CrowdTool::hundleAddAgents(float n,float* populationTab)
{
	if (!m_sample) return;
	if (!m_state) return;
	InputGeom* geom = m_sample->getInputGeom();
	if (!geom) return;
	dtCrowd* crowd = m_sample->getCrowd();
	if (!crowd) return;
	dtNavMeshQuery* m_navQuery = m_sample->getNavMeshQuery();

	int notFamPopSize=int(n*populationTab[0]/100);
	if (notFamPopSize)
	{
		float h=0;float t[3];
		float* randPointN = (float*) malloc(sizeof(float)*notFamPopSize*2);
		int not=randPoints2D(m_navQuery,randPointN,notFamPopSize,&h);
		for (int i = 0; i < not; i++)
			{
				t[0]=randPointN[i*2];
				t[1]=h;
				t[2]=randPointN[i*2+1];
				m_state->addAgent(t,DT_CROWDAGENT_NOT_FAMILIAR,0);
			}
		free(randPointN);
	}

	int hFamPopSize=int(n*populationTab[1]/100);
	if (hFamPopSize)
	{
		float h=0;float t[3];
		float* randPointH = (float*) malloc(sizeof(float)*hFamPopSize*2);
		int hi=randPoints2D(m_navQuery,randPointH,hFamPopSize,&h);
		for (int i = 0; i < hi; i++)
		{
			t[0]=randPointH[i*2];
			t[1]=h;
			t[2]=randPointH[i*2+1];
			m_state->addAgent(t,DT_CROWDAGENT_FAMILIAR,1);
		}
		free(randPointH);
	}

	int mFamPopSize=int(n*populationTab[2]/100);
	if (mFamPopSize)
	{
		float h=0;float t[3];
		float* randPointM = (float*) malloc(sizeof(float)*mFamPopSize*2);
		int me=randPoints2D(m_navQuery,randPointM,mFamPopSize,&h);
		for (int i = 0; i < me; i++)
		{
			t[0]=randPointM[i*2];
			t[1]=h;
			t[2]=randPointM[i*2+1];
			m_state->addAgent(t,DT_CROWDAGENT_FAMILIAR,2);
		}
		free(randPointM);
	}

	int lFamPopSize=int(n*populationTab[3]/100);
	if (lFamPopSize)
	{																															
		float h=0;float t[3];
		float* randPointL = (float*) malloc(sizeof(float)*lFamPopSize*2);
		int lo=randPoints2D(m_navQuery,randPointL,lFamPopSize,&h);
		for (int i = 0; i < lo; i++)
		{
			t[0]=randPointL[i*2];
			t[1]=h;
			t[2]=randPointL[i*2+1];
			m_state->addAgent(t,DT_CROWDAGENT_FAMILIAR,3);
		}
		free(randPointL);
	}
	
}

void CrowdTool::hundleAddDestinations(bool grouped,int n)
{
	
	if (!m_sample) return;
	if (grouped)
	{
		m_randDestination = (float*) malloc(sizeof(float)*n*3);
		if (m_randDestination)
		{
			dtNavMeshQuery* navquery = m_sample->getNavMeshQuery();
			int nrandPoints=randPoints(navquery,m_randDestination,n);
			if (nrandPoints!=0)
				m_state->setMoveTarget(m_randDestination,false,nrandPoints);
			else
				return;
		}
		else
		{
			return;
		}
	}
	else
	{
		m_randDestination = (float*) malloc(sizeof(float)*MAX_RAND_POINTS*3);
		if (m_randDestination)
		{
			dtNavMeshQuery* navquery = m_sample->getNavMeshQuery();
			int nrandPoints=randPoints(navquery,m_randDestination,MAX_RAND_POINTS);
			if (nrandPoints!=0)
				m_state->setMoveTarget(m_randDestination,false,nrandPoints);
			else
				return;
		}
		else
		{
			return;
		}
	}
}

void CrowdTool::hundleSave()
{
	dtCrowd* crowd = m_sample->getCrowd();
	
	FILE* fpr = fopen("results.txt", "w");
	if (!fpr) return;

	bcBucketMesh* bucMesh = m_sample->getBucketMesh();
	BcGraph& bucg = bucMesh->getBucGraph();
	bucket_property_t& bc_prop=bucMesh->getBc_prop();
	line_property_t& ln_prop=bucMesh->getLine_prop();
	b_edge_iterator_t ei, ei_end;
	std::map<int,int> roadCounter = crowd->getRoadCounter();

	fprintf(fpr,"areaId;roadCounter;betweeness\n");
	for (boost::tie(ei, ei_end) = edges(bucg); ei != ei_end; ++ei)
	{
		struct bcLine* L=(ln_prop)[*ei];
		int areaId=L->bcLineRef;
		fprintf(fpr,"%d;%d;%f\n",areaId,roadCounter[areaId],L->betweeness);	
	}	
				
	fclose(fpr);
}

void CrowdTool::handleMenu()
{
	if (!m_state)
		return;
	CrowdToolParams* params = m_state->getToolParams();

	imguiLabel("Population Params");

	if (imguiSlider("Size", &m_Population_Size, 0.0f, 3000.0f, 1.0f,m_Population_Size_set))
		{
			m_mode = TOOLMODE_CREATE;
			m_Population_Size_set=true;
		}

	imguiSeparator();
	imguiLabel("Familiarity Level");
	

	limH=100-(m_PopulationTab[2]+m_PopulationTab[3]+m_PopulationTab[0]);
	if (imguiSlider("High", &m_PopulationTab[1], 0.0f, limH, 1.0f,m_Agents_HighFam_set))
		{
			m_mode = TOOLMODE_CREATE;
			m_Agents_HighFam_set=true;
		}
	
	limM=100-(m_PopulationTab[1]+m_PopulationTab[3]+m_PopulationTab[0]);
	if (imguiSlider("Medium", &m_PopulationTab[2], 0.0f, limM, 1.0f,m_Agents_MediumFam_set))
		{
			m_mode = TOOLMODE_CREATE;
			m_Agents_MediumFam_set=true;
		}
	float d;
	limL=100-(m_PopulationTab[2]+m_PopulationTab[1]+m_PopulationTab[0]);
	if (imguiSlider("Low", &m_PopulationTab[3], 0.0f, limL, 1.0f,m_Agents_LowFam_set))
		{
			m_mode = TOOLMODE_CREATE;
			m_Agents_LowFam_set=true;
		}
	
	limU=100-(m_PopulationTab[2]+m_PopulationTab[3]+m_PopulationTab[1]);
	if (imguiSlider("Unfamiliar", &m_PopulationTab[0], 0.0f, limU, 1.0f,m_Agents_Not_Fam_set))
		{
			m_mode = TOOLMODE_CREATE;
			m_Agents_Not_Fam_set=true;
		}

	if (imguiButton("Add agents"))
	{
		if ((m_Agents_HighFam_set)||(m_Agents_MediumFam_set)||(m_Agents_LowFam_set)||(m_Agents_Not_Fam_set))
			hundleAddAgents(int(m_Population_Size),&m_PopulationTab[0]);
		
	}
	

	imguiSeparatorLine();
	imguiLabel("Destinations");

	if (imguiCheck("Random Destinations", m_mode == TOOLMODE_MOVE_ARBITRARY))
	{
		m_mode = TOOLMODE_MOVE_ARBITRARY;
		hundleAddDestinations(false,int(m_Population_Size));
	}

	if (imguiSlider("Random Common Destinations", &m_Destinations_Numb, 0.0f, 50.0f, 1.0f,m_Destinations_Numb_set))
	{
		m_mode = TOOLMODE_MOVE_DESTINATION;
		m_Destinations_Numb_set=true;
	}
	imguiIndent();
	if (imguiButton("Generate"))
	{
		if (m_Destinations_Numb_set)
			hundleAddDestinations(m_Destinations_Numb_set,int(m_Destinations_Numb));
		
	}
	imguiUnindent();
	imguiSeparatorLine();

	if (imguiCollapse("Results", 0, params->m_expandSelectedDebugDraw))
		params->m_expandSelectedDebugDraw = !params->m_expandSelectedDebugDraw;
		
	if (params->m_expandSelectedDebugDraw)
	{
		imguiIndent();
		if (imguiCheck("Show frequented roads pattern", params->m_showFrequentedRoads))
			params->m_showFrequentedRoads = !params->m_showFrequentedRoads;
		imguiUnindent();
	}

	imguiSeparatorLine();

	if (imguiButton("Save Results"))
	{
			hundleSave();
	}

	if (imguiButton("Reset"))
	{
			reset();
	}
}

void CrowdTool::handleClick(const float* s, const float* p, bool shift)
{
	if (!m_sample) return;
	if (!m_state) return;
	InputGeom* geom = m_sample->getInputGeom();
	if (!geom) return;
	dtCrowd* crowd = m_sample->getCrowd();
	if (!crowd) return;

	if (m_mode == TOOLMODE_CREATE)
	{
		if (shift)
		{
			// Delete
			int ahit = m_state->hitTestAgents(s,p);
			if (ahit != -1)
				m_state->removeAgent(ahit);
		}
		else
		{
			// Add
			m_state->addAgent(p,DT_CROWDAGENT_FAMILIAR,1);
		}
	}
	else if (m_mode == TOOLMODE_MOVE_TARGET)
	{
		m_state->setMoveTarget(p, shift,1);
	}
	else if (m_mode == TOOLMODE_SELECT)
	{
		// Highlight
		int ahit = m_state->hitTestAgents(s,p);
		m_state->hilightAgent(ahit);
	}
	else if (m_mode == TOOLMODE_TOGGLE_POLYS)
	{
		dtNavMesh* nav = m_sample->getNavMesh();
		dtNavMeshQuery* navquery = m_sample->getNavMeshQuery();
		if (nav && navquery)
		{
			dtQueryFilter filter;
			const float* ext = crowd->getQueryExtents();
			float tgt[3];
			dtPolyRef ref;
			navquery->findNearestPoly(p, ext, &filter, &ref, tgt);
			if (ref)
			{
				unsigned short flags = 0;
				if (dtStatusSucceed(nav->getPolyFlags(ref, &flags)))
				{
					flags ^= SAMPLE_POLYFLAGS_DISABLED;
					nav->setPolyFlags(ref, flags);
				}
			}
		}
	}

}

void CrowdTool::handleStep()
{
	if (!m_state) return;
	
	const float dt = 1.0f/20.0f;
	m_state->updateTick(dt);

	m_state->setRunning(false);
}

void CrowdTool::handleToggle()
{
	if (!m_state) return;
	m_state->setRunning(!m_state->isRunning());
}

void CrowdTool::handleUpdate(const float dt)
{
	rcIgnoreUnused(dt);
}

void CrowdTool::handleRender()
{
}

void CrowdTool::handleRenderOverlay(double* proj, double* model, int* view)
{
	rcIgnoreUnused(model);
	rcIgnoreUnused(proj);

	// Tool help
	const int h = view[3];
	int ty = h-40;
	
	if (m_mode == TOOLMODE_MOVE_TARGET)
	{
		imguiDrawText(280, ty, IMGUI_ALIGN_LEFT, "LMB: set move target.  Shift+LMB: adjust set velocity.", imguiRGBA(255,255,255,192));	
		ty -= 20;
		imguiDrawText(280, ty, IMGUI_ALIGN_LEFT, "Setting velocity will move the agents without pathfinder.", imguiRGBA(255,255,255,192));	
	}
	else if (m_mode == TOOLMODE_SELECT)
	{
		imguiDrawText(280, ty, IMGUI_ALIGN_LEFT, "LMB: select agent.", imguiRGBA(255,255,255,192));	
	}
	ty -= 20;
	imguiDrawText(280, ty, IMGUI_ALIGN_LEFT, "SPACE: Run/Pause simulation.  1: Step simulation.", imguiRGBA(255,255,255,192));	
	ty -= 20;

	if (m_state && m_state->isRunning())
		imguiDrawText(280, ty, IMGUI_ALIGN_LEFT, "- RUNNING -", imguiRGBA(255,32,16,255));	
	else 
		imguiDrawText(280, ty, IMGUI_ALIGN_LEFT, "- PAUSED -", imguiRGBA(255,255,255,128));	
}
