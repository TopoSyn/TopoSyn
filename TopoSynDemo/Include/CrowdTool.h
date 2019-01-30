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

#ifndef CROWDTOOL_H
#define CROWDTOOL_H

#include "Sample.h"
#include "DetourNavMesh.h"
#include "DetourObstacleAvoidance.h"
#include "ValueHistory.h"
#include "DetourCrowd.h"

// Tool to create crowds.

struct CrowdToolParams
{
	bool m_expandSelectedDebugDraw;
	bool m_showCorners;
	bool m_showCollisionSegments;
	bool m_showFrequentedRoads;
	bool m_showPath;
	bool m_showBetweenness;
	bool m_showGlobalIntegration;
	bool m_showVO;
	bool m_showOpt;
	bool m_showNeis;
	
	bool m_expandDebugDraw;
	bool m_showLabels;
	bool m_showGrid;
	bool m_showNodes;
	bool m_showPerfGraph;
	bool m_showDetailAll;
	
	bool m_expandOptions;
	bool m_anticipateTurns;
	bool m_optimizeVis;
	bool m_optimizeTopo;
	bool m_obstacleAvoidance;
	float m_obstacleAvoidanceType;
	bool m_separation;
	float m_separationWeight;
};

class CrowdToolState : public SampleToolState
{
	Sample* m_sample;
	dtNavMesh* m_nav;
	dtCrowd* m_crowd;
	bcBucketMesh* m_bucMesh;
	
	float m_targetPos[3];
	dtPolyRef m_targetRef;

	dtCrowdAgentDebugInfo m_agentDebug;
	dtObstacleAvoidanceDebugData* m_vod;
	
	static const int MAX_AGENTS = 3000;

	ValueHistory m_crowdTotalTime;
	ValueHistory m_crowdSampleCount;

	CrowdToolParams m_toolParams;

	std::map<int,int> m_roadCounter;

	bool m_run;

	static const int MAX_POLYS =256;
	dtPolyRef m_polys[MAX_POLYS];
	float m_straightPath[MAX_POLYS*3];
	unsigned char m_straightPathFlags[MAX_POLYS];
	dtPolyRef m_straightPathPolys[MAX_POLYS];
	int m_nstraightPath;
	int m_npolys;

public:
	CrowdToolState();
	virtual ~CrowdToolState();
	
	virtual void init(class Sample* sample);
	virtual void reset();
	virtual void handleRender();
	virtual void handleRenderOverlay(double* proj, double* model, int* view);
	virtual void handleUpdate(const float dt);

	inline bool isRunning() const { return m_run; }
	inline void setRunning(const bool s) { m_run = s; }
	
	void addAgent(const float* pos,dtCrowdAgentCategory c,unsigned int f);
	void removeAgent(const int idx);
	void removeAgents();
	void hilightAgent(const int idx);
	void updateAgentParams();
	int hitTestAgents(const float* s, const float* p);
	void setMoveTarget(const float* p, bool adjust,int n);
	void updateTick(const float dt);

	inline CrowdToolParams* getToolParams() { return &m_toolParams; }
};


class CrowdTool : public SampleTool
{
	Sample* m_sample;
	CrowdToolState* m_state;
	float m_PopulationTab[4];
	float m_Population_Size;
	float m_Agents_HighFam,m_Agents_MediumFam,m_Agents_LowFam;
	float m_Agents_Not_Fam;
	bool m_Agents_HighFam_set,m_Agents_MediumFam_set,m_Agents_LowFam_set;
	bool m_Agents_Not_Fam_set,m_Agents_Fam_set;
	bool m_Population_Size_set;
	float limH,limM,limL,limU;
	float m_Destinations_Numb;
	bool m_Destinations_Numb_set;
	dtQueryFilter m_filter;
	static const int MAX_RAND_POINTS = 3000;
	float* m_randDestination;
	
	enum ToolMode
	{
		TOOLMODE_CREATE,
		TOOLMODE_MOVE_TARGET,
		TOOLMODE_MOVE_ARBITRARY,
		TOOLMODE_MOVE_DESTINATION,
		TOOLMODE_SELECT,
		TOOLMODE_TOGGLE_POLYS,
	};
	ToolMode m_mode;
	
	void updateAgentParams();
	void updateTick(const float dt);
	int randPoints(dtNavMeshQuery* navQuery,float* randPointsTab,float nrandPoints);
	int randPoints2D(dtNavMeshQuery* navQuery,float* randPointsTab,float nrandPoints,float* h);
	float* randPoint(dtNavMeshQuery* navQuery);
	
public:
	CrowdTool();
	virtual ~CrowdTool();
	void resetCommonSettings();
	virtual int type() { return TOOL_CROWD; }
	virtual void init(Sample* sample);
	virtual void reset();
	virtual void handleMenu();
	virtual void handleClick(const float* s, const float* p, bool shift);
	virtual void handleToggle();
	virtual void handleStep();
	virtual void handleUpdate(const float dt);
	virtual void handleRender();
	virtual void handleRenderOverlay(double* proj, double* model, int* view);
	void hundleAddAgents(float n,float* populationTab);
	void hundleAddDestinations(bool grouped,int n);
	void hundleSave();
};

#endif // CROWDTOOL_H
