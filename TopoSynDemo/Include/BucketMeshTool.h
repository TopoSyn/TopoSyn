//

//

#ifndef BUCKETMESHTOOL_H
#define BUCKETMESHTOOL_H

#include "Sample.h"
#include "DetourNavMesh.h"
#include "DetourNavMeshQuery.h"
#include "VoroDiag.h"
#include "BucketBuilder.h"

struct BucketToolParams
{
	bool m_showConnectivity;
	bool m_showIntegration;
	bool m_showPath;
	bool m_expandDebugDraw;
	bool m_showLabels;
	bool m_showNodes;
	bool m_showGraph;
	bool m_showDetailAll;
};

class BucketMeshTool : public SampleTool
{
	Sample* m_sample;
	
	dtNavMesh* m_navMesh;
	bcBucketMesh* m_bucMesh;
	dtNavMeshQuery* m_navQuery;
	BcGraph* m_bucg;
	dtQueryFilter m_filter;
	dtStatus m_pathFindStatus;


	enum ToolMode
	{
		TOOLMODE_ISOVIST_POLY,
		TOOLMODE_ISOVIST_RIDGE,
		TOOLMODE_RAYCAST,
		TOOLMODE_PATHFIND,
	};
	
	ToolMode m_toolMode;

	int m_PathFindStrategy;

	BucketToolParams m_toolParams;
	
	static const int MAX_POLYS = 256;
	static const int MAX_SMOOTH = 2048;
	
	
	dtPolyRef m_queryPointRef;
	dtPolyRef m_dirRef;
	dtPolyRef m_startRef;
	dtPolyRef m_endRef;
	dtPolyRef m_polys[MAX_POLYS];
	dtPolyRef m_parent[MAX_POLYS];
	int m_npolys;
	float m_straightPath[MAX_POLYS*3];
	unsigned char m_straightPathFlags[MAX_POLYS];
	dtPolyRef m_straightPathPolys[MAX_POLYS];
	int m_nstraightPath;
	float m_polyPickExt[3];
	float m_smoothPath[MAX_SMOOTH*3];
	int m_nsmoothPath;
	float m_queryPoly[4*3];
	
	float m_qpos[3];
	float m_dpos[3];
	bool m_qposSet;
	bool m_dposSet;
	float m_spos[3];
	float m_epos[3];
	bool m_sposSet;
	bool m_eposSet;
	bool m_Agent_FamSet;
	unsigned int m_Agent_Fam;


	int m_nverts;
	float* m_verts;

	float m_hitPosMax[3];
	float m_hitPosMin[3];
	float m_hitNormal[3];
	bool m_hitResult;
	float m_distanceToWall;

	float m_ridgeV0[3];
	float m_ridgeV1[3];
	
	MaGraph m_vg;
	std::vector<polygon_t> m_isovist_poly;
	int m_pathIterNum;
	const dtPolyRef* m_pathIterPolys; 
	int m_pathIterPolyCount;
	float m_prevIterPos[3], m_iterPos[3], m_steerPos[3], m_targetPos[3];
	
	static const int MAX_STEER_POINTS = 10;
	float m_steerPoints[MAX_STEER_POINTS*3];
	int m_steerPointCount;
	dtPolyRef* m_path;

	bucket_property_t* m_bcv_prop;
	line_property_t* m_line_prop;
	
	
public:
	BucketMeshTool();
	~BucketMeshTool();

	virtual int type() { return TOOL_BUCKETMESH; }
	virtual void init(Sample* sample);
	virtual void reset();
	virtual void handleMenu();
	virtual void handleClick(const float* s, const float* p, bool shift);
	virtual void handleToggle();
	virtual void handleStep();
	virtual void handleUpdate(const float dt);
	virtual void handleRender();
	virtual void handleRenderOverlay(double* proj, double* model, int* view);

	void drawQueryPoint(const float* pos, float r, float h, float c,const unsigned int col);
	void recalc();
	
};

#endif // BUCKETMESHTOOL_H