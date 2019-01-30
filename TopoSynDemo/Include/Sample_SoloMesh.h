//

//

#ifndef RECASTSAMPLESOLOMESH_H
#define RECASTSAMPLESOLOMESH_H

#include "Sample.h"
#include "DetourNavMesh.h"
#include "Recast.h"

#include "VoroDiag.h"
#include "BucketBuilder.h"

class Sample_SoloMesh : public Sample
{
protected:
	bool m_keepInterResults;
	float m_totalBuildTimeMs;
	float m_totalBucketBuildTimeS;

	unsigned char* m_triareas;
	rcHeightfield* m_solid;
	rcCompactHeightfield* m_chf;
	rcContourSet* m_cset;
	rcPolyMesh* m_pmesh;
	rcConfig m_cfg;	
	rcPolyMeshDetail* m_dmesh;
	int m_dim;
	float m_ground;
	
	VoronoiDiagram m_vd;
	polygon_t m_poly;
	Graph m_graph;
	MaGraph m_ma_graph;
	MaGraph m_ma_graph_f;
	
	std::vector<PointType> m_point_data;
	std::vector<SegmentType> m_segment_data;
	bcBucketSet* m_bset;
	double m_maxBet[3];

	compMap_t m_compMap;
	coor_property_t m_coorMap;
	BcGraph m_BucketGraph;
	BcGraph m_VisibilityGraph;
	ConGraph m_ConnectivityGraph;
	bucket_property_t m_bc_prop;
	bucket_property_t m_bcv_prop;
	line_property_t m_line_prop;
	contourVector_t m_borderCont;
	std::vector<contourVector_t> m_holesConts;
	int m_numborderComponents;

	enum DrawMode
	{
		DRAWMODE_NAVMESH,
		DRAWMODE_NAVMESH_TRANS,
		DRAWMODE_NAVMESH_INVIS,
		DRAWMODE_MESH,
		DRAWMODE_BUCKETS,
		DRAWMODE_CONNECTIVITY,
		DRAWMODE_VISIBILITY,
		DRAWMODE_COLORMAP,
		DRAWMODE_NAVMESH_CONTOURS,
		MAX_DRAWMODE
	};
	
	DrawMode m_drawMode;
	
	void cleanup();
		
public:
	Sample_SoloMesh();
	virtual ~Sample_SoloMesh();
	
	virtual void handleSettings();
	virtual void handleTools();
	virtual void handleDebugMode();
	
	virtual void handleRender();
	virtual void handleRenderOverlay(double* proj, double* model, int* view);
	virtual void handleMeshChanged(class InputGeom* geom);
	virtual bool handleBuild();
	virtual bool handleBuildBuckets();

	std::vector<ma_vertex_t> m_ma_vertex_vector;
};


#endif // RECASTSAMPLESOLOMESHSIMPLE_H
