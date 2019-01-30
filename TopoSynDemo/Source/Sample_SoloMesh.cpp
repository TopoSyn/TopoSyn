//

//

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "SDL.h"
#include "SDL_opengl.h"
#include "imgui.h"
#include "InputGeom.h"
#include "Sample.h"
#include "Sample_SoloMesh.h"
#include "Recast.h"
#include "RecastDebugDraw.h"
#include "RecastDump.h"
#include "DetourNavMesh.h"
#include "DetourNavMeshBuilder.h"
#include "DetourDebugDraw.h"
#include "NavMeshTesterTool.h"
#include "NavMeshPruneTool.h"
#include "OffMeshConnectionTool.h"
#include "ConvexVolumeTool.h"
#include "CrowdTool.h"

#include "BucketDebugDraw.h"        
#include "BucketMeshTool.h"
#include "BucketBuilder.h"

#ifdef WIN32
#	define snprintf _snprintf
#endif


Sample_SoloMesh::Sample_SoloMesh() :
	m_keepInterResults(false),
	m_totalBuildTimeMs(0),
	m_totalBucketBuildTimeS(0),
	m_triareas(0),
	m_solid(0),
	m_chf(0),
	m_cset(0),
	m_pmesh(0),
	m_dmesh(0),
	m_dim(0),
	m_bset(0),
	m_drawMode(DRAWMODE_NAVMESH)
{
	setTool(new BucketMeshTool);
}
		
Sample_SoloMesh::~Sample_SoloMesh()
{
	cleanup();
}
	
void Sample_SoloMesh::cleanup()
{
	delete [] m_triareas;
	m_triareas = 0;
	rcFreeHeightField(m_solid);
	m_solid = 0;
	rcFreeCompactHeightfield(m_chf);
	m_chf = 0;
	rcFreeContourSet(m_cset);
	m_cset = 0;
	rcFreePolyMesh(m_pmesh);
	m_pmesh = 0;
	rcFreePolyMeshDetail(m_dmesh);
	m_dmesh = 0;
	dtFreeNavMesh(m_navMesh);
	m_navMesh = 0;
	m_graph.clear();
	m_vd.clear();
	m_ma_graph.clear();
	bcFreeBucketSet(m_bset);
	bcFreeBucketMesh(m_bucMesh);
	m_bset=0;
	m_bucMesh=0;
	m_BucketGraph.clear();
	m_ConnectivityGraph.clear();
	m_VisibilityGraph.clear();

	
}
			
void Sample_SoloMesh::handleSettings()
{
	Sample::handleCommonSettings();

	imguiSeparator();
	
	char msg[64];
	snprintf(msg, 64, "Build Time: %.1fms", m_totalBuildTimeMs);
	imguiLabel(msg);
	
	imguiSeparator();

	char msgb[64];
	snprintf(msgb, 64, "Bucket Build Time: %.1fs", m_totalBucketBuildTimeS);
	imguiLabel(msgb);
	
	imguiSeparator();
}

void Sample_SoloMesh::handleTools()
{
	int type = !m_tool ? TOOL_NONE : m_tool->type();
	
	if (imguiCheck("Create Crowds", type == TOOL_CROWD))
	{
		setTool(new CrowdTool);
	}

	if (imguiCheck("Test Bucket Mesh", type == TOOL_BUCKETMESH))
	{
		setTool(new BucketMeshTool);
	}
	
	imguiSeparatorLine();

	imguiIndent();

	if (m_tool)
		m_tool->handleMenu();

	imguiUnindent();

}

void Sample_SoloMesh::handleDebugMode()
{
	// Check which modes are valid.
	bool valid[MAX_DRAWMODE];
	for (int i = 0; i < MAX_DRAWMODE; ++i)
		valid[i] = false;

	if (m_geom)
	{
		valid[DRAWMODE_MESH] = true;
		valid[DRAWMODE_NAVMESH] = m_navMesh != 0;
		valid[DRAWMODE_BUCKETS] = m_BucketGraph.m_edges.size() != 0;  
		valid[DRAWMODE_COLORMAP] = m_BucketGraph.m_edges.size() != 0;  
	}
	
	int unavail = 0;
	for (int i = 0; i < MAX_DRAWMODE; ++i)
		if (!valid[i]) unavail++;

	if (unavail == MAX_DRAWMODE)
		return;

	imguiLabel("Draw");
	if (imguiCheck("Input Mesh", m_drawMode == DRAWMODE_MESH, valid[DRAWMODE_MESH]))
		m_drawMode = DRAWMODE_MESH;
	if (imguiCheck("Navmesh", m_drawMode == DRAWMODE_NAVMESH, valid[DRAWMODE_NAVMESH]))
		m_drawMode = DRAWMODE_NAVMESH;
	if (imguiCheck("Buckets mesh", m_drawMode == DRAWMODE_BUCKETS, valid[DRAWMODE_BUCKETS]))    
		m_drawMode = DRAWMODE_BUCKETS;
	if (imguiCheck("Buckets colored map", m_drawMode == DRAWMODE_COLORMAP, valid[DRAWMODE_COLORMAP]))    
		m_drawMode = DRAWMODE_COLORMAP;
}

void Sample_SoloMesh::handleRender()
{
	if (!m_geom || !m_geom->getMesh())
		return;
	
	DebugDrawGL dd;
	
	glEnable(GL_FOG);
	glDepthMask(GL_TRUE);

	const float texScale = 1.0f / (m_cellSize * 10.0f);
	
	// Draw mesh
	if (m_drawMode != DRAWMODE_COLORMAP)
	{
		duDebugDrawTriMeshSlope(&dd, m_geom->getMesh()->getVerts(), m_geom->getMesh()->getVertCount(),
								m_geom->getMesh()->getTris(), m_geom->getMesh()->getNormals(), m_geom->getMesh()->getTriCount(),
								m_agentMaxSlope, texScale);
		m_geom->drawOffMeshConnections(&dd);
	}

	glDisable(GL_FOG);
	glDepthMask(GL_FALSE);

	// Draw bounds
	const float* bmin = m_geom->getMeshBoundsMin();
	const float* bmax = m_geom->getMeshBoundsMax();
	duDebugDrawBoxWire(&dd, bmin[0],bmin[1],bmin[2], bmax[0],bmax[1],bmax[2], duRGBA(255,255,255,128), 1.0f);
	dd.begin(DU_DRAW_POINTS, 5.0f);
	dd.vertex(bmin[0],bmin[1],bmin[2],duRGBA(255,255,255,128));
	dd.end();
	
	if (m_navMesh && m_navQuery)
	{
		if (m_drawMode == DRAWMODE_NAVMESH)
			duDebugDrawNavMeshWithClosedList(&dd, *m_navMesh, *m_navQuery, m_navMeshDrawFlags);

		if ((m_drawMode == DRAWMODE_BUCKETS) && m_BucketGraph.m_edges.size())
			bcDebugDrawBucketMesh(&dd, *m_navMesh, *m_navQuery, m_navMeshDrawFlags);
		
		duDebugDrawNavMeshPolysWithFlags(&dd, *m_navMesh, SAMPLE_POLYFLAGS_DISABLED, duRGBA(0,0,0,128));
	}
		
	glDepthMask(GL_TRUE);

	if ((m_drawMode == DRAWMODE_COLORMAP)&& (m_BucketGraph.m_edges.size()))
	{
			glDepthMask(GL_FALSE);
			bcDebugDrawBucketGraphColored(&dd,m_bucMesh,m_BucketGraph,bmin[1],m_maxBet[2]);
			glDepthMask(GL_TRUE);
	}
	
	m_geom->drawConvexVolumes(&dd);

	if (m_tool)
		m_tool->handleRender();
	renderToolStates();

	glDepthMask(GL_TRUE);
}

void Sample_SoloMesh::handleRenderOverlay(double* proj, double* model, int* view)
{
	if (m_tool)
		m_tool->handleRenderOverlay(proj, model, view);
	renderOverlayToolStates(proj, model, view);
}

void Sample_SoloMesh::handleMeshChanged(class InputGeom* geom)
{
	Sample::handleMeshChanged(geom);

	dtFreeNavMesh(m_navMesh);
	bcFreeBucketSet(m_bset);
	bcFreeBucketMesh(m_bucMesh);
	m_navMesh = 0;
	m_bset=0;
	m_bucMesh=0;
	m_vd.clear();
	m_graph.clear();
	m_ma_graph.clear();
	m_BucketGraph.clear();
	m_ConnectivityGraph.clear();
	m_VisibilityGraph.clear();

	if (m_tool)
	{
		m_tool->reset();
		m_tool->init(this);
	}
	resetToolStates();
	initToolStates(this);
}

bool Sample_SoloMesh::handleBuild()
{
	if (!m_geom || !m_geom->getMesh())
	{
		m_ctx->log(RC_LOG_ERROR, "buildNavigation: Input mesh is not specified.");
		return false;
	}
	
	cleanup();
	
	const float* bmin = m_geom->getMeshBoundsMin();
	const float* bmax = m_geom->getMeshBoundsMax();
	const float* verts = m_geom->getMesh()->getVerts();
	const int nverts = m_geom->getMesh()->getVertCount();
	const int* tris = m_geom->getMesh()->getTris();
	const int ntris = m_geom->getMesh()->getTriCount();
	
	//
	// Step 1. Initialize build config.
	//
	
	// Init build configuration from GUI
	memset(&m_cfg, 0, sizeof(m_cfg));
	m_cfg.cs = m_cellSize;
	m_cfg.ch = m_cellHeight;
	m_cfg.walkableSlopeAngle = m_agentMaxSlope;
	m_cfg.walkableHeight = (int)ceilf(m_agentHeight / m_cfg.ch);
	m_cfg.walkableClimb = (int)floorf(m_agentMaxClimb / m_cfg.ch);
	m_cfg.walkableRadius = (int)ceilf(m_agentRadius / m_cfg.cs);
	m_cfg.maxEdgeLen = (int)(m_edgeMaxLen / m_cellSize);
	m_cfg.maxSimplificationError = m_edgeMaxError;
	m_cfg.minRegionArea = (int)rcSqr(m_regionMinSize);		// Note: area = size*size
	m_cfg.mergeRegionArea = (int)rcSqr(m_regionMergeSize);	// Note: area = size*size
	m_cfg.maxVertsPerPoly = (int)m_vertsPerPoly;
	m_cfg.detailSampleDist = m_detailSampleDist < 0.9f ? 0 : m_cellSize * m_detailSampleDist;
	m_cfg.detailSampleMaxError = m_cellHeight * m_detailSampleMaxError;
	
	// Set the area where the navigation will be build.
	// Here the bounds of the input mesh are used, but the
	// area could be specified by an user defined box, etc.
	rcVcopy(m_cfg.bmin, bmin);
	rcVcopy(m_cfg.bmax, bmax);
	rcCalcGridSize(m_cfg.bmin, m_cfg.bmax, m_cfg.cs, &m_cfg.width, &m_cfg.height);

	// Reset build times gathering.
	m_ctx->resetTimers();

	// Start the build process.	
	m_ctx->startTimer(RC_TIMER_TOTAL);
	
	m_ctx->log(RC_LOG_PROGRESS, "Building navigation:");
	m_ctx->log(RC_LOG_PROGRESS, " - %d x %d cells", m_cfg.width, m_cfg.height);
	m_ctx->log(RC_LOG_PROGRESS, " - %.1fK verts, %.1fK tris", nverts/1000.0f, ntris/1000.0f);
	
	//
	// Step 2. Rasterize input polygon soup.
	//
	
	// Allocate voxel heightfield where we rasterize our input data to.
	m_solid = rcAllocHeightfield();
	if (!m_solid)
	{
		m_ctx->log(RC_LOG_ERROR, "buildNavigation: Out of memory 'solid'.");
		return false;
	}
	if (!rcCreateHeightfield(m_ctx, *m_solid, m_cfg.width, m_cfg.height, m_cfg.bmin, m_cfg.bmax, m_cfg.cs, m_cfg.ch))
	{
		m_ctx->log(RC_LOG_ERROR, "buildNavigation: Could not create solid heightfield.");
		return false;
	}
	
	// Allocate array that can hold triangle area types.
	// If you have multiple meshes you need to process, allocate
	// and array which can hold the max number of triangles you need to process.
	m_triareas = new unsigned char[ntris];
	if (!m_triareas)
	{
		m_ctx->log(RC_LOG_ERROR, "buildNavigation: Out of memory 'm_triareas' (%d).", ntris);
		return false;
	}
	
	// Find triangles which are walkable based on their slope and rasterize them.
	// If your input data is multiple meshes, you can transform them here, calculate
	// the are type for each of the meshes and rasterize them.
	memset(m_triareas, 0, ntris*sizeof(unsigned char));
	rcMarkWalkableTriangles(m_ctx, m_cfg.walkableSlopeAngle, verts, nverts, tris, ntris, m_triareas);
	rcRasterizeTriangles(m_ctx, verts, nverts, tris, m_triareas, ntris, *m_solid, m_cfg.walkableClimb);

	if (!m_keepInterResults)
	{
		delete [] m_triareas;
		m_triareas = 0;
	}
	
	//
	// Step 3. Filter walkables surfaces.
	//
	
	// Once all geoemtry is rasterized, we do initial pass of filtering to
	// remove unwanted overhangs caused by the conservative rasterization
	// as well as filter spans where the character cannot possibly stand.
	rcFilterLowHangingWalkableObstacles(m_ctx, m_cfg.walkableClimb, *m_solid);
	rcFilterLedgeSpans(m_ctx, m_cfg.walkableHeight, m_cfg.walkableClimb, *m_solid);
	rcFilterWalkableLowHeightSpans(m_ctx, m_cfg.walkableHeight, *m_solid);


	//
	// Step 4. Partition walkable surface to simple regions.
	//

	// Compact the heightfield so that it is faster to handle from now on.
	// This will result more cache coherent data as well as the neighbours
	// between walkable cells will be calculated.
	m_chf = rcAllocCompactHeightfield();
	if (!m_chf)
	{
		m_ctx->log(RC_LOG_ERROR, "buildNavigation: Out of memory 'chf'.");
		return false;
	}
	if (!rcBuildCompactHeightfield(m_ctx, m_cfg.walkableHeight, m_cfg.walkableClimb, *m_solid, *m_chf))
	{
		m_ctx->log(RC_LOG_ERROR, "buildNavigation: Could not build compact data.");
		return false;
	}
	
	if (!m_keepInterResults)
	{
		rcFreeHeightField(m_solid);
		m_solid = 0;
	}
		
	// Erode the walkable area by agent radius.
	if (!rcErodeWalkableArea(m_ctx, m_cfg.walkableRadius, *m_chf))
	{
		m_ctx->log(RC_LOG_ERROR, "buildNavigation: Could not erode.");
		return false;
	}

	// (Optional) Mark areas.
	const ConvexVolume* vols = m_geom->getConvexVolumes();
	for (int i  = 0; i < m_geom->getConvexVolumeCount(); ++i)
		rcMarkConvexPolyArea(m_ctx, vols[i].verts, vols[i].nverts, vols[i].hmin, vols[i].hmax, (unsigned char)vols[i].area, *m_chf);

	
	// Partition the heightfield so that we can use simple algorithm later to triangulate the walkable areas.
	// There are 3 martitioning methods, each with some pros and cons:
	// 1) Watershed partitioning
	//   - the classic Recast partitioning
	//   - creates the nicest tessellation
	//   - usually slowest
	//   - partitions the heightfield into nice regions without holes or overlaps
	//   - the are some corner cases where this method creates produces holes and overlaps
	//      - holes may appear when a small obstacles is close to large open area (triangulation can handle this)
	//      - overlaps may occur if you have narrow spiral corridors (i.e stairs), this make triangulation to fail
	//   * generally the best choice if you precompute the nacmesh, use this if you have large open areas
	// 2) Monotone partioning
	//   - fastest
	//   - partitions the heightfield into regions without holes and overlaps (guaranteed)
	//   - creates long thin polygons, which sometimes causes paths with detours
	//   * use this if you want fast navmesh generation
	// 3) Layer partitoining
	//   - quite fast
	//   - partitions the heighfield into non-overlapping regions
	//   - relies on the triangulation code to cope with holes (thus slower than monotone partitioning)
	//   - produces better triangles than monotone partitioning
	//   - does not have the corner cases of watershed partitioning
	//   - can be slow and create a bit ugly tessellation (still better than monotone)
	//     if you have large open areas with small obstacles (not a problem if you use tiles)
	//   * good choice to use for tiled navmesh with medium and small sized tiles
	
	if (m_partitionType == SAMPLE_PARTITION_WATERSHED)
	{
		// Prepare for region partitioning, by calculating distance field along the walkable surface.
		if (!rcBuildDistanceField(m_ctx, *m_chf))
		{
			m_ctx->log(RC_LOG_ERROR, "buildNavigation: Could not build distance field.");
			return false;
		}
		
		// Partition the walkable surface into simple regions without holes.
		if (!rcBuildRegions(m_ctx, *m_chf, 0, m_cfg.minRegionArea, m_cfg.mergeRegionArea))
		{
			m_ctx->log(RC_LOG_ERROR, "buildNavigation: Could not build watershed regions.");
			return false;
		}
	}
	else if (m_partitionType == SAMPLE_PARTITION_MONOTONE)
	{
		// Partition the walkable surface into simple regions without holes.
		// Monotone partitioning does not need distancefield.
		if (!rcBuildRegionsMonotone(m_ctx, *m_chf, 0, m_cfg.minRegionArea, m_cfg.mergeRegionArea))
		{
			m_ctx->log(RC_LOG_ERROR, "buildNavigation: Could not build monotone regions.");
			return false;
		}
	}
	else // SAMPLE_PARTITION_LAYERS
	{
		// Partition the walkable surface into simple regions without holes.
		if (!rcBuildLayerRegions(m_ctx, *m_chf, 0, m_cfg.minRegionArea))
		{
			m_ctx->log(RC_LOG_ERROR, "buildNavigation: Could not build layer regions.");
			return false;
		}
	}
	
	//
	// Step 5. Trace and simplify region contours.
	//
	
	// Create contours.
	m_cset = rcAllocContourSet();
	if (!m_cset)
	{
		m_ctx->log(RC_LOG_ERROR, "buildNavigation: Out of memory 'cset'.");
		return false;
	}
	if (!rcBuildContours(m_ctx, *m_chf, m_cfg.maxSimplificationError, m_cfg.maxEdgeLen, *m_cset))
	{
		m_ctx->log(RC_LOG_ERROR, "buildNavigation: Could not create contours.");
		return false;
	}
	
	//
	// Step 6. Build polygons mesh from contours.
	//
	
	// Build polygon navmesh from the contours.
	m_pmesh = rcAllocPolyMesh();
	if (!m_pmesh)
	{
		m_ctx->log(RC_LOG_ERROR, "buildNavigation: Out of memory 'pmesh'.");
		return false;
	}
	if (!rcBuildPolyMesh(m_ctx, *m_cset, m_cfg.maxVertsPerPoly, *m_pmesh))
	{
		m_ctx->log(RC_LOG_ERROR, "buildNavigation: Could not triangulate contours.");
		return false;
	}
	
	//
	// Step 7. Create detail mesh which allows to access approximate height on each polygon.
	//
	
	m_dmesh = rcAllocPolyMeshDetail();
	if (!m_dmesh)
	{
		m_ctx->log(RC_LOG_ERROR, "buildNavigation: Out of memory 'pmdtl'.");
		return false;
	}

	if (!rcBuildPolyMeshDetail(m_ctx, *m_pmesh, *m_chf, m_cfg.detailSampleDist, m_cfg.detailSampleMaxError, *m_dmesh))
	{
		m_ctx->log(RC_LOG_ERROR, "buildNavigation: Could not build detail mesh.");
		return false;
	}

	if (!m_keepInterResults)
	{
		rcFreeCompactHeightfield(m_chf);
		m_chf = 0;
		rcFreeContourSet(m_cset);
		m_cset = 0;
	}

	// At this point the navigation mesh data is ready, you can access it from m_pmesh.
	// See duDebugDrawPolyMesh or dtCreateNavMeshData as examples how to access the data.
	
	//
	// (Optional) Step 8. Create Detour data from Recast poly mesh.
	//
	
	// The GUI may allow more max points per polygon than Detour can handle.
	// Only build the detour navmesh if we do not exceed the limit.
	if (m_cfg.maxVertsPerPoly <= DT_VERTS_PER_POLYGON)
	{
		unsigned char* navData = 0;
		int navDataSize = 0;

		// Update poly flags from areas.
		for (int i = 0; i < m_pmesh->npolys; ++i)
		{
			if (m_pmesh->areas[i] == RC_WALKABLE_AREA)
				m_pmesh->areas[i] = SAMPLE_POLYAREA_GROUND;
				
			if (m_pmesh->areas[i] == SAMPLE_POLYAREA_GROUND ||
				m_pmesh->areas[i] == SAMPLE_POLYAREA_ROAD)
			{
				m_pmesh->flags[i] = SAMPLE_POLYFLAGS_WALK;
			}
			else if (m_pmesh->areas[i] == SAMPLE_POLYAREA_CROSSROAD)
			{
				m_pmesh->flags[i] = SAMPLE_POLYFLAGS_WALK | SAMPLE_POLYFLAGS_DECIDE;
			}
		}


		dtNavMeshCreateParams params;
		memset(&params, 0, sizeof(params));
		params.verts = m_pmesh->verts;
		params.vertCount = m_pmesh->nverts;
		params.polys = m_pmesh->polys;
		params.polyAreas = m_pmesh->areas;
		params.polyFlags = m_pmesh->flags;
		params.polyCount = m_pmesh->npolys;
		params.nvp = m_pmesh->nvp;
		params.detailMeshes = m_dmesh->meshes;
		params.detailVerts = m_dmesh->verts;
		params.detailVertsCount = m_dmesh->nverts;
		params.detailTris = m_dmesh->tris;
		params.detailTriCount = m_dmesh->ntris;
		params.offMeshConVerts = m_geom->getOffMeshConnectionVerts();
		params.offMeshConRad = m_geom->getOffMeshConnectionRads();
		params.offMeshConDir = m_geom->getOffMeshConnectionDirs();
		params.offMeshConAreas = m_geom->getOffMeshConnectionAreas();
		params.offMeshConFlags = m_geom->getOffMeshConnectionFlags();
		params.offMeshConUserID = m_geom->getOffMeshConnectionId();
		params.offMeshConCount = m_geom->getOffMeshConnectionCount();
		params.walkableHeight = m_agentHeight;
		params.walkableRadius = m_agentRadius;
		params.walkableClimb = m_agentMaxClimb;
		rcVcopy(params.bmin, m_pmesh->bmin);
		rcVcopy(params.bmax, m_pmesh->bmax);
		params.cs = m_cfg.cs;
		params.ch = m_cfg.ch;
		params.buildBvTree = true;
		
		if (!dtCreateNavMeshData(&params, &navData, &navDataSize))
		{
			m_ctx->log(RC_LOG_ERROR, "Could not build Detour navmesh.");
			return false;
		}
		
		m_navMesh = dtAllocNavMesh();
		if (!m_navMesh)
		{
			dtFree(navData);
			m_ctx->log(RC_LOG_ERROR, "Could not create Detour navmesh");
			return false;
		}
		
		dtStatus status;
		
		status = m_navMesh->init(navData, navDataSize, DT_TILE_FREE_DATA);
		if (dtStatusFailed(status))
		{
			dtFree(navData);
			m_ctx->log(RC_LOG_ERROR, "Could not init Detour navmesh");
			return false;
		}
		
		status = m_navQuery->init(m_navMesh, 2048);
		if (dtStatusFailed(status))
		{
			m_ctx->log(RC_LOG_ERROR, "Could not init Detour navmesh query");
			return false;
		}
	}
	
	m_ctx->stopTimer(RC_TIMER_TOTAL);

	// Show performance stats.
	duLogBuildTimes(*m_ctx, m_ctx->getAccumulatedTime(RC_TIMER_TOTAL));
	m_ctx->log(RC_LOG_PROGRESS, ">> Polymesh: %d vertices  %d polygons", m_pmesh->nverts, m_pmesh->npolys);
	
	m_totalBuildTimeMs = m_ctx->getAccumulatedTime(RC_TIMER_TOTAL)/1000.0f;

	return true;
}

bool Sample_SoloMesh::handleBuildBuckets()
{
	if (!m_navMesh)
	{
		m_ctx->log(RC_LOG_ERROR, "BuildBuckets: Navigation mesh is not specified.");
		return false;
	}
	// Reset build times gathering.
	//m_ctx->resetTimers();
	// Start the build process.	
	m_ctx->startTimer(RC_TIMER_TOTAL);
	// Start the build process.
	m_ctx->startTimer(RC_TIMER_BUILD_BUCKETS_TOTAL);
	m_ctx->log(RC_LOG_PROGRESS, "Building Bucket graph:");

	//
	////
	//// Step 1. Extract contours *************************************************************************************
	////
	//

	const float* bmin = m_geom->getMeshBoundsMin();
	const float* bmax = m_geom->getMeshBoundsMax();
	m_ground = m_navQuery->findCityGroundPolys(m_navMesh);
	
	long int num=0;
	for (int i = 0; i < m_navMesh->getMaxTiles(); ++i)
	{
		const dtMeshTile* tile = m_navMesh->getTile(i);
	    num=num+bcBuildGraphContours(tile,*m_navMesh,m_graph,m_coorMap);
	}
	if (num==0)
	{
		m_ctx->log(RC_LOG_ERROR, "BuildBucket: failed to build Graph Contours.");
		return false;
	}
	
	bcExtractContours(m_graph,m_compMap);
	long int borderIdx=-1;
	int vectHolesContsSize=0;

	borderIdx = bcBuildBorderContour(m_graph,m_coorMap,m_compMap,bmin,bmax,m_agentRadius,m_borderCont);
	if (borderIdx==-1)
	{
		m_ctx->log(RC_LOG_ERROR, "BuildBucket: failed to build Border Contour.");
		return false;
	}

	long int holesIdx = bcBuilHolesContours(m_graph, m_coorMap, m_compMap,borderIdx, m_holesConts);
	if (m_holesConts.size()==0)
	{
		m_ctx->log(RC_LOG_ERROR, "BuildBucket: failed to build Holes Contours.");
		return false;
	}

	bcBuildWallsPoly(m_borderCont,m_holesConts,m_poly);
	int a=boost::geometry::area(m_poly);

	if (a==0)
	{
		m_ctx->log(RC_LOG_ERROR, "BuildBucket:: failed to build environment polygons.");
		return false;
	}

	m_bucMesh=bcAllocBucketMesh(m_BucketGraph,m_VisibilityGraph);
	if (!m_bucMesh)
	{
		m_ctx->log(RC_LOG_ERROR, "BuildBucket: Out of memory 'BucketMesh'.");
		return false;
	}

	m_bset = m_bucMesh->bcAllocBucketSet();
	if (!m_bset)
	{
		m_ctx->log(RC_LOG_ERROR, "BuildBucket: Out of memory 'BucketSet'.");
		return false;
	}

	bool suc=false;
	for (int i = 0; i < m_navMesh->getMaxTiles(); ++i)
	{
		const dtMeshTile* tile = m_navMesh->getTile(i);
		suc=m_bucMesh->bcBuildBucket(m_navQuery,tile,m_navMesh,*m_bset);
	}

	if (suc==false)
	{
		m_ctx->log(RC_LOG_ERROR, "BuildBucket: Could not build road buckets.");
		return false;
	}
	
	dtStatus status=m_bucMesh->bcSetBucketType(m_navMesh,*m_bset);
	if (status!=DT_SUCCESS)
	{
		m_ctx->log(RC_LOG_ERROR, "BuildBucket: Could not build crossroad buckets.");
		return false;
	}

	//
	////
	//// Step 4. Building Bucket Graph **********************************************************
	////
	//
	
	m_bucMesh->bcBuildBucketGraph(m_navMesh,m_navQuery,*m_bset,m_BucketGraph);
		
	if (boost::num_vertices(m_BucketGraph)>0)	
	{
		m_ctx->log(RC_LOG_ERROR, "Bucket graph %d edges %d vertices ",boost::num_edges(m_BucketGraph),boost::num_vertices(m_BucketGraph));
	}
	else 
	{
		m_ctx->log(RC_LOG_ERROR, "BuildBucket: Could not build Bucket graph ");
		return false;
	}

	//
	////
	//// Step 5. Building Connectivity Graph **********************************************************
	////
	//
	
	int n=m_bucMesh->bcBuildConnectivityGraph(m_BucketGraph,m_ConnectivityGraph);
	m_dim=boost::num_vertices(m_ConnectivityGraph);
	
	if (m_dim>0)
	{
		m_ctx->log(RC_LOG_ERROR, "Connectivity graph %d edges %d vertices",boost::num_edges(m_ConnectivityGraph), m_dim);
	}
	else 
	{
		m_ctx->log(RC_LOG_ERROR, "BuildBucket: Could not build Connectivity graph ");
		return false;
	}

	//
	////
	//// Step 6. Distance matrix*******************************************************************************
	////
	//
	
	typedef std::vector< std::vector<int> > matrix;
	matrix m_DistMat(m_dim, std::vector<int>(m_dim));

	for (int i=0;i<m_dim;++i)
		for (int j=0;j<m_dim;++j)
			m_DistMat[i][j]=0;

	int k=m_bucMesh->buildDistMat(m_ConnectivityGraph,m_DistMat);
	
	if (k>0)
	{
		//m_ctx->log(RC_LOG_PROGRESS, "Distance matrix dimension %d",k);
	}
	else 
	{
		m_ctx->log(RC_LOG_ERROR, "BuildBucket: Could not calculate Distance matrix ");
		return false;
	}

	//
	////
	//// Step 7. Space Syntax**************************************************************************************
	////
	//
	
	int nbr=m_bucMesh->bcCalcSSParam(m_ConnectivityGraph,m_BucketGraph,m_DistMat,m_dim,m_poly);//,m_b,m_c,m_i);
	m_bucMesh->bcCalcIntersectionBucParam(m_BucketGraph);
	
	////
	//////
	////// Step 8. Building Visibility Graph **********************************************************
	//////
	////
	//

	int nm=m_bucMesh->bcBuildVisibilityGraph(m_VisibilityGraph,m_BucketGraph); 
	int dv=boost::num_vertices(m_VisibilityGraph);
	if (dv!=0)
	{
		m_ctx->log(RC_LOG_ERROR, "Visibility graph %d edges %d vertices",boost::num_edges(m_VisibilityGraph),nm);
	}
	else 
	{
		m_ctx->log(RC_LOG_ERROR, "BuildBucket: Could not build Visibility graph ");
		return false;
	}


	////
	//////
	////// Step 9. Space Syntax parameters of the Visibility Graph **********************************************************
	//////
	////
	//

	
	m_bucMesh->bcCalcSSParammVG(m_navMesh,m_VisibilityGraph);
	
	//
	////
	//////
	////// Step 10. Building Roads *****************************************************************************
	//////
	////

	for (int i = 0; i < m_navMesh->getMaxTiles(); ++i)
	{
		const dtMeshTile* tile = m_navMesh->getTile(i);
		status=m_bucMesh->bcBuildRoads(m_navMesh,tile,m_BucketGraph);
	}

	if (status!=DT_SUCCESS)
	{
		m_ctx->log(RC_LOG_ERROR, "BuildBucket: Could not build roads buckets.");
		return false;
	}
	
	//// Free Distance Matrix 
	
	m_DistMat.clear();

	//
	////
	//// Show performance stats **************************************************************************************
	////
	//

	m_ctx->stopTimer(RC_TIMER_BUILD_BUCKETS_TOTAL);
	m_ctx->stopTimer(RC_TIMER_TOTAL);
	duLogBuildTimes(*m_ctx, m_ctx->getAccumulatedTime(RC_TIMER_BUILD_BUCKETS_TOTAL));
	m_ctx->log(RC_LOG_PROGRESS, ">> Bucket graph: %d vertices  %d edges", boost::num_edges(m_BucketGraph),boost::num_vertices(m_BucketGraph));
	m_totalBucketBuildTimeS = m_ctx->getAccumulatedTime(RC_TIMER_BUILD_BUCKETS_TOTAL)/1000000.0f;

	// Init tools ****************************************************************************************************

	bool classesStatus=m_bucMesh->bcBuildRoadClasses(m_BucketGraph,m_maxBet,3);
	if (classesStatus==false)
		return false;

	if (m_tool)
		m_tool->init(this);
	initToolStates(this);
	return true;

}
