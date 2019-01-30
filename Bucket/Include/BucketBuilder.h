
//
#ifndef BUCKETBUILDER_H
#define BUCKETBUILDER_H

#include "DetourNavMesh.h"
#include "DetourNavMeshQuery.h"
#include "VoroDiag.h"
#include "BucketAlloc.h"

struct bcNeiLine
{
	unsigned short int neiLineRef;
	double farDistIso;
	float areaIso;
};

struct bcBucket
{
	unsigned short int bcBucRef;
	float centre[3];	
	int npolys;
	dtPolyRef* bucketPoly;
	int nportals;
	dtPolyRef* portals;
	dtPolyRef* portalsNeigh;
	bcNeiLine* neiLine;
	float* leftVert;
	float* rightVert;
	int connectivity;
	float globalIntegration;
	float localIntegration;
	int nverts;            
	float* verts;			
	float area;           
};

struct bcBucketSet
{
	bcBucket* buckets;	   
	int nbuckets;			
};

struct bcLine
{
	unsigned short int bcLineRef;
	unsigned int interval;
	bcBucket* bucketL;
	bcBucket* bucketR;
	dtPolyRef polyL;
	dtPolyRef polyR;
	int connectivity;
	float globalIntegration;
	float localIntegration;
	double betweeness;
	bool visibility;
	int npolys;
	dtPolyRef* linePoly;
	float area;  
	float midpoint[3];
};

struct bcLineStep
{
	unsigned short int bcLineRef;
	dtPolyRef interBuc;
};

typedef std::map<b_vertex_t,struct bcBucket*> bucket_property_t;
typedef std::map<b_edge_t,struct bcLine*> line_property_t;

typedef std::vector<float> polyVerts;

typedef std::vector< std::vector<int> > matrix;

void bcFreeBucketSet(bcBucketSet* bset);


class bcBucketMesh
{
public:

	bcBucketMesh(BcGraph& bcg,BcGraph& vbcg);
	~bcBucketMesh();

	void setBucGraph(BcGraph& bcg);
	void setVisiBucGraph(BcGraph& bcg);
	
	BcGraph& getBucGraph();
	bucket_property_t& bcBucketMesh::getBc_prop ();
	bucket_property_t& bcBucketMesh::getVisiBc_prop ();
	line_property_t& bcBucketMesh::getLine_prop ();

	void bcView(dtNavMesh* mesh,MaGraph& vg,const float* portalApex,dtPolyRef poly,dtPolyRef parentpoly,
		        float* A,float* B,float* Right,float* Left);

	int bcCalcIsovist(dtNavMeshQuery* navQuery,dtNavMesh* mesh,dtPolyRef queryPointRef,const float* queryPointPos,MaGraph& vg);

	void bcFindHit(dtPolyRef curRef,int seg,float& sa,float smax);

	int bcRaycast1(dtNavMesh* mesh,dtPolyRef centerRef, const float* startPos, const float* endPos, 
		           const dtQueryFilter* filter,double& sa);

	int bcRaycast2(dtNavMesh* mesh,dtPolyRef centerRef, const float* startPos, const float* endPos, 
		           const dtQueryFilter* filter,double& si);

	void bcHitPos(const float* spos,float* epos,float* hitPosMax, float* hitPosMin,double smin, double smax);

	void bcBuildIsovistPoly(MaGraph vg,std::vector<polygon_t>& polysVect);

	bool bcBuildBucket(dtNavMeshQuery* navQuery,const dtMeshTile* tile,dtNavMesh* mesh,bcBucketSet& bset);

	vertex_t bcvertex(vertex_t vparent,vertex_t v,const MaGraph& mg);

	void bcBuildBucketGraph(dtNavMesh* mesh,dtNavMeshQuery* navQuery,bcBucketSet& bcset,BcGraph& cg);
	
	dtStatus bcSetBucketType(dtNavMesh* mesh,bcBucketSet& bcset);
	
	int bcBuildConnectivityGraph(BcGraph& bcg,ConGraph& cg);

	int bcBuildJusGraph(BcGraph& bcg,ConGraph& cg);
	
	dtStatus bcBuildRoads(dtNavMesh* mesh,const dtMeshTile* tile,BcGraph& bcg);

	dtStatus bcfindRoad(dtNavMesh* mesh,dtPolyRef polyRef,dtPolyRef neiRef,unsigned short int* resultAreaId,
					    dtPolyRef* path, int* pathCount, const int maxPath,dtPolyRef* startRef,dtPolyRef* endRef);

	void bcneiiter(dtNavMesh* mesh,dtPolyRef polyRef,dtPolyRef parentRef,dtPolyRef* path,int* j,unsigned short int* areaId);

	bcBucketSet* bcAllocBucketSet();

	int buildDistMat(ConGraph& cg,matrix& distMat);

	int TotalDepth(int node,matrix& distMat,int k,bool t );

	float MeanDepth(int node,matrix& distMat, int k,bool t );

	float Integration(int node,matrix& distMat, int k,bool t );

	int Connectivity(int node,ConGraph& cg);

	int VConnectivity(int node,BcGraph& bcg);

	bool bcLinePolyIntersection(struct bcBucket* L,struct bcBucket* R,polygon_t& poly);

	void bcCalcBetweeness(ConGraph& cg);

	int bcCalcSSParam(ConGraph& cg,BcGraph& bcg,matrix& distMat, int k,polygon_t& poly);

	void bcCalcIntersectionBucParam(BcGraph& cgv);

	int bcCalcSSParamm(ConGraph& cg,BcGraph& bcg,polygon_t& poly);

	int bcBuildVisibilityGraph(BcGraph& cgv,BcGraph& bcg);

	int buildDistMatV(BcGraph& cgv,matrix& distMat);

	void bcCalcSSParamVG(dtNavMesh* mesh,BcGraph& cgv,matrix& distMat, int k);

	void bcCalcSSParammVG(dtNavMesh* mesh,BcGraph& cgv);

	void bcCalcPath(dtNavMesh* mesh,bcBucket* bc,std::map<int, int> mapDensity,dtPolyRef parentRef,unsigned short parentLineRef,float* startPos,float* endPos,
					dtPolyRef endRef,unsigned int fam, unsigned int expl,struct bcLineStep* path,int* j,unsigned short p,unsigned short* t);

	void bcCalcChoice(dtNavMesh* mesh,bcBucket* bc,std::map<int, int> mapDensity,dtPolyRef parentRef,float* startPos,float* endPos,
	                          dtPolyRef endRef,struct bcLineStep* path,int* j,unsigned short p,unsigned short* t);
	
	dtStatus bcFindPathLight(dtNavMesh* mesh,BcGraph& bcg,dtPolyRef parentRef,dtPolyRef startRef, dtPolyRef endRef,float* startPos,
	                              float* endPos,unsigned int fam, unsigned int expl,struct bcLineStep* path, int* pathCount,unsigned short* t);
	void bcCalcPathLight(dtNavMesh* mesh,bcBucket* bc,dtPolyRef parentRef,unsigned short parentLineRef,float* startPos,float* endPos,
	                          dtPolyRef endRef,unsigned int fam, unsigned int expl,struct bcLineStep* path,int* j,unsigned short p,unsigned short* t);

	dtStatus bcFindPath(dtNavMesh* mesh,BcGraph& bcg,std::map<int, int> mapDensity,dtPolyRef parentRef,dtPolyRef startRef, dtPolyRef endRef,
						float* startPos,float* endPos,unsigned int fam, unsigned int expl,struct bcLineStep* path, int* pathCount,unsigned short* t);

	dtStatus bcChoosePath(dtNavMesh* mesh,BcGraph& bcg,std::map<int, int> mapDensity,dtPolyRef parentRef,dtPolyRef startRef, dtPolyRef endRef,float* startPos,
	                              float* endPos,struct bcLineStep* path, int* pathCount,unsigned short* t);

	bool bcBuildRoadClasses(BcGraph& bcg, double pts[], int n);

private:

	void addEdge(MaGraph& graph,const float* v0,const float* v1);

	int getWalls(dtNavMesh* mesh,dtPolyRef queryPointRef,int* e);

	char bcSegSegIntersect2d(const float* a,const float* b,const float* c,const float* d,double& s, double& t);

	char bcVerSegSegIntersect2d(const float* a,const float* b,const float* c,const float* d,float& s, float& t);

	bool bcIntersectLinePoly2D(const float* p0, const float* p1,const float* verts, int nverts,double& smin, double& smax,
							  int& segMin, int& segMax);

	bool bcIntersectVerLinePoly2D(const float* p0, const float* p1,const float* verts, int nverts,
							  float& ymin, float& ymax,int& segMin, int& segMax);

	int getNeighbours(dtNavMesh* mesh,dtPolyRef queryPointRef,dtPolyRef* neighRef);
	
	bool getABPointsRef(dtNavMesh* mesh,dtPolyRef from, dtPolyRef to, float* left, float* right,
				          unsigned char& fromType, unsigned char& toType);
	bool getABPoints(dtPolyRef from, const dtPoly* fromPoly, const dtMeshTile* fromTile,
								dtPolyRef to, const dtPoly* toPoly, const dtMeshTile* toTile,
								float* left, float* right);

	float bcDistancePtPtSqr2D(float* A, float* B);

	float bcCalcArea(const dtMeshTile* tile,const dtPoly* poly);

	float bcCalcRoadArea(dtNavMesh* mesh,dtPolyRef* linePoly,int n);

	float bcCalcAngle(float* p1,float* p2,float* p3);

	float bcCalcdist(float* p1,float* p2,float* p3);

	int bcCalcCost(std::vector<float>& costVect,std::vector<float> angles,std::vector<float> dists,
		           unsigned int fam, unsigned int expl);

	int bcCalcBenefit(std::vector<float>& costVect,std::vector<float> angles,std::vector<float> dists);

	int bcBuildPartialIso(dtNavMesh* mesh,const dtMeshTile* tile,dtPolyRef bcPolyRef,float* verts,int nv,std::list<polygon_t>& output);

	bool bcAffectIsoToNei(dtNavMesh* mesh,const dtMeshTile* tile,float* queryPoint,dtPolyRef bcPolyRef,std::list<polygon_t> polygons,float* a,double* d);

	BcGraph& m_bucketGraph;
	BcGraph& m_visibilityBucketGraph;
	bucket_property_t m_bct_prop;
	bucket_property_t m_bctv_prop;
	line_property_t m_ln_prop;
	static const unsigned short SeuilMax = 3;
	static const unsigned short SeuilMin = 0.3;
	
};

bcBucketMesh* bcAllocBucketMesh(BcGraph& bcg,BcGraph& vbcg);
void bcFreeBucketMesh(bcBucketMesh* bucmesh);

#endif // BUCKETBUILDER_H