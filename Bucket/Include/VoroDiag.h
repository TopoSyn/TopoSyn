
//
#ifndef VORODIAG_H
#define VORODIAG_H

#include <cstdio>
#include <vector>




#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/polygon.hpp>

#include <boost/polygon/voronoi.hpp>
#include "boost/polygon/segment_data.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#include <iostream>
#include <vector>
#include "DetourNavMesh.h"

struct VertexProperties 
{ 
   double x;
   double y;
   VertexProperties() : x(0.0), y(0.0) {}
   VertexProperties(double m, double n) : x(m), y(n) {}
};

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor vertex_t;
typedef boost::graph_traits<Graph>::edge_descriptor edge_t;
typedef boost::graph_traits<Graph>::vertex_iterator vertex_iterator_t;
typedef boost::graph_traits<Graph>::edge_iterator edge_iterator_t;
typedef boost::graph_traits<Graph>::adjacency_iterator adjacency_iterator_t;


typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,VertexProperties> MaGraph;
typedef boost::graph_traits<MaGraph>::vertex_descriptor ma_vertex_t;
typedef boost::graph_traits<MaGraph>::edge_descriptor ma_edge_t;
typedef boost::graph_traits<MaGraph>::vertex_iterator ma_vertex_iterator_t;
typedef boost::graph_traits<MaGraph>::edge_iterator ma_edge_iterator_t;
typedef boost::graph_traits<MaGraph>::out_edge_iterator ma_out_edge_iterator_t;
typedef boost::graph_traits < MaGraph >::adjacency_iterator ma_adjacency_iterator_t;


typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,boost::no_property,boost::property<boost::edge_weight_t,int>> BcGraph;
typedef boost::graph_traits<BcGraph>::vertex_descriptor b_vertex_t;
typedef boost::graph_traits<BcGraph>::edge_descriptor b_edge_t;
typedef boost::graph_traits<BcGraph>::vertex_iterator b_vertex_iterator_t;
typedef boost::graph_traits<BcGraph>::edge_iterator b_edge_iterator_t;
typedef boost::graph_traits<BcGraph>::out_edge_iterator b_out_edge_iterator_t;
typedef boost::graph_traits<BcGraph>::in_edge_iterator b_in_edge_iterator_t;
typedef boost::graph_traits<BcGraph>::adjacency_iterator b_adjacency_iterator_t;

typedef struct ConVertexProperties 
{ 
   b_vertex_t source;
   b_vertex_t target;
   ConVertexProperties() : source(0),target(0) {}
   ConVertexProperties(b_vertex_t s,b_vertex_t t) : source(s),target(t) {}
};

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,ConVertexProperties,boost::property<boost::edge_weight_t,int>> ConGraph;
typedef boost::graph_traits<ConGraph>::vertex_descriptor c_vertex_t;
typedef boost::graph_traits<ConGraph>::edge_descriptor c_edge_t;
typedef boost::graph_traits<ConGraph>::vertex_iterator c_vertex_iterator_t;
typedef boost::graph_traits<ConGraph>::edge_iterator c_edge_iterator_t;
typedef boost::graph_traits<ConGraph>::adjacency_iterator c_adjacency_iterator_t;



typedef std::map<vertex_t,float*> coor_property_t;
typedef std::map<vertex_t,size_t> compMap_t;
typedef std::vector <int> vectcont_t;
typedef boost::graph_traits<Graph>::edge_iterator edge_iter_t;
typedef std::vector <float*> contourVector_t;



typedef boost::geometry::model::d2::point_xy<double> point_t;
typedef boost::geometry::model::polygon<point_t> polygon_t;
typedef boost::geometry::model::linestring<point_t> line_t;
typedef boost::geometry::model::segment<point_t> segment_t;
typedef std::vector<point_t> pointsVect_t;

struct bcContour
{
	int index;
	int contSize;
	vertex_t firstVert;
	vertex_t SecondVert;
};


using boost::polygon::voronoi_diagram;
typedef double CoordinateType;
typedef boost::polygon::point_data< CoordinateType > PointType;
typedef boost::polygon::segment_data< CoordinateType > SegmentType;
typedef boost::polygon::rectangle_data<CoordinateType> rect_type;
typedef boost::polygon::voronoi_diagram< CoordinateType > VoronoiDiagram;


void calcDirection(int compIndex,compMap_t compMap,const Graph g,coor_property_t m_coor,vertex_t V1,vertex_t& V2);


int bcBuildGraphContours(const dtMeshTile* tile,const dtNavMesh& mesh,Graph& g,coor_property_t& coor);
long int bcBuildBorderContour(const Graph g,coor_property_t coor,compMap_t compMap,const float* bmin,const float* bmax,
	                        float walkableRadius,contourVector_t& borderContour);
int bcBuilHolesContours(const Graph g,coor_property_t coor,compMap_t compMap,long int borderContourIdx,std::vector<contourVector_t>& holesContours);
void bcExtractContours(const Graph g,compMap_t& compMap);
void bcBuildWallsPoly (contourVector_t& borderContour,std::vector<contourVector_t>& holesContours,polygon_t& poly);

ma_vertex_t bcneig(ma_vertex_t vparent,ma_vertex_t v,const MaGraph& mg);


static const std::size_t EXTERNAL_COLOR = 1;
static PointType shift_;
static std::vector<PointType> point_data_;
static std::vector<SegmentType> segment_data_;
static SegmentType Segment;
static rect_type brect_;
static bool primary_edges_only;
static bool internal_edges_only;
static bool brect_initialized_;
 
#endif // VORODIAG_H

