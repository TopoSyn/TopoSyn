
#include "DetourDebugDraw.h"
#include "VoroDiag.h"
#include "DetourNavMesh.h"
#include "BucketAlloc.h"
#include <new>

#ifdef WIN32
#	define snprintf _snprintf
#endif

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

void calcDirection(int compIndex,compMap_t compMap,const Graph g,coor_property_t coor,vertex_t V1,vertex_t& V2)
{
	float* B;float* C;float* A;
	vertex_t rv,cVert,aVert;
	std::map<vertex_t,size_t>::iterator itc;
	std::map<vertex_t,float*>::iterator itg;
	bool top;

	for(itc=compMap.begin(); itc!=compMap.end(); ++itc)
	   {
            if (itc->second==compIndex)
			{
				rv=itc->first;
				break;
			}
	   }

	B=coor[rv]; 
	edge_iterator_t ei, ei_end;

    for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
	{
        if (source(*ei, g)==rv)
		{
			cVert=target(*ei, g);
			C=coor[target(*ei, g)];
        }
		if (target(*ei, g)==rv)
		{
			aVert=source(*ei, g);
			A=coor[source(*ei, g)];
		}
	}
	
	top=calcVectProd(A,B,C);
    if (top)
	{
		V2=cVert;
	}
	else
	{
		V2=aVert;
	}
		
}

int bcBuildGraphContours(const dtMeshTile* tile,const dtNavMesh& mesh,Graph& g,coor_property_t& coor)
{
	
    boost::associative_property_map<coor_property_t> coor_map(coor);
	vertex_t vg,vg0,vg1;
	float* fx;
	bool exist;
	std::pair<vertex_iterator_t, vertex_iterator_t> vp;

	for (int i = 0; i < tile->header->polyCount; ++i)
	{
		const dtPoly* p = &tile->polys[i];
		if (p->type != DT_POLYTYPE_CITY_GROUND) continue;
		
		for (int j = 0, nj = (int)p->vertCount; j < nj; ++j)
		{
			if (p->neis[j] != 0) continue;
			exist= false;
		    float* v0 = &tile->verts[p->verts[j]*3];
            for (vp = vertices(g); vp.first != vp.second; ++vp.first)
	        {
		      fx=coor_map[*vp.first];
       
		      if (fx==v0)
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
				vg0 = boost::add_vertex(g);
				put(coor_map,vg0,v0);
			}


			exist= false;
			float* v1 = &tile->verts[p->verts[(j+1)%nj]*3];
			for (vp = vertices(g); vp.first != vp.second; ++vp.first)
	        {
		      fx=coor_map[*vp.first];
       
		      if (fx==v1)
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
				vg1 = boost::add_vertex(g);
				put(coor_map,vg1,v1 );
			}

			boost::add_edge(vg0, vg1, g);
			boost::add_edge(vg1, vg0, g);
			
		}
	}
int numV=boost::num_vertices(g);

return numV;

}

void bcExtractContours(const Graph g,compMap_t& compMap)
{
    
    boost::associative_property_map<compMap_t> componentMap(compMap); 
	int num = boost::connected_components(g, componentMap); 
	
}

long int bcBuildBorderContour(const Graph g,coor_property_t coor,compMap_t compMap,const float* bmin,const float* bmax,
	                        float walkableRadius,contourVector_t& borderContour)
{
	struct bcContour borderCont;
	long int numbElem=0,idx=-1;
	float* fx=NULL;
	vertex_t rv;
	vertex_t borderVert,V1,V2;
	std::map<vertex_t,size_t>::iterator itc;
	std::map<vertex_t,float*>::iterator itg;
	bool exist=false;
	std::pair<vertex_iterator_t, vertex_iterator_t> vp;
	
	static const float EPS = 0.00000001f;

	for (itg=coor.begin(); itg!=coor.end(); ++itg)
	{
		fx = itg->second;

		if ((ceilf(fx[0])-ceilf(bmin[0]+walkableRadius)<EPS)||
			(ceilf(fx[0])-ceilf(bmax[0]-walkableRadius)<EPS)||
			(ceilf(fx[2])-ceilf(bmin[2]+walkableRadius)<EPS)||
			(ceilf(fx[2])-ceilf(bmax[2]-walkableRadius)<EPS))
	        {
		      rv=itg->first;
			  for(itc=compMap.begin(); itc!=compMap.end(); ++itc)
			  {
				  if (itc->first==rv )
				  {
					idx=itc->second;
					exist=true;
					break;
				  }
			  }

			}
		if (exist)
		{
			break;
		}
	}
	
	if (exist==false)
	{
		printf("bcBuildContours: Failed to build border contour.");
		return -1;
	}
	else
	{
	    
		for(itc=compMap.begin(); itc!=compMap.end(); ++itc)
		{
			if (itc->second==idx)
				numbElem=numbElem+1;
		}

		V1=rv;
		calcDirection(idx,compMap,g,coor,V1,V2);
		std::map<vertex_t,float*>::iterator itg;
	    edge_iterator_t ei, ei_end;
		float* A=NULL;
		float* B=NULL;

		for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
	    {
			if ((source(*ei, g)==V1)&&(target(*ei, g)==V2))
			{
				A=coor[source(*ei, g)];
				B=coor[target(*ei, g)];
				borderContour.push_back(A);
				borderContour.push_back(B);
				break;
			}
	    }

		for (int j=2;j<numbElem;++j)
		{
			for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
			{
				if ((source(*ei, g)==V2)&&(target(*ei, g)!=V1))
				{
					B=coor[target(*ei, g)];
					borderContour.push_back(B);
					V1=source(*ei, g);
					V2=target(*ei, g);
					break;
				}
			}
		}
		return idx;
	}
}

void bcBuildWallsPoly (contourVector_t& borderContour,std::vector<contourVector_t>& holesContours,polygon_t& poly)
{
   static const float thr = 0.01f*0.01f;
   float* A=NULL;
   float* B=NULL;
   float* C=NULL;
   
   //****************create border contour poly points***********************
   std::vector<point_t> borderPoints;
   A=borderContour[0];
   point_t point1(A[0], A[2]);
   borderPoints.push_back(point1);
   for (int i=1;i<borderContour.size();++i)
    {
	  B=borderContour[i];
	  C=borderContour[(i+1)%borderContour.size()];
	  if (distPtLine2d(C,A,B) < thr)
	   {
		 continue;
	   }
	  else
	   {
           point_t point(B[0], B[2]);
		   borderPoints.push_back(point);
		   A=B;
	   }
    }
   A=borderContour[0];
   point_t point2(A[0], A[2]);
   borderPoints.push_back(point2);

   boost::geometry::append(poly, borderPoints);
   boost::geometry::correct(poly);

//*******************************create inner polygons************************
   std::vector<point_t> inner_points;
   for (int k=0;k<holesContours.size();++k)
	{
	   A=holesContours[k][0];
	   point_t point1(A[0], A[2]);
	   inner_points.push_back(point1);
	   for (int i=1;i<holesContours[k].size();++i)
		{
		  B=holesContours[k][i];
		  C=holesContours[k][(i+1)%holesContours[k].size()];
		  if (distPtLine2d(C,A,B) < thr)
		   {
			 continue;
		   }
		  else
		   {
			   point_t point(B[0], B[2]);
			   inner_points.push_back(point);
			   A=B;
		   }
		 }

	   A=holesContours[k][0];
	   point_t point2(A[0], A[2]);
	   inner_points.push_back(point2);
	   polygon_t hull;
	   polygon_t con_hull;
	   boost::geometry::append(hull, inner_points);
       boost::geometry::convex_hull(hull, con_hull);
	   inner_points.clear();

	   for(auto it = boost::begin(boost::geometry::exterior_ring(con_hull)); 
		        it != boost::end(boost::geometry::exterior_ring(con_hull));++it)
	   {
			inner_points.push_back(*it);
	   }
	  
	   poly.inners().resize(k+1);
       boost::geometry::model::ring<point_t>& inner = poly.inners().back();
       boost::geometry::assign_points(inner, inner_points);
	   inner_points.clear();
	   
    }
   boost::geometry::correct(poly);
   int a=boost::geometry::area(poly); 
  
}

int bcBuilHolesContours(const Graph g,coor_property_t coor,compMap_t compMap,long int borderContourIdx,
	                    std::vector<contourVector_t>& holesContours)
{
	struct bcContour holeCont;
	int holesContoursSize=0;
	int numbElem;
	bool exist;
	std::vector<int> visitedIdx;
    visitedIdx.push_back(borderContourIdx);
	vertex_t borderVert,V1,V2;
	std::map<vertex_t,size_t>::iterator itc;
	std::map<vertex_t,float*>::iterator itg;
	std::map<vertex_t,size_t>::iterator itcS;
	edge_iterator_t ei, ei_end;
    float* A=NULL;
    float* B=NULL;

	for(itc=compMap.begin(); itc!=compMap.end(); ++itc)
	{
			exist=false;
		    for(int y=0; y<visitedIdx.size(); y++)
            {
                if (visitedIdx[y]==itc->second)
				{
				exist=true;  
				break;
				}
            }
			if (exist)
			{
				continue;
			}
			else
			{
				visitedIdx.push_back(itc->second);
				holesContoursSize=holesContoursSize+1;
				for (itg=coor.begin(); itg!=coor.end(); ++itg)
	            {
		          if (itg->first==itc->first)
				  {
					V1 = itg->first;
				  }
				}
		        calcDirection(itc->second,compMap,g,coor,V1,V2);
		        holeCont.index=itc->second;
				numbElem=0;
				for(itcS=compMap.begin(); itcS!=compMap.end(); ++itcS)
		        {
			        if (itcS->second==itc->second)
				        numbElem=numbElem+1;
		        }
	
                contourVector_t holeContour;
               
	            for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
	            {
                       if ((source(*ei, g)==V1)&&(target(*ei, g)==V2))
		               {
			              A=coor[source(*ei, g)];
			              B=coor[target(*ei, g)];
			              holeContour.push_back(A);
			              holeContour.push_back(B);
			              break;
                       }
	             }
	            for (int j=2;j<numbElem;++j)
	            {
		             for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
		             {
			              if ((source(*ei, g)==V2)&&(target(*ei, g)!=V1))
			              {
				              B=coor[target(*ei, g)];
				              holeContour.push_back(B);
				              V1=source(*ei, g);
	                          V2=target(*ei, g);
				              break;
			               }
		             }
	            }
                holesContours.push_back(holeContour);
			}
    }
return holesContours.size();
}

ma_vertex_t bcneig(ma_vertex_t vparent,ma_vertex_t v,const MaGraph& mg)
{
	int i=out_degree(v, mg);
	if (i<3)
		{
			if (i==1)
			{
				return v; 
			}
			if (i==2)
			{
				std::pair<ma_adjacency_iterator_t, ma_adjacency_iterator_t> neighbors = boost::adjacent_vertices(v,mg);
				for(; neighbors.first != neighbors.second; ++neighbors.first)
				{
					if (*neighbors.first!=vparent)
					{
						return bcneig(v,*neighbors.first,mg);
					}
				}
			}
		}
	else
		{
		return v;
		}
}
 
