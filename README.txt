1.	Toposyn Project description:
The software is a software application which allows generating a representation of the navigable areas in an urban environment. This representation enables realistic simulation of pathfinding behaviour in a familiar environment through the software. Simulating pathfinding in an unfamiliar environment is our next topic of research.

2.	Features of Toposyn
Toposyn creates a data structure that represents the walkable surfaces in virtual urban environments. We refer to this structure as a bucket mesh. This bucket mesh is used to endow agents with mental maps and simulate the pathfinding behaviour.
This code is based on the Original Recast Demo framework. Recast is a navigation mesh construction toolset for games. Recast is accompanied by Detour and Crowd Detour. Detour is a spatial reasoning toolkit and Crowd Detour is a crowd management module. 

3.	Building Toposyn
This program is developed in Microsoft Visual Studio 2010. To compile the project, you need steps follow:
•	Install MPFR and GMP for C++ on visual studio 2010.
•	Install Boost 1.60 so that the path ‘local/boost_1_60_0’ is valid. 
•	Download SDL development library (version 1.2.15 or later) and unzip it in ‘TopoSynDemo\Contrib’. Rename the SDL folder such that the path ‘TopoSynDemo\Contrib\SDL\lib\x86’ is valid.
•	Open the solution, build, and run.

4.	Getting started with TopoSyn
There are three tasks that you are very likely to encounter when using TopoSyn: (a) Building a navigation mesh, (b) building a bucket mesh, and (c) testing the created bucket mesh or run a crowd simulation. This section explains how to perform these three tasks:
a)	Building navigation mesh:
It consists in creating a navigation mesh which will be used to create the bucket mesh. It is important to understand the various configuration parameters that affect the final result:

Configuration Parameter	
Cell Size: The width and depth resolution used when sampling the source geometry. The width and depth of the cell columns that make up voxel fields.
Cell Height: The height resolution used when sampling the source geometry. The height of the voxels in voxel fields.
Max Climb: Represents the maximum ledge height that is considered to still be traversable.
Max Slope: The maximum slope that is considered traversable. (In degrees.)
Min Region Size: The minimum region size for unconnected (island) regions.
Merged Region Size: Any regions smaller than this size will, if possible, be merged with larger regions.
Max Edge Length: The maximum length of polygon edges that represent the border of meshes.
Max Edge Error: The maximum distance the edges of meshes may deviate from the source geometry.
Verts Per Poly: The maximum number of vertices per polygon.
Sample Distance: Sets the sampling distance to use when matching the detail mesh to the surface of the original geometry.
Max Sample Error: The maximum distance the surface of the detail mesh may deviate from the surface of the original geometry.

b)	Build bucket mesh
This step takes the navigation mesh created in the previous step as input and generates the bucket mesh which represents the navigable surfaces of the environment by a set of polygons. These polygons are of two types: crossroad buckets or road buckets.

c)	Testing the created bucket mesh or run a crowd simulation
Once the environment description (buckets mesh) is generated two tools are available: 
•	Testing the buckets mesh: calculating paths according to the shortest path algorithm or the proposed incremental algorithm (implementation of the strategy applied by pedestrians familiar with their environment).
•	Running a crowd simulation: simulating the pathfinding behaviour of pedestrians according to their level of familiarity with their environment.

5.	Special Thanks To
•	Mikko Mononen for making Recast available for study and use by everyone. More information at: https://github.com/recastnavigation/recastnavigation
•	John McCullock for its implementation in C++ of the K-Means clustering algorithm. More information at: http://mnemstudio.org/clustering-k-means-example-1.htm
