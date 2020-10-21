/*
 * RaceTrap implementation based on RaceTrap.java
 *
 * Created on 22. juni 2000, 13:48
 * 
 * Brian Vinter
 * 
 * Modified by John Markus Bjørndalen, 2008-12-04, 2009-10-15. 
 * Modified by Sergiusz Michalik, 2018-09-20
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <limits.h>
#include "StopWatch.h"
#include <omp.h>

typedef struct {
    int x; 
    int y; 
} Coord; 

typedef struct {
    int           visited[0]; // visited array for each thread.
    double        length;     // Length of the current path (distance)
    unsigned char nPlaced;    // Number of bags currently placed
    unsigned char path[0];    // Array of vertex/bag numbers in the path (see comment in
                              // Alloc_RouteDefinition())
} RouteDefinition; 

int      nBags = 0;             // Number of grain-bags
Coord   *bagCoords;             // Coordinates for the grain-bags
double **distanceTable;         // Table of distances between any two grain-bags
double   maxRouteLen = 10E100;  // Initial best distance, must be longer than any possible route
double   globalBest  = 10E100;  // Bounding variable
int     *bestPath;              // Array of vertex/bag numbers in the path

FILE *fp;
int DO_DUMP = 0; // true if we want to dump the iterations to the file

void dump_data(int *path) 
{
  int i;
  if (!DO_DUMP)
    return;

  if (fp == NULL) {
      remove("data/racetrap.data");
    fp = fopen("data/racetrap.data", "a+");
  }
  /* Stores the data as a Python datastructure for easy inspection and plotting
   */
  printf("dumping data\n");
  for (i = 0; i<nBags; i++) {
      fprintf(fp, "%d ", path[i]);
  }
  fprintf(fp, "\n");
}

static inline RouteDefinition* Alloc_RouteDefinition()
{
    if (nBags <= 0) 
    {
        fprintf(stderr, "Error: Alloc_RouteDefinition called with invalid nBags (%d)\n", nBags); 
        exit(-1); 
    }
    // NB: The +nBags*sizeof.. trick "expands" the path[0] array in RouteDefintion
    // to a path[nBags] array.
    RouteDefinition *def = NULL;  

    return (RouteDefinition*) malloc(sizeof(RouteDefinition) + nBags * sizeof(def->path[0]));
}

// copy the optimal path to the bestPath array and used it to save in dump_data if needed.
void copyToBestPath(int currrent_path[]) {
    int i;
    for (i=0; i<nBags; i++)
        bestPath[i] = currrent_path[i];
}

//calculate the cost of first tour edges adjacent to a city
double cost_of_first_edge(int city) {
    double min = maxRouteLen;
    int k;
    for (k=0; k<nBags; k++){
        if (distanceTable[city][k]<min && city != k){
            min = distanceTable[city][k];
        }
    }
    return min;
}

//calculate the cost of second tour edges adjacent to a city
double cost_of_second_edge(int city) {
    double first = maxRouteLen, second = maxRouteLen;
    int j;
    for (j=0; j<nBags; j++) {
        if (city == j)
            continue;
 
        if (distanceTable[city][j] <= first) {
            second = first;
            first = distanceTable[city][j];
        }
        else if (distanceTable[city][j] <= second && distanceTable[city][j] != first){
            second = distanceTable[city][j];
        }
    }
    return second;
}

//calculte the initial lowerBound for start node
double calculate_initial_lowerBound(){
    int i;
    double result = 0.0;

    // use reduction in order to incrementing a result variable in parallel mode
    #pragma omp parallel for reduction(+:result)
    for (i=0; i<nBags; i++){
        result += (cost_of_first_edge(i) + cost_of_second_edge(i));
    }
    result =  result/2;

    return result;
}

//calculte the initial lowerBound for the first level
double calculate_firstLevel_lowerBound(int firstPath, int secondPath){
    int i;
    double result = 0.0;
    result = ((cost_of_first_edge(secondPath) + cost_of_first_edge(firstPath))/2);
    return result;
}

//calculte the lowerBound
double calculate_lowerBound(int firstPath, int secondPath){
    int i;
    double result = 0.0;
    result = ((cost_of_second_edge(secondPath) + cost_of_first_edge(firstPath))/2);
    return result;
}

/* 
 * A recursive Traveling Salesman Solver using branch-and-bound. 
 * 
 * Returns the shortest roundtrip path based on the starting path described in "route". 
 * 
 * The returned path is "route" if "route" already includes all the bags
 * in the route. If not, route will be freed, and a new path will be returned. 
 * 
 * NB: this function is destructive - it calls free() on route if it finds
 *     a better route than the provided route! 
 * 
 * NB2: there is a slight problem with the below code: ShortestRoute will return a 
 *      semi-intialized bestRoute if all the new permitations are longer than 
 *      globalBest. It shouldn't cause any problems though, as that route 
 *      will be thrown away. 
 */ 
void ShortestRoute(double lowerBound, double length, int level, int route_path[], int visited_city[]) {
    if (level==nBags) {
        if (distanceTable[route_path[level-1]][route_path[0]] != 0) {
            double currLength = length + distanceTable[route_path[level-1]][route_path[0]];
            
            //critial area that any thread need to accesss it.
            #pragma omp critical
            if (currLength < globalBest) {
                copyToBestPath(route_path);
                globalBest = currLength;
                dump_data(route_path);
            }
        }
        return;
    }
 
    int i;
    #pragma omp parallel for default(none) firstprivate(lowerBound, length, level)\
        shared(route_path,visited_city,distanceTable,globalBest,nBags)
    for (i=0; i<nBags; i++) {
        int local_visited_city[nBags];
        int local_route_path[nBags];
        int j;
        for (j=0;j<nBags;j++) {local_visited_city[j] = visited_city[j];}
        for (j=0;j<nBags;j++) local_route_path[j] = route_path[j];
        if (distanceTable[local_route_path[level-1]][i] != 0 && local_visited_city[i] == 0) {
            double tempBound = lowerBound;
            length += distanceTable[local_route_path[level-1]][i];
            if (level==1)
              lowerBound -= calculate_firstLevel_lowerBound(i, local_route_path[level-1]);
            else
              lowerBound -= calculate_lowerBound(i,local_route_path[level-1]);
            if (lowerBound + length < globalBest) {
                local_route_path[level] = i;
                local_visited_city[i] = 1;
                ShortestRoute(lowerBound, length, level+1, local_route_path,local_visited_city);
            }
            length -= distanceTable[local_route_path[level-1]][i];
            lowerBound = tempBound;
        }
    }
}

// In the desert, the shortest route is a straight line :)
double EuclidDist(Coord *from, Coord *to)
{ 
    double dx = fabs(from->x - to->x);
    double dy = fabs(from->y - to->y);
    return sqrt(dx*dx + dy*dy);
}

// Reads coordinates from a file and generates a distance-table
static void ReadRoute()
{ 
    FILE *file = fopen("./route.dat", "r");
    int i,j;
    
    // Read how many bags there are
    if (fscanf(file, "%d", &nBags) != 1) 
    {
        printf("Error: couldn't read number of bags from route definition file.\n");
        exit(-1);
    }
    
    // Allocate array of bag coords. 
    bagCoords = (Coord*) malloc(nBags * sizeof(Coord)); 
    
    // Read the coordinates of each grain bag
    #pragma omp parallel for
    for (i = 0; i < nBags; i++)
    {
        if (fscanf(file,"%d %d", &bagCoords[i].x, &bagCoords[i].y) != 2) 
        {
            printf("Error: missing or invalid definition of coordinate %d.\n", i);
            exit(-1);
        }
    }
    
    // Allocate distance table 
    distanceTable = (double**) malloc(nBags * sizeof(double*));
    #pragma omp parallel for
    for (i = 0; i < nBags; i++)
        distanceTable[i] = (double*) malloc(nBags * sizeof(double));
    
    // Compute the distances between each of the grain bags.
    for (i = 0; i < nBags; i++)	  
        for (j = 0; j < nBags; j++)	  
            distanceTable[i][j] = EuclidDist(&bagCoords[i], &bagCoords[j]);
}

int main (int argc, char **argv) 
{
    //get the number of threads to run the parallel program
    int numThreads = strtol(argv[1], NULL, 10);

    //set the OpenMP parameters in order to set the global number of threads in parallel for clause
    omp_set_dynamic(0);
    omp_set_num_threads(numThreads);

    char buf[256];

    // ./RaceTrap dump
    if (argc == 3 && strcmp("dump", argv[2]) == 0) {
        DO_DUMP = 1;
    }

    ReadRoute();

    int currrent_path[nBags];
    double initial_lowerBound = 0.0;
    int visited[nBags];

    bestPath = (int*)malloc(nBags*sizeof(int));
    memset(currrent_path, -1,sizeof(currrent_path));
    memset(visited, 0, sizeof(visited));

    initial_lowerBound = calculate_initial_lowerBound();
    visited[0] = 1;
    currrent_path[0] = 0;

    sw_init();
    sw_start();

    ShortestRoute(initial_lowerBound, 0, 1, currrent_path,visited);

    sw_stop();
    sw_timeString(buf);
    
    printf("Route length is %lf it took %s\n", globalBest, buf);

    if (DO_DUMP) {
        fclose(fp);
    }

    return 0;
}