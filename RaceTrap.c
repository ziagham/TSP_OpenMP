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
#include "StopWatch.h"


typedef struct {
    int x; 
    int y; 
} Coord; 

typedef struct {
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

FILE *fp;
int DO_DUMP = 0; // true if we want to dump the iterations to the file

void dump_data(RouteDefinition *route) 
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
  for (i = 0; i < nBags; i++) {
      fprintf(fp, "%d ", route->path[i]);
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
RouteDefinition *ShortestRoute(RouteDefinition *route)
{ 
    int i;
    RouteDefinition *bestRoute;
    
    if (route->nPlaced == nBags)
    {
            // If all grain bags have been placed we simply add the distance to get back home (closing the loop)
        route->length += distanceTable[route->path[route->nPlaced-1]][route->path[0]];
        return route;
    }
    
        // Provide an initial subtree that will be replaced with a better
        // solution. We set the length to something that will be longer
        // than any possible path to make sure that it will be replaced with
        // the first possible solution.
    bestRoute = Alloc_RouteDefinition(); 
    bestRoute->length = maxRouteLen;
    
        // Try a number of permutations:
        // For each 'i', use bag 'i' as the next bag in the route (at
        // position route->nPlaced). To do this, we copy 'route' to
        // newRoute and swap the bags in position 'i' and route->nPlaced
        // in newRoute.
    for (i = route->nPlaced; i < nBags; i++)
    {
            // newlength is the current route length + the path to bag number 'i'.
        double newLength = route->length + distanceTable[route->path[route->nPlaced-1]][route->path[i]];
        RouteDefinition *newRoute; 
        
            // Before we set up the new branch, check whether the partial path 
            // that we construct is longer than the current best complete route. 
            // There is no point in branching along that path if it is longer, as 
            // the total path will never be the shortest path. 
        if (newLength >= globalBest)
            continue;  // New path was longer, try another
        
            // The new partial route is shorter than the currently
            // best known complete route, so it's a candidate. We search
            // further for routes starting with this partial route.
        newRoute = Alloc_RouteDefinition();  
        
        memcpy(newRoute->path, route->path, nBags);   // Copy current route from route
        
            // Swaps the position of bag # 'i' and bag # 'nPlaced' from route
        newRoute->path[route->nPlaced] = route->path[i];
        newRoute->path[i]              = route->path[route->nPlaced]; 
        newRoute->nPlaced = route->nPlaced + 1;
        newRoute->length  = newLength;
        
            // Recurse further into this tree and return the shortest path along this permutation
        newRoute = ShortestRoute(newRoute);
        
        if (newRoute->length < bestRoute->length) 
        { 
                // Found a new best route, free the old one and replace with the new one
            free(bestRoute);
            bestRoute = newRoute; 
            if (bestRoute->length < globalBest)
            {
                globalBest = bestRoute->length;
                dump_data(bestRoute);
            }
        }
        else 
        {
                // The new route was longer than the current best solution. 
            free(newRoute);
        }
    }
    free(route);
    
    return bestRoute;
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
    for (i = 0; i < nBags; i++)
        distanceTable[i] = (double*) malloc(nBags * sizeof(double));
    
        // Compute the distances between each of the grain bags.
    for (i = 0; i < nBags; i++)	  
        for (j = 0; j < nBags; j++)	  
            distanceTable[i][j] = EuclidDist(&bagCoords[i], &bagCoords[j]);
}


int main (int argc, char **argv) 
{
    RouteDefinition *originalRoute, *res;
    int i;
    char buf[256];
    
    // ./RaceTrap dump
    if (argc == 2 && strcmp("dump", argv[1]) == 0) {
        DO_DUMP = 1;
    }

    ReadRoute();
    
        // Set up an initial path that goes through each bag in turn. 
    originalRoute = Alloc_RouteDefinition(); 
    for (i = 0; i < nBags; i++)
        originalRoute->path[i] = (unsigned char) i;
    
    originalRoute->length = 0.0;
    originalRoute->nPlaced = 1;
    dump_data(originalRoute);

    sw_init();
    sw_start();
        // Find the best route
    res = ShortestRoute(originalRoute);  
    sw_stop();
    sw_timeString(buf);
    
    printf("Route length is %lf it took %s\n", res->length, buf);
    
    if (DO_DUMP) {
        fclose(fp);
    }

    return 0;
}
