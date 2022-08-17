// application.cpp <Starter Code>
// <Abeer Fatima>
//
// University of Illinois at Chicago
// CS 251: Fall 2021
// Project #7 - Openstreet Maps
//
// References:
// TinyXML: https://github.com/leethomason/tinyxml2
// OpenStreetMap: https://www.openstreetmap.org
// OpenStreetMap docs:
//   https://wiki.openstreetmap.org/wiki/Main_Page
//   https://wiki.openstreetmap.org/wiki/Map_Features
//   https://wiki.openstreetmap.org/wiki/Node
//   https://wiki.openstreetmap.org/wiki/Way
//   https://wiki.openstreetmap.org/wiki/Relation
//

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iomanip> /*setprecision*/
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <stack>
#include <string>
#include <vector>

#include "dist.h"
#include "graph.h"
#include "osm.h"
#include "tinyxml2.h"

using namespace std;
using namespace tinyxml2;

// global variable for infinity
const double INF = numeric_limits<double>::max();

//
// Implement your creative component application here
// TO DO: add arguments
//
void creative() {}

// PRIORITIZE CLASS
class prioritize {
 public:
  bool operator()(const pair<long long, double>& p1,
                  const pair<long long, double>& p2) const {
    return p1.second > p2.second;
  }
};

// Helper function to  print the nodes
void PrintNodes(long long near1, long long near2, long long centerNear,
                map<long long, Coordinates>& Nodes) {
  cout << endl;
  cout << "Nearest P1 node:" << endl;
  cout << " " << near1 << endl;
  cout << " (" << Nodes[near1].Lat << ", " << Nodes[near1].Lon << ")" << endl;

  cout << "Nearest P2 node:" << endl;
  cout << " " << near2 << endl;
  cout << " (" << Nodes[near2].Lat << ", " << Nodes[near2].Lon << ")" << endl;

  cout << "Nearest destination node:" << endl;
  cout << " " << centerNear << endl;
  cout << " (" << Nodes[centerNear].Lat << ", " << Nodes[centerNear].Lon << ")"
       << endl;
  cout << endl;
}

// helper function to print the buildings
void PrintBuilding(BuildingInfo building1, BuildingInfo building2,
                   BuildingInfo buildingCenter) {
  cout << endl;
  cout << "Person 1's point:" << endl;
  cout << " " << building1.Fullname << endl;
  cout << " (" << building1.Coords.Lat << ", " << building1.Coords.Lon << ")"
       << endl;
  cout << "Person 2's point:" << endl;
  cout << " " << building2.Fullname << endl;
  cout << " (" << building2.Coords.Lat << ", " << building2.Coords.Lon << ")"
       << endl;
  cout << "Destination Building:" << endl;
  cout << " " << buildingCenter.Fullname << endl;
  cout << " (" << buildingCenter.Coords.Lat << ", " << buildingCenter.Coords.Lon
       << ")" << endl;
}

// helper function
// takes in the string and searches abbreviations and returns
// MILESTONE 7
BuildingInfo searchBuilding(string query, vector<BuildingInfo>& Buildings) {
  // MILESTONE 7. search the buildings
  for (auto it : Buildings) {
    if (it.Abbrev == query ||
        it.Fullname == query) {  // compare fullname and abbreviation
      return it;
    }
    // this is the partial searching edge case
    if (it.Abbrev.find(query) != string::npos ||
        it.Fullname.find(query) != string::npos) {
      return it;
    }
  }

  BuildingInfo exit;  // blank building to be returned
  exit.Fullname = "empty";
  exit.Abbrev = "empty";
  return exit;
}

BuildingInfo ModifiedNearest(Coordinates midpoint,
                             vector<BuildingInfo>& Buildings,
                             set<string>& unvisitedB) {
  BuildingInfo near;
  // vector<BuildingInfo> &blackListed
  double min = INF;
  // iterate through the buildings vector
  for (auto it : Buildings) {
    if (unvisitedB.count(it.Fullname) > 0) {
      continue;
    } else {
      double distance = distBetween2Points(midpoint.Lat, midpoint.Lon,
                                           it.Coords.Lat, it.Coords.Lon);
      if (distance < min) {  // if the un
        min = distance;
        near = it;
      }
    }
  }
  return near;
}

// MILESTONE 8
BuildingInfo nearestBuilding(Coordinates midpoint,
                             vector<BuildingInfo>& Buildings) {
  BuildingInfo near;
  // vector<BuildingInfo> &blackListed
  double min = INF;
  // iterate through the buildings vector
  for (auto it : Buildings) {
    double distance = distBetween2Points(midpoint.Lat, midpoint.Lon,
                                         it.Coords.Lat, it.Coords.Lon);
    if (distance < min) {
      min = distance;
      near = it;
    }
  }

  // Do something about the unreachable buildings
  return near;
}

// Milestone 9
// return the nearest node
long long nearestNode(BuildingInfo b, vector<FootwayInfo>& Footways,
                      map<long long, Coordinates>& Nodes) {
  double min = INF;
  long long nearest;
  // iterate through the footways vector
  for (auto iter : Footways) {
    for (size_t i = 0; i < iter.Nodes.size(); i++) {
      long long NodeNear = iter.Nodes[i];
      // calculate distance
      double distance = distBetween2Points(
          Nodes[NodeNear].Lat, Nodes[NodeNear].Lon, b.Coords.Lat, b.Coords.Lon);
      // implement the min algorithm
      if (distance < min) {
        min = distance;
        nearest = NodeNear;
      }
    }
  }
  return nearest;  // return the nearest node here.
}

// write getPath, pass in the predecessor and the end point

// DIJKSTRA'S ALGORITHM FROM THE LAB 11
void Dijkstra(graph<long long, double>& G, long long startV,
              map<long long, double>& distances,
              map<long long, long long>& predecessor) {
  vector<long long> visited;

  // create the queue here
  priority_queue<pair<long long, double>, vector<pair<long long, double>>,
                 prioritize>
      pq;

  vector<long long> allNodes = G.getVertices();
  for (auto vertex : allNodes) {
    distances[vertex] = INF;
    // set the predecessor to the map
    predecessor[vertex] = 0;
    // push to the queue
    pq.push(make_pair(vertex, INF));
  }
  // change all strings to long long, and change ints to doubles
  pq.push(make_pair(startV, 0));
  // set start v to zero
  distances[startV] = 0;
  // set the predecessor to zero too
  predecessor[startV] = 0;

  while (pq.size() > 0) {
    pair<long long, double> thisNode = pq.top();
    pq.pop();
    //     if(thisNode.second == INF){
    //       break;
    //     }
    if (distances[thisNode.first] == INF) {
      break;
    }

    bool found = false;
    for (auto node : visited) {
      if (thisNode.first == node) {
        found = true;
      }
    }
    if (found) {
      continue;
    }
    visited.push_back(thisNode.first);
    set<long long> neighbors = G.neighbors(thisNode.first);
    for (auto neighbor : neighbors) {
      double edgeWeight = 0.0;
      G.getWeight(thisNode.first, neighbor, edgeWeight);

      double pathDist = thisNode.second + edgeWeight;

      if (pathDist < distances[neighbor]) {
        distances[neighbor] = pathDist;
        // set predecessor to current vertex
        predecessor[neighbor] = thisNode.first;  // the
        pq.push(make_pair(neighbor, pathDist));  // push to the queue here
      }
    }
  }
  // return visited; // don't return since it's a void function
}

// Milestone 11 helper function
vector<long long> getPath(map<long long, long long>& predecessor,
                          long long destination) {
  vector<long long> result;

  long long curV = destination;
  stack<long long> S;
  while (curV != 0) {
    S.push(curV);  // push to the stack
    curV = predecessor[curV];
  }
  // iterate and push to vector
  while (!S.empty()) {
    curV = S.top();
    S.pop();
    result.push_back(curV);
  }

  return result;
}

// Implement your standard application here
// TO DO: add a parameter for the graph you make.
//
void application(map<long long, Coordinates>& Nodes,
                 vector<FootwayInfo>& Footways, vector<BuildingInfo>& Buildings,
                 graph<long long, double>& G) {
  string person1Building, person2Building;

  cout << endl;
  cout << "Enter person 1's building (partial name or abbreviation), or #> ";
  getline(cin, person1Building);

  while (person1Building != "#") {
    cout << "Enter person 2's building (partial name or abbreviation)> ";
    getline(cin, person2Building);

    // MILESTONE 7. search the buildings
    BuildingInfo person1 = searchBuilding(person1Building, Buildings);
    BuildingInfo person2 = searchBuilding(person2Building, Buildings);

    // edge case to check if the building does not exist
    if (person1.Fullname == "empty") {
      cout << "Person 1's building not found\n";
      cout << "Enter person 1's building (partial name or abbreviation), or #> ";
      getline(cin, person1Building);
      continue;
    }
    // edge case to check if the building is empty
    if (person2.Fullname == "empty") {
      cout << "Person 2's building not found\n";
      cout
          << "Enter person 1's building (partial name or abbreviation), or #> ";
      getline(cin, person1Building);
      continue;
    }

    // Milestone 8
    // Locate Center building
    Coordinates midpoint =
        centerBetween2Points(person1.Coords.Lat, person1.Coords.Lon,
                             person2.Coords.Lat, person2.Coords.Lon);
    BuildingInfo buildingCenter;
    // EDGE CASE
    if (person1.Fullname == person2.Fullname) {
      buildingCenter = person1;
    } else {
      // use the nearestBuilding function to find the center building
      buildingCenter = nearestBuilding(midpoint, Buildings);
    }

    // Milestone 9
    // call the nearest node 3 times
    long long nearest1 = nearestNode(person1, Footways, Nodes);
    long long nearest2 = nearestNode(person2, Footways, Nodes);
    long long nearestCenter = nearestNode(buildingCenter, Footways, Nodes);

    // Milestone 10
    // call Dijkstra's algo twice for ms10 and create 4 maps
    map<long long, long long> Pred_Node1;  // predecessor map
    map<long long, long long> Pred_Node2;  // predecessor map 2

    map<long long, double> distance_map1;  // distance map
    map<long long, double> distance_map2;  // distance map two

    // calling Dijkstra's algorithm here
    Dijkstra(G, nearest1, distance_map1, Pred_Node1);
    Dijkstra(G, nearest2, distance_map2, Pred_Node2);
    //////////////////////////////////////////////////////

    // Print the buildings HERE
    PrintBuilding(person1, person2, buildingCenter);
    // Print the Nodes HERE
    PrintNodes(nearest1, nearest2, nearestCenter, Nodes);

    // this is the SET for last two test cases
    set<string> unvisitedB;  // set for unvisitable buildings.

    // check if the first person's vertex is the same as the second person's
    if (nearest1 == nearest2) {
      // set the predecessor at the center vertex to zero
      distance_map1[nearestCenter] = 0;
      distance_map2[nearestCenter] = 0;
    }
    // ANOTHER EDGE CASE
    if (distance_map1[nearestCenter] >= INF ||
        distance_map2[nearestCenter] >= INF) {
      if (distance_map1[nearest2] >= INF) {
        cout << "Sorry, destination unreachable." << endl;
        cout << "Enter person 1's building (partial name or abbreviation), or "
                "#> ";
        getline(cin, person1Building);
        continue;
      } else {  // for the next closest

        cout << endl;
        cout << "At least one person was unable to reach the destination "
                "building. Finding next closest building..."
             << endl;
        cout << endl;
        unvisitedB.insert(
            buildingCenter
                .Fullname);  // insert into the set if unvisitable buildings

        // loop and add
        while (true) {
          cout << endl;
          // cout << "Please print something" << endl;
          unvisitedB.insert(buildingCenter.Fullname);
          buildingCenter = ModifiedNearest(midpoint, Buildings,
                                           unvisitedB);  // reset the center
          nearestCenter = nearestNode(buildingCenter, Footways,
                                      Nodes);  // set the center to the nearest
          ///////////////////////////
          cout << "New destination building:" << endl;
          cout << " " << buildingCenter.Fullname << endl;
          cout << " (" << buildingCenter.Coords.Lat << ", "
               << buildingCenter.Coords.Lon << ")" << endl;
          cout << "Nearest destination node:" << endl;
          cout << " " << nearestCenter << endl;
          cout << " (" << Nodes[nearestCenter].Lat << ", "
               << Nodes[nearestCenter].Lon << ")" << endl;
          cout << endl;
          // do a check here
          if ((distance_map1[nearestCenter] < INF) ||
              (distance_map2[nearestCenter] < INF)) {
            break;
          }

          cout << "At least one person was unable to reach the destination "
                  "building. Finding next closest building..."
               << endl;
          cout << endl;
        }
      }
    }
    cout << "Person 1's distance to dest: " << distance_map1[nearestCenter]
         << " miles" << endl;
    // cout << endl;
    // Person 1's path
    vector<long long> Path_1 = getPath(Pred_Node1, nearestCenter);
    // loop through and print the vector
    cout << "Path: ";
    for (size_t i = 0; i < Path_1.size() - 1; i++) {
      cout << Path_1[i] << "->";
    }
    cout << Path_1[Path_1.size() - 1] << endl;
    cout << endl;

    cout << "Person 2's distance to dest: " << distance_map2[nearestCenter]
         << " miles" << endl;
    // Person 2's Path
    vector<long long> Path_2 = getPath(Pred_Node2, nearestCenter);
    cout << "Path: ";
    for (size_t j = 0; j < Path_2.size() - 1; j++) {
      cout << Path_2[j] << "->";
    }
    cout << Path_2[Path_2.size() - 1] << endl;

    //
    //
    cout << endl;
    cout << "Enter person 1's building (partial name or abbreviation), or #> ";
    getline(cin, person1Building);
  }
}

int main() {
  // maps a Node ID to it's coordinates (lat, lon)
  map<long long, Coordinates> Nodes;
  // info about each footway, in no particular order
  vector<FootwayInfo> Footways;
  // info about each building, in no particular order
  vector<BuildingInfo> Buildings;
  XMLDocument xmldoc;

  // this is the graph pls work pls
  graph<long long, double> G;

  cout << "** Navigating UIC open street map **" << endl;
  cout << endl;
  cout << std::setprecision(8);

  string def_filename = "map.osm";
  string filename;

  cout << "Enter map filename> ";
  getline(cin, filename);

  if (filename == "") {
    filename = def_filename;
  }

  //
  // Load XML-based map file
  //
  if (!LoadOpenStreetMap(filename, xmldoc)) {
    cout << "**Error: unable to load open street map." << endl;
    cout << endl;
    return 0;
  }

  //
  // Read the nodes, which are the various known positions on the map:
  //
  int nodeCount = ReadMapNodes(xmldoc, Nodes);

  //
  // Read the footways, which are the walking paths:
  //
  int footwayCount = ReadFootways(xmldoc, Footways);

  //
  // Read the university buildings:
  //
  int buildingCount = ReadUniversityBuildings(xmldoc, Nodes, Buildings);

  //
  // Stats
  //
  assert(nodeCount == (int)Nodes.size());
  assert(footwayCount == (int)Footways.size());
  assert(buildingCount == (int)Buildings.size());

  // add the vertices
  // loop through the nodes graph with for each
  for (auto iter : Nodes) {
    G.addVertex(iter.first);  // add the vertex to the graph
  }

  // MILESTONE 6
  // add the edges
  for (auto edges : Footways) {  // loop through the vector
    for (size_t i = 0; i < edges.Nodes.size() - 1; i++) {
      Coordinates n1 = Nodes[edges.Nodes[i]], n2 = Nodes[edges.Nodes[i + 1]];
      // calculate the distance with the distBetween2points function
      double dist = distBetween2Points(n1.Lat, n1.Lon, n2.Lat, n2.Lon);
      // add the edges to the graph now
      G.addEdge(edges.Nodes[i], edges.Nodes[i + 1], dist);
      G.addEdge(edges.Nodes[i + 1], edges.Nodes[i], dist);
    }
  }

  cout << endl;
  cout << "# of nodes: " << Nodes.size() << endl;
  cout << "# of footways: " << Footways.size() << endl;
  cout << "# of buildings: " << Buildings.size() << endl;
  cout << "# of vertices: " << G.NumVertices() << endl;
  cout << "# of edges: " << G.NumEdges() << endl;

  cout << endl;
  string userInput;
  cout << "Enter \"a\" for the standard application or "
       << "\"c\" for the creative component application> ";
  getline(cin, userInput);
  if (userInput == "a") {
    // TO DO: add argument for the graph you make.
    application(Nodes, Footways, Buildings, G);
  } else if (userInput == "c") {
    // TO DO: add arguments
    creative();
  }
  cout << "** Done **" << endl;
  return 0;
}
