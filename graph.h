// graph.h <Starter Code>
//<Abeer Fatima>
//
// Basic graph class using adjacency matrix representation.  Currently
// limited to a graph with at most 100 vertices.
//
// University of Illinois at Chicago
// CS 251: Fall 2021
// Project #7 - Openstreet Maps
//

#pragma once

#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

using namespace std;

template <typename VertexT, typename WeightT>
class graph {
 private:
  map<VertexT, map<VertexT, WeightT>> AdjList;
  int NumberEdges = 0;

 public:
  //
  // constructor:
  //
  // Constructs an empty graph where n is the max # of vertices
  // you expect the graph to contain.
  //
  // NOTE: the graph is implemented using an adjacency matrix.
  // If n exceeds the dimensions of this matrix, an exception
  // will be thrown to let you know that this implementation
  // will not suffice.
  //
  graph() {}

  //
  // NumVertices
  //
  // Returns the # of vertices currently in the graph.
  int NumVertices() const {
    return (int)AdjList.size();  // the size is just returned
  }

  //
  // NumEdges
  //
  // Returns the # of edges currently in the graph.
  //
  int NumEdges() const {
    // iterate thru the map here
    return NumberEdges;
  }

  //
  // addVertex
  //
  // Adds the vertex v to the graph if there's room, and if so
  // returns true.  If the graph is full, or the vertex already
  // exists in the graph, then false is returned.
  //
  bool addVertex(VertexT v) {
    if (AdjList.count(v) == 0) {  // check if vertex even exsists
      AdjList.emplace(
          v,
          map<VertexT, WeightT>());
      return true;
    }
    return false;
  }

  //
  // addEdge
  //
  // Adds the edge (from, to, weight) to the graph, and returns
  // true.  If the vertices do not exist or for some reason the
  // graph is full, false is returned.
  //
  // NOTE: if the edge already exists, the existing edge weight
  // is overwritten with the new edge weight.
  //
  bool addEdge(VertexT from, VertexT to, WeightT weight) {
    if (AdjList.count(to) == 0 || AdjList.count(from) == 0) {
      return false;  // vertex does not exist, so return false
    }

    if (AdjList.at(from).count(to) != 0) {
      AdjList.at(from).at(to) = weight;  // update the weight at this "index"
    } else {

      AdjList.at(from).emplace(to, weight);
      NumberEdges += 1;  // increase the number of edges
    }

    return true;  // return true after updating the weight
  }

  //
  // getWeight
  //
  // Returns the weight associated with a given edge.  If
  // the edge exists, the weight is returned via the reference
  // parameter and true is returned.  If the edge does not
  // exist, the weight parameter is unchanged and false is
  // returned.
  //
  bool getWeight(VertexT from, VertexT to, WeightT& weight) const {
    if (AdjList.count(from) == 0) {
      return false;
    } else if (AdjList.count(to) == 0) {
      return false;
    } else if (AdjList.at(from).count(to) == 0) {
      return false;
    } else {
      weight = AdjList.at(from).at(to);
    }
    return true;
  }

  // Returns a set containing the neighbors of v, i.e. all
  // vertices that can be reached from v along one edge.
  // Since a set is returned, the neighbors are returned in
  // sorted order; use foreach to iterate through the set.
  //
  set<VertexT> neighbors(VertexT v) const {
    set<VertexT> S;
    if (AdjList.count(v) == 0) {
      return S;  // return the empty set if the vertex does not exist
    }
    for (auto iter : AdjList.at(v)) {
      S.emplace(iter.first);
    }
    return S;
  }

  //
  // getVertices
  //
  // Returns a vector containing all the vertices currently in
  // the graph.
  //
  vector<VertexT> getVertices() const {
    vector<VertexT> vertV;
    for (auto iter : AdjList) {
      vertV.push_back(iter.first);  // append the vertices to the vector
    }

    return vertV;
  }

  //
  // dump
  //
  // Dumps the internal state of the graph for debugging purposes.
  //
  // Example:
  //    graph<string,int>  G(26);
  //    ...
  //    G.dump(cout);  // dump to console
  //
  void dump(ostream& output) const {
         output << "***************************************************" <<
         endl; output << "********************* GRAPH ***********************"
         << endl;

         output << "**Num vertices: " << this->NumVertices() << endl;
         output << "**Num edges: " << this->NumEdges() << endl;

         output << endl;
         output << "**Vertices:" << endl;
         for (int i = 0; i < this->NumVertices(); ++i) {
           output << " " << i << ". " << this->Vertices[i] << endl;
         }

         output << endl;
         output << "**Edges:" << endl;
         for (int row = 0; row < this->NumVertices(); ++row) {
           output << " row " << row << ": ";

           for (int col = 0; col < this->NumVertices(); ++col) {
             if (this->AdjMatrix[row][col].EdgeExists == false) {
               output << "F ";
             } else {
               output << "(T,"
                 << this->AdjMatrix[row][col].Weight
                 << ") ";
             }
           }
           output << endl;
         }
         output << "**************************************************" <<
         endl;
  }
};
