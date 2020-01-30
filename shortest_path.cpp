/**
 * @file shortest_path.cpp
 * Test script for using our templated Graph to determine shortest paths.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <vector>
#include <fstream>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"

#include "Graph.hpp"

// Define our types
using GraphType = Graph<double>;
using NodeType = typename GraphType::node_type;
using NodeIter = typename GraphType::node_iterator;

/** Find the node with the minimum euclidean distance to a point.
 * @param g  The graph of nodes to search.
 * @param point  The point to use as the query.
 * @return An iterator to the node of @a g with the minimun Eucliean
 *           distance to @a point.
 *           graph.node_end() if graph.num_nodes() == 0.
 *
 * @post For all i, 0 <= i < graph.num_nodes(),
 *          norm(point - *result) <= norm(point - g.node(i).position())
 */

struct distanceFunc
{
  distanceFunc(const Point Refpoint) : _refpoint(Refpoint) {}

  bool operator()(GraphType::node_type n1, GraphType::node_type n2)
  {
    auto point1 = n1.position();
    auto point2 = n2.position();
    auto temppoint1 = point1 - _refpoint;
    auto dist1 = norm(temppoint1);
    auto temppoint2 = point2 - _refpoint;
    auto dist2 = norm(temppoint2);

    return (dist1 > dist2);
  }
  const Point _refpoint;
};

struct colorFunc
{
  colorFunc() {}

  CME212::Color operator()(GraphType::node_type n1)
  {
    return (CME212::Color::make_heat(n1.value()));
  }
};

NodeIter nearest_node(const GraphType &g, const Point &point)
{
  // HW1 #3: YOUR CODE HERE
  distanceFunc distance(point);
  auto result = std::min_element(g.node_begin(), g.node_end(), distance);
  return result;
}

/** Update a graph with the shortest path lengths from a root node.
 * @param[in,out] g     Input graph
 * @param[in,out] root  Root node to start the search.
 * @return The maximum path length found.
 *
 * @post root.value() == 0
 * @post Graph has modified node values indicating the minimum path length
 *           to the root.
 * @post Graph nodes that are unreachable from the root have value() == -1.
 *
 * This sets all nodes' value() to the length of the shortest path to
 * the root node. The root's value() is 0. Nodes unreachable from
 * the root have value() -1.
 */

unsigned int minDistance(std::vector<double> &dist, std::vector<bool> &sptSet, GraphType &g)
{
  // Initialize min value
  double min = std::numeric_limits<double>::max();
  unsigned int min_index;

  for (auto it = g.node_begin(); it != g.node_end(); ++it)
    if (sptSet[(*it).index()] == false && dist[(*it).index()] <= min)
      min = dist[(*it).index()], min_index = (*it).index();
  return min_index;
}

double distancebetween(GraphType::Node n1, GraphType::Node n2)
{
  auto p1 = n1.position();
  auto p2 = n2.position();
  auto temppoint = p2 - p1;
  return norm(temppoint);
}

int shortest_path_lengths(GraphType &g, NodeType &root)
{
  //Dijkstra's Algorithm

  unsigned int root_id = root.index();
  std::vector<double> distances(g.num_nodes(), std::numeric_limits<double>::max()); // Distance vector initialized to inf except for root
  std::vector<bool> sptset(g.num_nodes(), false);
  distances[root_id] = 0.;

  for (auto it = g.node_begin(); it != g.node_end(); ++it)
  {
    int currentnodeindex = minDistance(distances, sptset, g);
    auto currentnode = g.node(currentnodeindex);
    sptset[currentnodeindex] = true;

    // std::cout << "Current Node " << currentnodeindex << std::endl;
    auto inc_iter_end = currentnode.edge_end();
    auto inc_iter = currentnode.edge_begin();

    for (; inc_iter != inc_iter_end; ++inc_iter)
    {
      auto adjnode = (*inc_iter).node2();
      double dist2currentnode = distancebetween(currentnode, adjnode);
      // std::cout << "Adjacent Node " << adjnode.index() << "Distance from current node " << dist2currentnode << std::endl;
    }

    // Iterate over the adjacent nodes and update their dist values only if they are not in sptset
    // and total weight of path from root is smaller than current value of distances[current node]
    for (auto inc_iter = currentnode.edge_begin(); inc_iter != currentnode.edge_end(); ++inc_iter)
    {
      auto adjnode = (*inc_iter).node2();
      double dist2currentnode = distancebetween(currentnode, adjnode);

      // std::cout << "Adjacent Node " << adjnode.index() << "Distance from current node " << dist2currentnode << std::endl;

      if ((sptset[adjnode.index()] == false) && (distances[currentnodeindex] < dist2currentnode + distances[adjnode.index()]))
      {
        distances[adjnode.index()] = distances[currentnode.index()] + dist2currentnode;
      }
    }
  }

  for (auto it = distances.begin(); it != distances.end(); ++it)
  {
    if ((*it) > 1.e200)
    {
      (*it) = 0.;
    }
  }

  auto maxdist = *(std::max_element(distances.begin(), distances.end()));
      // std::cout <<"Max Dist is ---> " << maxdist<< std::endl;

  for (auto it = distances.begin(); it != distances.end(); ++it)
  {
    (*it) = (*it) / maxdist;
    // std::cout << (*it) << std::endl;
  }

  for (auto it = g.node_begin(); it != g.node_end(); ++it)
  {
    g.set_value((*it).index(), distances[(*it).index()]);
    // std::cout <<distances[(*it).index()] << std::endl;
  }
  return 0;
}

void colorize(GraphType &g)
{
  for (auto it = g.node_begin(); it != g.node_end(); ++it)
  {
  }
}

int main(int argc, char **argv)
{
  // Check arguments
  if (argc < 3)
  {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a Graph
  GraphType graph;
  std::vector<GraphType::node_type> nodes;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  while (CME212::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int, 4> t;
  while (CME212::getline_parsed(tets_file, t))
    for (unsigned i = 1; i < t.size(); ++i)
      for (unsigned j = 0; j < i; ++j)
        graph.add_edge(nodes[t[i]], nodes[t[j]]);

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SFML_Viewer
  CME212::SFML_Viewer viewer;

  // HW1 #3: YOUR CODE HERE
  // Use nearest_node and shortest_path_lengths to set the node values
  // Construct a Color functor and view with the SFML_Viewer

  auto root = nearest_node(graph, Point(-1, 0, 1));
  auto rootnode = (*root);
  auto i = shortest_path_lengths(graph, rootnode);

  auto node_map = viewer.empty_node_map(graph);
  colorFunc colorfunctor;

  viewer.add_nodes(graph.node_begin(), graph.node_end(), colorfunctor, node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  // Center the view and enter the event loop for interactivity
  viewer.center_view();
  viewer.event_loop();

  return 0;
}
