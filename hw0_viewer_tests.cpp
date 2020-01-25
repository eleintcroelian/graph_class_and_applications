/*
 * Test viewer runtimes
 */

#include <fstream>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"

#include "Graph.hpp"


int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Define our types
  using GraphType = Graph;
  using NodeType  = typename GraphType::node_type;
  using size_type = typename GraphType::size_type;


  // Construct a Graph
  GraphType graph;
  std::vector<NodeType> nodes;
  const clock_t begin_time = clock();
  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  while (CME212::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CME212::getline_parsed(tets_file, t))
    for (unsigned i = 1; i < t.size(); ++i)
      for (unsigned j = 0; j < i; ++j)
        graph.add_edge(nodes[t[i]], nodes[t[j]]);


  // Iterate through nodes via index
  size_type num_nodes = graph.num_nodes();
  for (size_type i = 0; i < num_nodes; ++i) {
    graph.node(i);
  }

  // Iterate through edges via index
  size_type num_edges = graph.num_edges();
  for (size_type i = 0; i < num_edges; ++i) {
    graph.edge(i);
  }
    
  std::cout << "Speed in seconds for " << argv[1] << ": "
    << float( clock () - begin_time ) /  CLOCKS_PER_SEC << std::endl;
  return 0;
}

