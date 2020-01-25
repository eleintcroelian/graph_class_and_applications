/**
 * Unit tests for CME212, hw0
 */

#include <fstream>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"

#include "Graph.hpp"

void print_desc(const char *desc){
  std::cout << "Test: " << desc << std::endl;
}

// Argument is the test number
int main(int argc, char** argv)
{
  // Check arguments
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " test-number\n";
    exit(1);
  }

  // Define our types
  using GraphType = Graph;
  using NodeType  = typename GraphType::node_type;
  using EdgeType  = typename GraphType::edge_type;
  
  int test_number = atoi(argv[1]);
  bool test_return = false;
  
  // Construct a Graph
  GraphType graph;
  
  // Create vector of 10 Points
  std::vector<Point> points;
  for(int i = 0; i < 10; i++)
    points.push_back(Point(i));
    
  //
  // BASIC TESTING
  //
  if(test_number == 0){
    // Test has_node function
    print_desc("Test has_node function");
    GraphType::node_type n0 = graph.add_node(points[0]);
    test_return = graph.has_node(n0);
    
  } else if(test_number == 1){
    // Test num nodes/size functions
    print_desc("Test num nodes/size functions");
    test_return = (graph.num_nodes() == graph.size() and
                    graph.size() == 0);
    
    graph.add_node(points[0]);
    graph.add_node(points[1]);
    
    test_return = (test_return and
                    graph.num_nodes() == graph.size() and
                    graph.size() == 2);
  
  } else if(test_number == 2){
    // Test node function
    print_desc("Test node function");
    NodeType n0 = graph.add_node(points[0]);
    test_return = (n0 == graph.node(0));
  
  } else if(test_number == 3){
    // Test index function
    print_desc("Test index function");
    graph.add_node(points[0]);
    graph.add_node(points[1]);
    test_return = (graph.node(0).index() == 0 and
                    graph.node(1).index() == 1);
      
  } else if(test_number == 4){
    // Test position function
    print_desc("Test position function");
    NodeType n0 = graph.add_node(points[0]);
    test_return = (n0.position() == points[0]);
    
  } else if(test_number == 5){
    // Verify trichotomy on nodes
    print_desc("Verify trichotomy on nodes");
    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph.add_node(points[1]);
    test_return = ((n0 < n1) ^ (n1 < n0));
      
  } else if(test_number == 6){
    // Test adding edge, has_edge function
    print_desc("Test add_edge, has_edge functions");
    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph.add_node(points[1]);
    graph.add_edge(n0, n1);
    test_return = (graph.has_edge(n0, n1) and
                    graph.has_edge(n1, n0));
    
  } else if(test_number == 7){
    // Test num_edges function
    print_desc("Test num_edges function");
    test_return = (graph.num_edges() == 0);
    
    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph.add_node(points[1]);
    NodeType n2 = graph.add_node(points[2]);
    
    graph.add_edge(n0, n1);
    graph.add_edge(n1, n2);
    
    test_return = test_return and (graph.num_edges() == 2);
    
  } else if(test_number == 8){
    // Test edge function
    print_desc("Test edge function");
    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph.add_node(points[1]);
    graph.add_node(points[2]);    

    EdgeType e0 = graph.add_edge(n0, n1);
    test_return = (e0 == graph.edge(0));
    
  } else if(test_number == 9){
    // Test node1, node2 functions
    print_desc("Test node1, node2 functions");
    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph.add_node(points[1]);
    EdgeType e0 = graph.add_edge(n0, n1);
    
    test_return = ((e0.node1() == n0 and e0.node2() == n1) or
                    (e0.node1() == n1 and e0.node2() == n0));
      
  } else if(test_number == 10){
    // Verify only one of e0 < e1 or e1 < e0 is true
    print_desc("Verify only one of e0 < e1 or e1 < e0 is true");
    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph.add_node(points[1]);
    NodeType n2 = graph.add_node(points[2]);
    
    EdgeType e0 = graph.add_edge(n0, n1);
    EdgeType e1 = graph.add_edge(n1, n2);
    
    test_return = ((e0 < e1) ^ (e1 < e0));
  } else if(test_number == 11){
    // Graph properly cleared
    print_desc("Graph properly cleared");
    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph.add_node(points[1]);
    NodeType n2 = graph.add_node(points[2]);
    
    graph.add_edge(n0, n1);
    graph.add_edge(n1, n2);
    
    graph.clear();
    
    test_return = (graph.num_nodes() == 0 and
                    graph.num_edges() == 0);
  } else if(test_number == 12){
    // Test Node size <= 16
    print_desc("Node size <= 16");
    NodeType n0 = graph.add_node(points[0]);
    test_return = (sizeof(n0) <= 16);

  } else if(test_number == 13){
    // Test Edge size <= 32
    print_desc("Edge size <= 32");
    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph.add_node(points[1]);
    
    EdgeType e0 = graph.add_edge(n0, n1);
    test_return = (sizeof(e0) <= 32);

  } else if(test_number == 14){
    // Add existing edge doesn't create another edge
    print_desc("Add existing edge doesn't create another edge");
    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph.add_node(points[1]);

    graph.add_edge(n0, n1);
    graph.add_edge(n0, n1);
    graph.add_edge(n1, n0);
    test_return = (graph.num_edges() == 1);

  }

  //
  // MEDIUM-LEVEL TEST CASES
  //
  int medium = 15;
  if(test_number == medium){
    // Nodes of different graphs aren't equal
    print_desc("Nodes of different graphs aren't equal");
    Graph graph2;
    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph2.add_node(points[0]);
    test_return = !(n0 == n1);
    
  } else if(test_number == medium+1){
    // Nodes at same position aren't equal
    print_desc("Nodes at same position aren't equal");
    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph.add_node(points[0]);
    test_return = !(n0 == n1);
    
  } else if(test_number == medium+2){
    // Invalid Nodes don't break
    print_desc("Invalid Nodes don't break");
    NodeType x;
    x = graph.add_node(points[0]);
    test_return = (x.index() == 0);
    
  } else if(test_number == medium+3){
    // add_edge returns correct existing edge
    print_desc("add_edge returns correct existing edge");
    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph.add_node(points[1]);
    NodeType n2 = graph.add_node(points[2]);
    
    EdgeType e0 = graph.add_edge(n0, n1);
    EdgeType e1 = graph.add_edge(n1, n2);
    test_return = (graph.add_edge(n0, n1) == e0 and
                   graph.add_edge(n1, n0) == e0 and
                   graph.add_edge(n1, n2) == e1 and
                   graph.add_edge(n2, n1) == e1);
    
  } else if(test_number == medium+4){
    // Clearing graph still allows interaction 1
    print_desc("Clearing graph still allows interaction 1");
    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph.add_node(points[1]);
    NodeType n2 = graph.add_node(points[2]);
    
    EdgeType e0 = graph.add_edge(n0, n1);
    EdgeType e1 = graph.add_edge(n1, n2);
    
    graph.clear();
    
    n0 = graph.add_node(points[0]);
    n1 = graph.add_node(points[1]);
    n2 = graph.add_node(points[2]);
    
    e0 = graph.add_edge(n0, n1);
    e1 = graph.add_edge(n1, n2);
    
    test_return = (graph.num_nodes() == 3 and
                    graph.num_edges() == 2);
    
  } else if(test_number == medium+5){
    // Clearing graph still allows interaction 2
    print_desc("Clearing graph still allows interaction 2");
    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph.add_node(points[1]);
    NodeType n2 = graph.add_node(points[2]);
    
    EdgeType e0 = graph.add_edge(n0, n1);
    EdgeType e1 = graph.add_edge(n1, n2);
    
    graph.clear();
    
    n0 = graph.add_node(points[0]);
    n1 = graph.add_node(points[1]);
    n2 = graph.add_node(points[2]);
    
    e0 = graph.add_edge(n0, n1);
    e1 = graph.add_edge(n1, n2);
    
    test_return = (graph.has_node(n0) and
                    graph.has_edge(n0, n1));
  }
  
  return (!test_return);
}

