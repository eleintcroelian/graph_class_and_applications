/**
 * Unit tests for CME212, hw1
 */

#include <fstream>

// #include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include <algorithm>
#include <vector>
#include <set>
#include <utility>

#include "Graph.hpp"
//#include "poisson.hpp"
//#include "mtl_test.hpp"

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
  using GraphType = Graph<int, int>;
  using NodeType  = typename GraphType::node_type;
  using EdgeType  = typename GraphType::edge_type;
  
  int test_number = atoi(argv[1]);
  bool test_return = false;
  
  // Construct a Graph
  GraphType graph;
  GraphType graphtwo;
  
  // Create vector of 10 Points
  std::vector<Point> points;
  for(int i = 0; i < 10; i++)
    points.push_back(Point(i));
    
  //
  // BASIC TESTING: HW 0 functionality
  //
  if(test_number == 0){
    // Test has_node function
    GraphType::node_type n0 = graph.add_node(points[0]);
    test_return = graph.has_node(n0);
    
  } else if(test_number == 1){
    // Test num nodes/size functions
    test_return = (graph.num_nodes() == graph.size() and
                    graph.size() == 0);
    
    graph.add_node(points[0]);
    graph.add_node(points[1]);
    
    test_return = (test_return and
                    graph.num_nodes() == graph.size() and
                    graph.size() == 2);
  
  } else if(test_number == 2){
    // Test node function
    NodeType n0 = graph.add_node(points[0]);
    test_return = (n0 == graph.node(0));
  
  } else if(test_number == 3){
    // Test index function
    graph.add_node(points[0]);
    graph.add_node(points[1]);
    test_return = (graph.node(0).index() == 0 and
                    graph.node(1).index() == 1);
      
  } else if(test_number == 4){
    // Test position function
    NodeType n0 = graph.add_node(points[0]);
    test_return = (n0.position() == points[0]);
    
  } else if(test_number == 5){
    // Verify trichotomy on nodes
    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph.add_node(points[1]);
    test_return = ((n0 < n1) ^ (n1 < n0));
      
  } else if(test_number == 6){
    // Test adding edge, has_edge function
    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph.add_node(points[1]);
    graph.add_edge(n0, n1);
    test_return = (graph.has_edge(n0, n1) and
                    graph.has_edge(n1, n0));
    
  } else if(test_number == 7){
    // Test num_edges function
    test_return = (graph.num_edges() == 0);
    
    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph.add_node(points[1]);
    NodeType n2 = graph.add_node(points[2]);
    
    graph.add_edge(n0, n1);
    graph.add_edge(n1, n2);
    
    test_return = test_return and (graph.num_edges() == 2);
    
  } else if(test_number == 8){
    // Test edge function
    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph.add_node(points[1]);
    graph.add_node(points[2]);    

    EdgeType e0 = graph.add_edge(n0, n1);
    test_return = (e0 == graph.edge(0));
    
  } else if(test_number == 9){
    // Test node1, node2 functions
    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph.add_node(points[1]);
    EdgeType e0 = graph.add_edge(n0, n1);
    
    test_return = ((e0.node1() == n0 and e0.node2() == n1) or
                    (e0.node1() == n1 and e0.node2() == n0));
      
  } else if(test_number == 10){
    // Verify only one of e0 < e1 or e1 < e0 is true
    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph.add_node(points[1]);
    NodeType n2 = graph.add_node(points[2]);
    
    EdgeType e0 = graph.add_edge(n0, n1);
    EdgeType e1 = graph.add_edge(n1, n2);
    
    test_return = ((e0 < e1) ^ (e1 < e0));
  } else if(test_number == 11){
    // Graph properly cleared
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
    NodeType n0 = graph.add_node(points[0]);
    test_return = (sizeof(n0) <= 16);

  } else if(test_number == 13){
    // Test Edge size <= 32
    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph.add_node(points[1]);
    
    EdgeType e0 = graph.add_edge(n0, n1);
    test_return = (sizeof(e0) <= 32);

  } else if(test_number == 14){
    // Add existing edge doesn't create another edge
    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph.add_node(points[1]);

    graph.add_edge(n0, n1);
    graph.add_edge(n0, n1);
    graph.add_edge(n1, n0);
    test_return = (graph.num_edges() == 1);

  } else if(test_number == 15){
    // Nodes of different graphs aren't equal
    GraphType graph2;
    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph2.add_node(points[0]);
    test_return = !(n0 == n1);
    
  } else if(test_number == 16){
    // Nodes at same position aren't equal
    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph.add_node(points[0]);
    test_return = !(n0 == n1);
    
  } else if(test_number == 17){
    // Invalid Nodes don't break
    NodeType x;
    x = graph.add_node(points[0]);
    test_return = (x.index() == 0);
    
  } else if(test_number == 18){
    // add_edge returns correct existing edge
    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph.add_node(points[1]);
    NodeType n2 = graph.add_node(points[2]);
    
    EdgeType e0 = graph.add_edge(n0, n1);
    EdgeType e1 = graph.add_edge(n1, n2);
    test_return = (graph.add_edge(n0, n1) == e0 and
                   graph.add_edge(n1, n0) == e0 and
                   graph.add_edge(n1, n2) == e1 and
                   graph.add_edge(n2, n1) == e1);
    
  } else if(test_number == 19){
    // Clearing graph still allows interaction 1
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
    
  } else if(test_number == 20){
    // Clearing graph still allows interaction 2
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
  else if(test_number == 21){
    // Test adding node with default value
    graph.add_node(points[0]);
    test_return = graph.node(0).value() == 0;
  } else if(test_number == 22){
    // Test adding nodes with provided values
    graph.add_node(points[0], 5);
    graph.add_node(points[1], 6);
    test_return = (graph.node(0).value() == 5 and graph.node(1).value() == 6);
    test_return = test_return or (graph.node(0).value() == 6 and graph.node(1).value() == 5);
  } else if(test_number == 23){
    // Test node iterator
    graph.add_node(points[0]);
    graph.add_node(points[1]);
    graph.add_node(points[2]);
    
    int iter = 0;
    for(auto ni = graph.node_begin(); ni != graph.node_end(); ++ni){
      ++iter;
    }
    test_return = (iter == 3);
  
  } else if(test_number == 24){
    // Test degree function
    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph.add_node(points[1]);
    NodeType n2 = graph.add_node(points[2]);
    graph.add_edge(n0, n1);
    test_return = (n2.degree() == 0 and n1.degree() == 1);

  } else if(test_number == 25){
    // Test incident iterator: node1 accurate?
    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph.add_node(points[1]);
    
    graph.add_edge(n0, n1);

    test_return = ((*(n1.edge_begin())).node1()==n1);
    
  } else if(test_number == 26){
    // Test incident iterator

    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph.add_node(points[1]);
    NodeType n2 = graph.add_node(points[2]);
    
    graph.add_edge(n0, n1);
    graph.add_edge(n1, n2);
    graph.add_edge(n0, n2);

    int iter = 0;
    for(auto ni = n1.edge_begin(); ni != n1.edge_end(); ++ni){
      ++iter;
    }
    test_return = (iter == 2);
      
  } else if(test_number == 27){
    // Test edge iterator
    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph.add_node(points[1]);
    NodeType n2 = graph.add_node(points[2]);
    graph.add_node(points[3]);
    
    graph.add_edge(n0, n1);
    graph.add_edge(n1, n2);
    graph.add_edge(n0, n2);

    // This last edge should be the same as e2
    graph.add_edge(n2, n0);

    int iter = 0;
    for(auto ni = graph.edge_begin(); ni != graph.edge_end(); ++ni){
      ++iter;
    }
    test_return = (iter == 3);
    
  } 
  
  // BEGIN HW2 Functionality
  else if (test_number == 28) {
    // Test modifiable node position
    NodeType n0 = graph.add_node(points[0]);
    n0.position() = points[1];
    test_return = (n0.position() == points[1]);

  } else if (test_number == 29) {
    // Test default edge value
    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph.add_node(points[1]);
    EdgeType e0 = graph.add_edge(n0, n1);
    test_return = (e0.value() == 0);

  } else if (test_number == 30) {
    // Test modifiable edge value
    NodeType n0 = graph.add_node(points[0]);
    NodeType n1 = graph.add_node(points[1]);
    EdgeType e0 = graph.add_edge(n0, n1);
    e0.value() = 1;
    test_return = (e0.value() == 1);

  } else if (test_number == 31) {
    // First node remove test from test_nodes
    NodeType n0 = graph.add_node(points[0]);
    graph.remove_node(n0);
    test_return = (graph.num_nodes() == 0);

  } else if (test_number == 32) {
    // Second node remove test from test_nodes
    for (int k = 0; k < 100; ++k) {
      graph.add_node(
        Point(CME212::random(), CME212::random(), CME212::random()), k);
    }
    for (unsigned k = 0; k < 50; ++k) {
      unsigned n = unsigned(CME212::random(0, graph.num_nodes()));
      graph.remove_node(graph.node(n));
    }
    test_return = (graph.num_nodes() == 50);

  } else if (test_number == 33) {
    // Third node remove test from test_nodes
    // Tests node invariant
    for (int k = 0; k < 100; ++k) {
      graph.add_node(
        Point(CME212::random(), CME212::random(), CME212::random()), k);
    }
    for (unsigned k = 0; k < 50; ++k) {
      unsigned n = unsigned(CME212::random(0, graph.num_nodes()));
      graph.remove_node(graph.node(n));
    }
    test_return = true;
    for (unsigned k = 0; test_return && k < 50; ++k) {
      NodeType node = graph.node(k);
      if (node.index() != k) {
        test_return = false;
      }
    }

  } else if (test_number == 34) {
    // Fourth node remove test
    // Tests node inequality between graphs after remove
    for (unsigned k = 0; k < 50; ++k) {
      Point p(CME212::random(), CME212::random(), CME212::random());
      graph.add_node(p);
      graphtwo.add_node(p);
    }
    graph.remove_node(graph.node(25));
    graphtwo.remove_node(graphtwo.node(25));
    test_return = (
      graphtwo.node(48) == graphtwo.node(48)
      && graphtwo.node(29) != graphtwo.node(48)
      && graphtwo.node(29) != graph.node(29)
      && graphtwo.node(48) != graph.node(29));

  } else if (test_number == 35) {
    // Test position after node removal
    for (unsigned k = 0; k < 50; ++k) {
      Point p(CME212::random(), CME212::random(), CME212::random());
      graph.add_node(p);
    }
    Point q = graph.node(25).position();
    graph.remove_node(graph.node(25));
    test_return = (graph.node(25).position() != q);

  } else if (test_number == 36) {
    // Test node removal with iterators
    for (int k = 0; k < 10; ++k) {
      graph.add_node(
        Point(CME212::random(), CME212::random(), CME212::random()), k);
    }
    auto it = graph.node_begin();
    for (int i = 0; i < 5; ++i) {
      ++it;
    }
    graph.remove_node(it);
    test_return = (graph.num_nodes() == 9);

  } else if (test_number == 37) {
    // Tougher node removal test with iterators
    for (int k = 0; k < 10; ++k) {
      graph.add_node(
        Point(CME212::random(), CME212::random(), CME212::random()), k);
    }
    auto it = graph.node_begin();
    for (int i = 0; i < 5; ++i) {
      it = graph.remove_node(it);
    }
    test_return = (graph.num_nodes() == 5);

  } else if (test_number == 38) {
    // Toughest node removal test with iterators
    std::vector<int> nodevec;
    for (int k = 0; k < 10; ++k) {
      NodeType n = graph.add_node(
        Point(CME212::random(), CME212::random(), CME212::random()), k);
      nodevec.push_back(n.value());
    }
    auto it = graph.node_begin();
    for (int i = 0; i < 5; ++i) {
      nodevec.erase(std::find(nodevec.begin(), nodevec.end(), (*it).value()));
      it = graph.remove_node(it);
    }

    test_return = true;
    for (unsigned int i = 0; i < 5; ++i) {
      auto it = graph.node_begin();
      for (; (*it).value() != nodevec[i]; ++it);
      if (graph.node(i).index() != i || graph.node((*it).index()).value() != (*it).value()) {
        test_return = false;
      }
    }

  } else if (test_number == 39) {
    // First remove edge test from test_edges
    Point p1(CME212::random(), CME212::random(), CME212::random());
    Point p2(CME212::random(), CME212::random(), CME212::random());

    graph.add_node(p1);
    graph.add_node(p2);
    graph.add_edge(graph.node(0), graph.node(1));
    graph.remove_edge(graph.node(0), graph.node(1));
    test_return = !graph.has_edge(graph.node(1), graph.node(0));

  } else if (test_number == 40) {
    // Second remove edge test from test_edges
    Point p1(CME212::random(), CME212::random(), CME212::random());
    Point p2(CME212::random(), CME212::random(), CME212::random());

    graph.add_node(p1);
    graph.add_node(p2);
    graph.add_edge(graph.node(0), graph.node(1));
    graph.remove_edge(graph.node(0), graph.node(1));
    test_return = !graph.remove_edge(graph.node(1), graph.node(0));

  } else if (test_number == 41) {
    // Third remove edge test from test_edges
    for (int k = 0; k < 100; ++k) {
      graph.add_node(
        Point(CME212::random(), CME212::random(), CME212::random()));
    }

    // Add edges
    for (unsigned k = 0; k < 100; ++k) {
      unsigned n1, n2;
      do {
        n1 = (unsigned) CME212::random(0, graph.num_nodes());
        n2 = (unsigned) CME212::random(0, graph.num_nodes());
      } while (n1 == n2 || graph.has_edge(graph.node(n1), graph.node(n2)));
      graph.add_edge(graph.node(n1), graph.node(n2));
    }

    // Remove edges
    for (unsigned k = 0; k < 20; ++k) {
      unsigned n1, n2;
      do {
        n1 = (unsigned) CME212::random(0, graph.num_nodes());
        n2 = (unsigned) CME212::random(0, graph.num_nodes());
      } while (!graph.has_edge(graph.node(n1), graph.node(n2)));
      graph.remove_edge(graph.node(n1), graph.node(n2));
    }

    unsigned count_edges = 0;
    for (unsigned k = 0; k < graph.num_nodes(); ++k) {
      for (unsigned j = k+1; j < graph.num_nodes(); ++j) {
        if (graph.has_edge(graph.node(k), graph.node(j))) {
          ++count_edges;
        }
      }
    }

    test_return = (graph.num_edges() == count_edges);

  } else if (test_number == 42) {
    // Test removal after changing node position
    for (int k = 0; k < 100; ++k) {
      NodeType n = graph.add_node(
        Point(CME212::random(), CME212::random(), CME212::random()));
      n.value() = k;
    }

    // Add edges
    std::set<std::pair<int, int>> edgeset;
    for (unsigned k = 0; k < 200; ++k) {
      unsigned n1, n2;
      do {
        n1 = (unsigned) CME212::random(0, graph.num_nodes());
        n2 = (unsigned) CME212::random(0, graph.num_nodes());
      } while (n1 == n2 || graph.has_edge(graph.node(n1), graph.node(n2)));
      NodeType node1 = graph.node(n1);
      NodeType node2 = graph.node(n2);
      graph.add_edge(node1, node2);
      edgeset.insert(std::make_pair(node1.value(), node2.value()));
    }

    // Move nodes
    for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
      (*it).position() += 0.1 * Point(
        CME212::random(), CME212::random(), CME212::random());
    }

    // Remove nodes
    for (unsigned k = 0; k < 10; ++k) {
      unsigned n = (unsigned) CME212::random(0, graph.num_nodes());
      graph.remove_node(graph.node(n));
    }

    // Move nodes
    for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
      (*it).position() += 0.1 * Point(
        CME212::random(), CME212::random(), CME212::random());
    }

    test_return = true;
    for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it) {
      if (
	  edgeset.find(std::make_pair((*it).node1().value(), (*it).node2().value())) == edgeset.end()
	  && edgeset.find(std::make_pair((*it).node2().value(), (*it).node1().value())) == edgeset.end()) {
        test_return = false;
      }
    }

  }
  return (!test_return);
}

