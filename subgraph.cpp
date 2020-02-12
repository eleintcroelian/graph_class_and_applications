/**
 * @file subgraph.cpp
 * Test script for viewing a subgraph from our Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <fstream>
#include <iterator>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"

#include "Graph.hpp"

using GraphType = Graph<double,double>;
using NodeType = typename GraphType::node_type;
using NodeIter = typename GraphType::node_iterator;
/** An iterator that skips over elements of another iterator based on whether
 * those elements satisfy a predicate.
 *
 * Given an iterator range [@a first, @a last) and a predicate @a pred,
 * this iterator models a filtered range such that all i with
 * @a first <= i < @a last and @a pred(*i) appear in order of the original range.
 */
template <typename Pred, typename It>
class filter_iterator : private equality_comparable<filter_iterator<Pred, It>>
{
public:
  // Get all of the iterator traits and make them our own
  using value_type = typename std::iterator_traits<It>::value_type;
  using pointer = typename std::iterator_traits<It>::pointer;
  using reference = typename std::iterator_traits<It>::reference;
  using difference_type = typename std::iterator_traits<It>::difference_type;
  using iterator_category = typename std::input_iterator_tag;

  // Constructor
  filter_iterator(const Pred &p, const It &first, const It &last)
      : p_(p), it_(first), end_(last)
  {
    // HW1 #4: YOUR CODE HERE
  }

  // HW1 #4: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  value_type operator*() const
  {
    return (*it_);
  };

  filter_iterator &operator++()
  {
    ++it_;

    while ((it_ != end_) && (!p_(*it_)))
    {
      ++it_;
    }
    return *this;
  };
  bool operator==(const filter_iterator &filtiter) const
  {
    if (*(filtiter.it_) == *it_)
    {
      return true;
    }
    return false;
  };

private:
  Pred p_;
  It it_;
  It end_;
};

/** Helper function for constructing filter_iterators. This deduces the type of
 * the predicate function and the iterator so the user doesn't have to write it.
 * This also allows the use of lambda functions as predicates.
 *
 * Usage:
 * // Construct an iterator that filters odd values out and keeps even values.
 * std::vector<int> a = ...;
 * auto it = make_filtered(a.begin(), a.end(), [](int k) {return k % 2 == 0;});
 */
template <typename Pred, typename Iter>
filter_iterator<Pred, Iter> make_filtered(const Iter &it, const Iter &end,
                                          const Pred &p)
{
  return filter_iterator<Pred, Iter>(p, it, end);
}

// HW1 #4: YOUR CODE HERE
// Specify and write an interesting predicate on the nodes.
// Explain what your predicate is intended to do and test it.
// If you'd like you may create new nodes and tets files.

/** Test predicate for HW1 #4 */
struct SlicePredicate
{
  template <typename NODE>
  bool operator()(const NODE &n)
  {
    return n.position().x < 0;
  }
};
struct Robocop
{
  template <typename NODE>
  bool operator()(const NODE &n)
  {
    return ((n.position().z > t / 2.) || (n.position().z < -t / 2.));
  }
  // Point center_ = Point(0.,0.5,0.);
  // double radius = 0.3;
  double t = 0.1;
};
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

int shortest_path_lengths(GraphType &g, NodeType &root)
{
  unsigned int root_id = root.index();

  std::vector<unsigned int> queue;
  std::vector<bool> visited(g.num_nodes(), false);
  std::vector<double> distances(g.num_nodes(), 0);

  visited[root_id] = true;
  queue.push_back(root_id);
  distances[root_id] = 0;

  while (!queue.empty())
  {
    auto currentnodeid = queue[0];
    auto currentnode = g.node(currentnodeid);

    queue.erase(queue.begin());

    for (auto kt = currentnode.edge_begin(); kt != currentnode.edge_end(); ++kt)
    {
      int adjnodeind;
      if ((*kt).node2().index() == currentnodeid)
      {
        adjnodeind=(*kt).node1().index();
      }
      else
      {
        adjnodeind=(*kt).node2().index();
      }
      if (!visited[adjnodeind])
      {
        visited[adjnodeind] = true;
        queue.push_back(adjnodeind);
        distances[adjnodeind] = distances[currentnodeid] + 1;
      }
    }
  }

  auto maxdist = *(std::max_element(distances.begin(), distances.end()));
  for (auto it = distances.begin(); it != distances.end(); ++it)
  {
    (*it) = (*it) / maxdist;
  }

  for (auto it = g.node_begin(); it != g.node_end(); ++it)
  {
    g.set_value((*it).index(), distances[(*it).index()]);
  }

  return 0;
}
int main(int argc, char **argv)
{
  // Check arguments
  if (argc < 3)
  {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Define our types
  using GraphType = Graph<double,double>;
  using NodeType = typename GraphType::node_type;

  // Construct a Graph
  GraphType graph;
  std::vector<NodeType> nodes;

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

  // HW1 #4: YOUR CODE HERE
  // Use the filter_iterator to plot an induced subgraph.
  Robocop robo;
  auto filtiter = make_filtered(graph.node_begin(), graph.node_end(), robo);
  auto filtiter_end = make_filtered(graph.node_end(), graph.node_end(), robo);

  auto root = nearest_node(graph, Point(0, 5, 0));
  auto rootnode = (*root);
  auto i = shortest_path_lengths(graph, rootnode);
  (void)i;
  colorFunc colorfunctor;

  auto node_map = viewer.empty_node_map(graph);
  viewer.add_nodes(filtiter, filtiter_end, colorfunctor,node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  // Center the view and enter the event loop for interactivity
  viewer.center_view();
  viewer.event_loop();

  return 0;
}
