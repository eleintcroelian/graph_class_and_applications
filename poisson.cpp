/**
 * @file poisson.cpp
 * Test script for treating the Graph as a MTL Matrix
 * and solving a Poisson equation.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles.
 * Second file: Eges (one per line) defined by 2 indices into the point list
 *              of the first file.
 *
 * Launches an SFML_Viewer to visualize the solution.
 */
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include <fstream>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"
#include "CME212/BoundingBox.hpp"

#include "Graph.hpp"
#include <math.h>
#include <algorithm>
using GraphType = Graph<char, char>; //<  DUMMY Placeholder
using NodeType = typename GraphType::node_type;
bool is_boundary(NodeType n)
{
  Box3D B(Point(-0.6, -0.2, -1), Point(0.6, 0.2, 1));
  auto p = n.position();
  if ((norm_inf(p) == 1) ||
      (norm_inf(p - Point(0.6, 0.6, 0)) < 0.2) ||
      (norm_inf(p - Point(-0.6, 0.6, 0)) < 0.2) ||
      (norm_inf(p - Point(0.6, -0.6, 0)) < 0.2) ||
      (norm_inf(p - Point(-0.6, -0.6, 0)) < 0.2) ||
      B.contains(p))
  {
    return true;
  }
  return false;
};

class g_func
{
public:
  g_func(){};
  double operator()(Point p) const
  {
    if (norm_inf(p) == 1)
    {
      return 0;
    }
    if ((norm_inf(p - Point(0.6, 0.6, 0)) < 0.2) ||
        (norm_inf(p - Point(-0.6, 0.6, 0)) < 0.2) ||
        (norm_inf(p - Point(0.6, -0.6, 0)) < 0.2) ||
        (norm_inf(p - Point(-0.6, -0.6, 0)) < 0.2))
    {
      return -0.2;
    }
    if (B.contains(p))
    {
      return 1;
    }
  };
  Box3D B = Box3D(Point(-0.6, -0.2, -1), Point(0.6, 0.2, 1));
};

class GraphSymmetricMatrix
{
  /** Helper function to perform multiplication . Allows for delayed
* evaluation of results .
* Assign :: apply (a, b) resolves to an assignment operation such as
* a += b, a -= b, or a = b.
* @pre @a size (v) == size (w) */
public:
  GraphSymmetricMatrix(const Graph<char, char> &g, std::vector<unsigned int> &boundary) : g_(&g), bound_(boundary){};

  template <typename VectorIn, typename VectorOut, typename Assign>
  void mult(const VectorIn &v, VectorOut &w, Assign) const
  {
    for (size_t i = 0; i < g_->num_nodes(); i++) //row
    {
      auto n_i = g_->node(i);
      double temp = 0.;
      if (is_boundary(n_i))
      {
        temp += v[i];
      }
      else
      {
        double deg = n_i.degree();
        temp -= deg * v[i];
        for (auto it = n_i.edge_begin(); it != n_i.edge_end(); ++it)
        {
          auto node2 = (*it).node2();
          if (!is_boundary(node2))
          {
            temp += v[node2.index()];
          }
        }
      }
      Assign::apply(w[i], temp);
      // Assign::apply(w[i], temp);
    }
  }
  /** Matvec forwards to MTL 's lazy mat_cvec_multiplier operator */
  template <typename Vector>
  mtl::vec::mat_cvec_multiplier<GraphSymmetricMatrix, Vector>
  operator*(const Vector &v) const
  {
    return mtl::vec::mat_cvec_multiplier<GraphSymmetricMatrix, Vector>(*this, v);
  }
  const Graph<char, char> *g_;
  std::vector<unsigned int> bound_;
};

inline std::size_t size(const GraphSymmetricMatrix &A) { return A.g_->num_nodes() * A.g_->num_nodes(); }
inline std::size_t num_rows(const GraphSymmetricMatrix &A) { return A.g_->num_nodes(); }
inline std::size_t num_cols(const GraphSymmetricMatrix &A) { return A.g_->num_nodes(); }

/** Traits that MTL uses to determine properties of our GraphSymmetricMatrix . */
namespace mtl
{
namespace ashape
{
/** Define GraphSymmetricMatrix to be a non - scalar type . */
template <>
struct ashape_aux<GraphSymmetricMatrix>
{
  typedef nonscal type;
};
} // end namespace ashape

/** GraphSymmetricMatrix implements the Collection concept
 * with value_type and size_type */
template <>
struct Collection<GraphSymmetricMatrix>
{
  typedef double value_type;
  typedef unsigned size_type;
};
} // namespace mtl

/** Remove all the nodes in graph @a g whose position is within Box3D @a bb.
 * @param[in,out] g  The Graph to remove nodes from
 * @param[in]    bb  The BoundingBox, all nodes inside this box will be removed
 * @post For all i, 0 <= i < @a g.num_nodes(),
 *        not bb.contains(g.node(i).position())
 */
void remove_box(GraphType &g, const Box3D &bb)
{
  g_func G;
  std::vector<NodeType> eraselist;
  for (auto it = g.node_begin(); it != g.node_end(); ++it)
  {
    auto p = (*it).position();
    if (bb.contains(p))
    {
      eraselist.push_back(*it);
    };
  };
  for (auto it = eraselist.begin(); it != eraselist.end(); ++it)
  {
    g.remove_node(*it);
  }
};
double f(Point p)
{
  return 5 * cos(norm_1(p));
};

void boundary_tag(GraphType *g, std::vector<unsigned int> &boundary)
{
  Box3D B(Point(-0.6, -0.2, -1), Point(0.6, 0.2, 1));
  for (auto it = g->node_begin(); it != g->node_end(); ++it)
  {
    auto p = (*it).position();
    if ((norm_inf(p) == 1) ||
        (norm_inf(p - Point(0.6, 0.6, 0)) < 0.2) ||
        (norm_inf(p - Point(-0.6, 0.6, 0)) < 0.2) ||
        (norm_inf(p - Point(0.6, -0.6, 0)) < 0.2) ||
        (norm_inf(p - Point(-0.6, -0.6, 0)) < 0.2) ||
        B.contains(p))
    {
      boundary.push_back((*it).index());
    }
  };
};
struct NodeColor
{
  NodeColor(mtl::dense_vector<double> &x)
  {
    x_ = x;
    maxel = *std::max_element(x_.begin(), x_.end());
    minel = *std::min_element(x_.begin(), x_.end());
    range = std::abs(maxel - minel + 1e-5);
  }

  CME212::Color operator()(GraphType::node_type n1)
  {
    double val = x_[n1.index()];
    return (CME212::Color::make_heat((val - minel + 1e-5) / range));
  }
  mtl::dense_vector<double> x_;
  double maxel;
  double minel;
  double range;
};

struct NodePosition
{
  NodePosition(mtl::dense_vector<double> &x) { x_ = x; }

  Point operator()(GraphType::node_type n1)
  {
    Point pos = n1.position();
    pos.z = x_[n1.index()];
    return pos;
  }
  mtl::dense_vector<double> x_;
};
template <class Real>
class visual_iteration : public itl::cyclic_iteration<Real>
{
public:
  typedef itl::cyclic_iteration<Real> super;
  typedef visual_iteration self;

  visual_iteration(CME212::SFML_Viewer *viewer, GraphType *graph, mtl::dense_vector<double> *x,
                   mtl::dense_vector<double> *b, std::map<NodeType, unsigned> &map) : super(*b, 500, 1.e-10, 0.0, 50)
  {
    viewer_ = viewer;
    graph_ = graph;
    x_ = x;
    map_ = map;
  };

  void plot_iteration()
  {
    // viewer_->clear();
    // auto node_map = viewer_->empty_node_map(*graph_);
    NodeColor colorfunc(*x_);
    NodePosition positionfunc(*x_);
    viewer_->add_nodes(graph_->node_begin(), graph_->node_end(), colorfunc, positionfunc, map_);
    // viewer_->add_edges(graph_->edge_begin(), graph_->edge_end(), node_map);
  }

  template <typename T>
  bool finished(const T &r)
  {
    bool ret = super::finished(r);
    plot_iteration();
    return ret;
  }

private:
  CME212::SFML_Viewer *viewer_;
  GraphType *graph_;
  mtl::dense_vector<double> *x_;
  std::map<NodeType, unsigned> map_;
};

int main(int argc, char **argv)
{
  // Check arguments
  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }
  // Define an empty Graph
  GraphType graph;
  std::vector<unsigned int> boundary;
  {
    // Create a nodes_file from the first input argument
    std::ifstream nodes_file(argv[1]);
    // Interpret each line of the nodes_file as a 3D Point and add to the Graph
    std::vector<NodeType> node_vec;
    Point p;
    while (CME212::getline_parsed(nodes_file, p))
      node_vec.push_back(graph.add_node(2 * p - Point(1, 1, 0)));

    // Create a tets_file from the second input argument
    std::ifstream tets_file(argv[2]);
    // Interpret each line of the tets_file as four ints which refer to nodes
    std::array<int, 4> t;
    while (CME212::getline_parsed(tets_file, t))
    {
      graph.add_edge(node_vec[t[0]], node_vec[t[1]]);
      graph.add_edge(node_vec[t[0]], node_vec[t[2]]);
      graph.add_edge(node_vec[t[1]], node_vec[t[3]]);
      graph.add_edge(node_vec[t[2]], node_vec[t[3]]);
    }
  }

  // Get the edge length, should be the same for each edge
  auto it = graph.edge_begin();
  assert(it != graph.edge_end());
  double h = norm((*it).node1().position() - (*it).node2().position());

  // Make holes in our Graph
  remove_box(graph, Box3D(Point(-0.8 + h, -0.8 + h, -1), Point(-0.4 - h, -0.4 - h, 1)));
  remove_box(graph, Box3D(Point(0.4 + h, -0.8 + h, -1), Point(0.8 - h, -0.4 - h, 1)));
  remove_box(graph, Box3D(Point(-0.8 + h, 0.4 + h, -1), Point(-0.4 - h, 0.8 - h, 1)));
  remove_box(graph, Box3D(Point(0.4 + h, 0.4 + h, -1), Point(0.8 - h, 0.8 - h, 1)));
  remove_box(graph, Box3D(Point(-0.6 + h, -0.2 + h, -1), Point(0.6 - h, 0.2 - h, 1)));

  boundary_tag(&graph, boundary);

  std::cout << boundary.size() << std::endl;

  // HW3: YOUR CODE HERE
  // Define b using the graph, f, and g.
  // Construct the GraphSymmetricMatrix A using the graph
  // Solve Au = b using MTL.
  mtl::dense_vector<double> x(graph.num_nodes()), b(graph.num_nodes());
  g_func G;

  // Building b
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it)
  {
    auto i = (*it).index();
    auto xi = (*it).position();
    if (is_boundary(*it))
    {
      b[i] = G(xi);
    }
    else
    {
      double sum = 0;
      for (auto it2 = (*it).edge_begin(); it2 != (*it).edge_end(); ++it2)
      {
        auto node2 = (*it2).node2();
        if (is_boundary(node2))
        {
          sum += G(node2.position());
        }
      }
      b[i] = (h * h) * f(xi) - sum;
    }
  };
  GraphSymmetricMatrix A(graph, boundary);
  itl::pc::identity<GraphSymmetricMatrix> P(A);
  x = b;

  CME212::SFML_Viewer viewer;
  NodeColor colorfunc(x);
  NodePosition positionfunc(x);

  auto node_map = viewer.empty_node_map(graph);
  viewer.add_nodes(graph.node_begin(), graph.node_end(), colorfunc, positionfunc, node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
  viewer.center_view();

  // itl::cyclic_iteration<double> iter(b, 300, 1.e-15, 0.0, 50);
  visual_iteration<double> vis_iter(&viewer, &graph, &x, &b, node_map);

  bool interrupt_sim_thread = false;

  auto sim_thread = std::thread([&]() { cg(A, x, b, P, vis_iter); });
  viewer.event_loop();

  interrupt_sim_thread = true;
  sim_thread.join();

  return 0;
}
