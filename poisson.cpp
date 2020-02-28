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
    double Lij;
    double Aij;
    for (size_t i = 0; i < g_->num_nodes(); i++) //row
    {
      unsigned int temp = 0;
      for (size_t j = 0; j < g_->num_nodes(); j++) //column
      {
        if (i == j)
        {
          Lij = -(double)g_->node(i).degree();
        }
        else if (g_->has_edge(g_->node(i), g_->node(j)))
        {
          Lij = 1;
        }
        else
        {
          Lij = 0;
        }
        {if ((i == j) && (is_boundary(g_->node(i))))
        {
          Aij = 1;
        }
        if ((i != j) && (is_boundary(g_->node(i)) || is_boundary(g_->node(j))))
        {
          Aij = 0;
        }
        else
        {
          Aij = Lij;
        }}
        if (Aij == 0)
        {
          continue;
        }
        temp += Aij * v[j];
        // A_row.push_back(Aij);
      }
      // temp = std::inner_product(A_row.begin(), A_row.end(), v.begin(), 0);
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
  for(auto it = eraselist.begin(); it != eraselist.end(); ++it)
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
  mtl::dense_vector<double> x(graph.num_nodes(), 1.0), b(graph.num_nodes());
  g_func G;

  // Building b
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it)
  {
    auto i = (*it).index();
    auto xi = (*it).position();
    if (boundary.end() != find(boundary.begin(), boundary.end(), i))
    {
      b[i] = G(xi);
    }
    else
    {
      double sum = 0;
      for (auto it2 = (*it).edge_begin(); it2 != (*it).edge_end(); ++it2)
      {
        auto node2 = (*it2).node2();
        sum += G(node2.position());
      }
      b[i] = (h * h) * f(xi) - sum;
    }
  };
  GraphSymmetricMatrix A(graph, boundary);
  itl::pc::identity<GraphSymmetricMatrix> P(A);
  x = 0;
  // std::cout << "b is " << b << std::endl;
  // std::cout << "x is " << x << std::endl;
  itl::cyclic_iteration<double> iter(b, 100, 1.e-11, 0.0, 5);
  cg(A, x, b, P, iter);
  // std::cout << "b is " << b << std::endl;
  // std::cout << "x is " << x << std::endl;
  return 0;
}
