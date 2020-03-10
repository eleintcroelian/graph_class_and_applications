/**
 * @file mass_spring.cpp
 * Implementation of mass-spring system using Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <fstream>
#include <chrono>
#include <thread>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"

// Gravity in meters/sec^2
static constexpr double grav = 9.81;
static constexpr double c = 0.01;

/** Custom structure of data to store with Nodes */
struct NodeData
{
  Point vel;   //< Node velocity
  double mass; //< Node mass
  NodeData() : vel(0), mass(1) {}
};

struct EdgeData
{
  double RestLength; //< Node mass
  EdgeData() : RestLength(0.) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 */
template <typename G, typename F, typename C>
double symp_euler_step(G &g, double t, double dt, F force, C constraint)
{

  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it)
  {
    auto n = *it;
    // if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))
    // {
    //   continue;
    // }
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }
  constraint(&g, t);
  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it)
  {
    auto n = *it;
    // if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))
    // {
    //   continue;
    // }

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }
  return t + dt;
}

/** Force function object for HW2 #1. */
struct Problem1Force
{
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t)
  {
    (void)t;
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))
    {
      return Point(0, 0, 0);
    }

    // std::cout<<"Node: "<< n.index()<< " connected to : ";
    Point f_spring(0, 0, 0);
    auto p1 = n.position();
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it)
    {
      auto current_edge = *it;
      auto n2 = current_edge.node2();
      // std::cout<< n2.index()<<" ";
      auto p2 = n2.position();
      Point direction = (p1 - p2) / norm(p2 - p1);
      double current_distance = norm(p2 - p1);
      f_spring = f_spring - 100 * direction * (current_distance - current_edge.value().RestLength);
    }
    // std::cout<< std::endl;
    Point f_gravity(0, 0, -1);
    return f_spring + (grav * n.value().mass) * f_gravity;
  }
};
class ZeroForce
{
public:
  virtual Point operator()(Node n, double t)
  {
    (void)t;
    (void)n;
    Point force(0, 0, 0);
    return force;
  }
};

class GravityForce : public ZeroForce
{
  virtual Point operator()(Node n, double t)
  {
    (void)t;
    Point f_gravity(0, 0, -1);
    return (grav * n.value().mass) * f_gravity;
  }
};
class MassSpringForce : public ZeroForce
{
  virtual Point operator()(Node n, double t)
  {
    (void)t;
    Point f_spring(0, 0, 0);
    auto p1 = n.position();
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it)
    {
      auto current_edge = *it;
      auto n2 = current_edge.node2();
      // std::cout<<n2.index()<<std::endl;
      auto p2 = n2.position();
      Point direction = (p1 - p2) / norm(p2 - p1);
      double current_distance = norm(p2 - p1);
      f_spring = f_spring - 100 * direction * (current_distance - current_edge.value().RestLength);
    }
    return f_spring;
  }
};

class DampingForce : public ZeroForce
{
  virtual Point operator()(Node n, double t)
  {
    (void)t;
    Point f_damping(1, 1, 1);
    return f_damping * (-c) * n.value().vel;
  }
};

struct CombinedForces
{
  CombinedForces(std::vector<ZeroForce *> inputforces) : inputforces_(inputforces){};

  Point operator()(Node n, double t)
  {
    (void)t;
    Point Sum(0, 0, 0);
    for (auto it = inputforces_.begin(); it != inputforces_.end(); it++)
    {
      Sum = Sum + (*(*it))(n, t);
    }
    return Sum;
  }
  std::vector<ZeroForce *> inputforces_;
};

CombinedForces make_combined_force(GravityForce gf, MassSpringForce mf)
{
  std::vector<ZeroForce *> vec;
  vec.push_back(&gf);
  vec.push_back(&mf);
  CombinedForces output(vec);
  return output;
};

class ZeroConstraint
{
public:
  virtual void operator()(GraphType *graph, double t)
  {
    (void)t;
    (void)graph;
  }
};

class PinConstraint : public ZeroConstraint
{
public:
  PinConstraint(GraphType *graph)
  {
    for (auto it = (*graph).node_begin(); it != (*graph).node_end(); ++it)
    {
      if ((*it).position() == Point(0, 0, 0))
      {
        node_1 = *it;
      }
      if ((*it).position() == Point(1, 0, 0))
      {
        node_2 = *it;
      }
    }
  };
  void operator()(GraphType *graph, double t)
  {
    (void)t;
    (void)graph;
    node_1.position() = Point(0, 0, 0);
    node_2.position() = Point(1, 0, 0);
    node_1.value().vel = Point(0, 0, 0);
    node_2.value().vel = Point(0, 0, 0);
  }
  Node node_1;
  Node node_2;
};

class PlaneConstraint : public ZeroConstraint
{
public:
  void operator()(GraphType *graph, double t)
  {
    (void)t;
    for (auto it = (*graph).node_begin(); it != (*graph).node_end(); ++it)
    {
      auto n = *it;
      if (n.position().z < -0.75)
      {
        n.position().z = -0.75;
        n.value().vel.z = 0;
      }
    }
  }
};
class SphereConstraint : public ZeroConstraint
{
public:
  void operator()(GraphType *graph, double t)
  {
    (void)t;
    for (auto it = (*graph).node_begin(); it != (*graph).node_end(); ++it)
    {
      auto n = *it;
      if (norm(n.position() - Point(0.5, 0.5, -0.5)) < 0.15)
      {
        Point center(0.5, 0.5, -0.5);
        Point R = (n.position() - center) / norm(center - n.position());
        n.position() = center + R * 0.15;
        n.value().vel = n.value().vel - (dot(n.value().vel, R) * R);
      }
    }
  }
};
class RemoveConstraint : public ZeroConstraint
{
public:
  void operator()(GraphType *graph, double t)
  {
    (void)t;
    std::vector<Node> eraselist;
    for (auto it = (*graph).node_begin(); it != (*graph).node_end(); ++it)
    {
      auto n = *it;
      if (norm(n.position() - Point(0.5, 0.5, -0.5)) < 0.15)
      {
        eraselist.push_back(n);
      }
    }
    for (auto el : eraselist)
    {
      // auto extnode = graph->find_external(el);
      // graph->remove_node(extnode);
      graph->remove_node(el);
    }
  }
};
struct CombinedConstraints
{
  CombinedConstraints(std::vector<ZeroConstraint *> inputconstraints) : inputconstraints_(inputconstraints){};

  void operator()(GraphType *g, double t)
  {
    (void)t;
    for (auto it = inputconstraints_.begin(); it != inputconstraints_.end(); it++)
    {
      (*(*it))(g, t);
    };
  }
  std::vector<ZeroConstraint *> inputconstraints_;
};

int main(int argc, char **argv)
{
  // Check arguments
  if (argc < 3)
  {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct an empty graph
  GraphType graph;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  std::vector<typename GraphType::node_type> nodes;
  while (CME212::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int, 4> t;
  while (CME212::getline_parsed(tets_file, t))
  {
    graph.add_edge(nodes[t[0]], nodes[t[1]]);
    graph.add_edge(nodes[t[0]], nodes[t[2]]);

    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);

    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.

  for (auto it = graph.node_begin(); it != graph.node_end(); ++it)
  {
    auto n = *it;
    n.value().mass = 1. / graph.num_nodes();
  }
  for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it)
  {
    auto e = *it;
    e.value().RestLength = e.length();
  }

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the Viewer
  CME212::SFML_Viewer viewer;
  auto node_map = viewer.empty_node_map(graph);

  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  viewer.center_view();

  // We want viewer interaction and the simulation at the same time
  // Viewer is thread-safe, so launch the simulation in a child thread
  bool interrupt_sim_thread = false;
  auto sim_thread = std::thread([&]() {
    // Begin the mass-spring simulation
    double dt = 0.0001;
    double t_start = 0;
    double t_end = 5.0;
    // double L = (*graph.edge_begin()).length();

    for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt)
    {
      //std::cout << "t = " << t << std::endl;
      std::vector<ZeroConstraint *> constraint_vector;
      PinConstraint p_c(&graph);
      SphereConstraint s_c;
      PlaneConstraint pl_c;
      RemoveConstraint r_c;
      constraint_vector.push_back(&p_c);
      // constraint_vector.push_back(&s_c);
      constraint_vector.push_back(&r_c);

      symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce()),
                      CombinedConstraints(constraint_vector));
      // viewer.clear();
      // node_map.clear();
      viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
      // viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
      // Update viewer with nodes' new positions
      viewer.set_label(t);
      // These lines slow down the animation for small graphs, like grid0_*.
      // Feel free to remove them or tweak the constants.
      if (graph.size() < 100)
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
  }); // simulation thread

  viewer.event_loop();
  // If we return from the event loop, we've killed the window.
  interrupt_sim_thread = true;
  sim_thread.join();

  return 0;
}
