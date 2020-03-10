#ifndef CME212_g_HPP
#define CME212_g_HPP

/** @file The peercode g_1481.hpp from hw2 have been highly taken as basis for
 * this Graph class for the purposes of hw3.
 * @brief An undirected graph type
 */

//--functionality_0
//--good job correcting Graph
//--END

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph
{
private:
  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

public:
  //
  // PUBLIC TYPE DEFINITIONS
  //

  typedef V node_value_type;
  typedef E edge_value_type;

  struct NodeContents
  {
    NodeContents(Point p, V v)
    {
      point = p;
      value = v;
    }
    Point point;
    V value;
  };
  struct EdgeContents
  {
    EdgeContents(E e, std::vector<unsigned int> n)
    {
      value = e;
      nodes = n;
    }
    E value;
    std::vector<unsigned> nodes;
  };

  /** Type of this graph. */
  using g_type = Graph;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  using node_iterator = NodeIterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  using edge_iterator = EdgeIterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  using incident_iterator = IncidentIterator;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
  {
  }

  /** Default destructor */
  ~Graph() = default;

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node : private totally_ordered<Node>
  {
  public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Graph::node_type x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */

    Node()
    {
    }

    /** Return this node's position. */
    const Point &position() const
    {
      return g_->NodeVector[idx_].point;
    }

    Point &position()
    {
      return g_->NodeVector[idx_].point;
    }

    /** Return this node's index, a number in the range [0, g_size). */
    size_type index() const
    {
      return g_->node_u2i_[idx_];
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /* @return This node's value. */
    node_value_type &value()
    {
      return g_->NodeVector[idx_].value;
    }
    /* @return This node's value. */
    const node_value_type &value() const
    {
      return g_->NodeVector[idx_].value;
    }

    /* @return The number of incident edges.*/
    size_type degree() const
    {
      return g_->EdgeMap[idx_].size();
    }

    /* @return An iterator pointing 
     * at the start of the edges incident to the node 
     */
    incident_iterator edge_begin() const
    {
      return incident_iterator(g_, idx_, true);
    }

    /* @return An iterator pointing 
     * at the end of the edges incident to the node
     */
    incident_iterator edge_end() const
    {
      return incident_iterator(g_, idx_, false);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node &n) const
    {
      if ((n.g_ == g_) && (n.idx_ == idx_))
        return true;
      else
        return false;
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node &n) const
    {
      std::less<Graph *> node_less;
      if (node_less(g_, n.g_))
        return true;
      else if ((g_ == n.g_) && (idx_ < n.idx_))
        return true;
      else
        return false;
    }

  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph *g_;
    size_type idx_; //in my code, all idx_ actually works in the way of uid
    Node(const Graph *graph, size_type idx)
        : g_(const_cast<Graph *>(graph)), idx_(idx)
    {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const
  {
    return node_i2u_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const
  {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point &position,
                const node_value_type &value = node_value_type())
  {
    NodeVector.push_back(NodeContents(position, value));
    node_i2u_.push_back(node_u2i_.size());
    node_u2i_.push_back(node_i2u_.size() - 1);
    return Node(this, NodeVector.size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node &n) const
  {
    if ((n.index() < node_i2u_.size()) && (n.g_ == this))
      return true;
    else
      return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const
  {
    assert(i >= 0);
    assert(i < num_nodes());
    return Node(this, node_i2u_[i]);
  }

  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge : private totally_ordered<Edge>
  {
  public:
    /** Construct an invalid Edge. */
    Edge()
    {
    }

    //HW2:
    /* @return This node's value. */
    edge_value_type &value()
    {
      return g_->EdgeVector[idx_].value;
    }
    /* @return This node's value. */
    const edge_value_type &value() const
    {
      return g_->EdgeVector[idx_].value;
    }

    /** Return a node of this Edge */
    Node node1() const
    {
      return g_->node(g_->node_u2i_[idx1_]);
    }

    /** Return the other node of this Edge */
    Node node2() const
    {
      return g_->node(g_->node_u2i_[idx2_]);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge &e) const
    {
      if ((e.g_ == g_) &&
          (e.node1() == node1()) && (e.node2() == node2()))
        return true;
      else if ((e.g_ == g_) &&
               (e.node2() == node1()) && (e.node1() == node2()))
        return true;
      else
        return false;
    }

    double length() const
    {
      Point diff = node1().position() - node2().position();
      double len = norm(diff);
      return len;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */

    bool operator<(const Edge &e) const
    {
      std::less<Graph *> edge_less;
      if (edge_less(g_, e.g_))
        return true;
      else if ((g_ == e.g_) && (idx_ < e.idx_))
        return true;
      else
        return false;
    }

  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph *g_;
    size_type idx1_; //node 1's uid
    size_type idx2_; //node 2's uid
    size_type idx_;  //edge's uid
    Edge(const Graph *graph, size_type idx1, size_type idx2, size_type idx)
        : g_(const_cast<Graph *>(graph)), idx1_(idx1), idx2_(idx2), idx_(idx)
    {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const
  {
    return edge_i2u_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const
  {
    assert(i >= 0);
    assert(i < EdgeVector.size());
    return Edge(this, EdgeVector[edge_i2u_[i]].nodes[0], EdgeVector[edge_i2u_[i]].nodes[1], edge_i2u_[i]);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node &a, const Node &b) const
  {
    assert(has_node(a));
    assert(has_node(b));

    auto iter = EdgeMap.find(node_i2u_[a.index()]);

    if (iter == EdgeMap.end())
    {
      return false;
    }
    else if (iter->second.find(node_i2u_[b.index()]) == iter->second.end())
    {
      return false;
    }
    return true;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node &a, const Node &b, const edge_value_type &value = edge_value_type())
  {
    assert(has_node(a));
    assert(has_node(b));
    assert(a.index() != b.index());

    size_type a_uid = node_i2u_[a.index()];
    size_type b_uid = node_i2u_[b.index()];

    if (has_edge(a, b))
    {
      size_type temp_idx = EdgeMap[a_uid][b_uid];
      return Edge(this, a_uid, b_uid, temp_idx);
    }

    std::vector<size_type> temp = {a_uid, b_uid};
    EdgeVector.push_back(EdgeContents(value, temp));
    edge_i2u_.push_back(edge_u2i_.size());
    edge_u2i_.push_back(edge_i2u_.size() - 1);
    EdgeMap[a_uid][b_uid] = EdgeVector.size() - 1;
    EdgeMap[b_uid][a_uid] = EdgeVector.size() - 1;
    return Edge(this, a_uid, b_uid, EdgeVector.size() - 1);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear()
  {
    EdgeMap.clear();
    EdgeVector.clear();
    NodeVector.clear();
    node_i2u_.clear();
    node_u2i_.clear();
    edge_i2u_.clear();
    edge_u2i_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator
  {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Node;                           // Element type
    using pointer = Node *;                            // Pointers to elements
    using reference = Node &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator()
    {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /* @brief Dereference operator
     * @return An node pointed by the iterator
     */
    Node operator*() const
    {
      return itr_g_->node(idx_);
    }

    /* @brief Increment to the next node
     * @return A node iterator pointing to the next node
     * @post The iterator points to the next node
     */
    NodeIterator &operator++()
    {
      idx_++;
      return *this;
    }

    /* @brief Defines equality between two iterators
     * @return True if two iterators points to the same node
     */
    bool operator==(const NodeIterator &node_itr) const
    {
      return ((itr_g_ == node_itr.itr_g_) && (idx_ == node_itr.idx_));
    }

    /* @brief Defines inequality between two iterators
     * @return True if two iterators don't point to the same node
     */
    bool operator!=(const NodeIterator &node_itr) const
    {
      return ((itr_g_ != node_itr.itr_g_) || (idx_ != node_itr.idx_));
    }

  private:
    friend class Graph;
    Graph *itr_g_;
    size_type idx_;
    NodeIterator(const Graph *itr_graph, size_type id)
        : itr_g_(const_cast<Graph *>(itr_graph)), idx_(id)
    {
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /* 
   * @return An iterator that points at the start of the nodes
   */
  node_iterator node_begin() const
  {
    return node_iterator(this, 0);
  }

  /*
   * @return An iterator that points at end end of the nodes
   */
  node_iterator node_end() const
  {
    return node_iterator(this, num_nodes());
  }
  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator
  {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                           // Element type
    using pointer = Edge *;                            // Pointers to elements
    using reference = Edge &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator()
    {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /* @brief Dereference operator
     * @return An edge pointed by the iterator
     */
    Edge operator*() const
    {
      size_type node2_idx_ = map_itr_->first;
      size_type edge_idx_ = map_itr_->second;
      return Edge(ii_g_, node_idx_, node2_idx_, edge_idx_);
    }

    /* @brief Increment to the next edge incident to the node
     * @return An edge iterator pointing to the next edge incident to the node
     * @post The iterator points to the next edge incident to the node
     */
    IncidentIterator &operator++()
    {
      map_itr_++;
      return *this;
    }

    /* @brief Defines equality between two iterators
     * @return True if two iterators points to the same edge 
     */
    bool operator==(const IncidentIterator &inc_itr) const
    {
      return ((ii_g_ == inc_itr.ii_g_) && (node_idx_ == inc_itr.node_idx_) && (map_itr_ == inc_itr.map_itr_));
    }

    /* @brief Defines inequality between two iterators
     * @return True if two iterators don't point to the same edge
     */
    bool operator!=(const IncidentIterator &inc_itr) const
    {
      return ((ii_g_ != inc_itr.ii_g_) || (node_idx_ != inc_itr.node_idx_) || (map_itr_ != inc_itr.map_itr_));
    }

  private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph *ii_g_;
    size_type node_idx_;
    std::map<size_type, size_type>::iterator map_itr_;
    IncidentIterator(const Graph *ii_graph, size_type node_idx, bool begin_)
        : ii_g_(const_cast<Graph *>(ii_graph)), node_idx_(node_idx)
    {
      if (begin_ == true)
      {
        map_itr_ = ii_g_->EdgeMap[node_idx].begin();
      }
      else
      {
        map_itr_ = ii_g_->EdgeMap[node_idx].end();
      }
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator
  {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                           // Element type
    using pointer = Edge *;                            // Pointers to elements
    using reference = Edge &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator()
    {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /* @brief Dereference operator
     * @return An edge pointed by the iterator
     */
    Edge operator*() const
    {
      return g_->edge(edge_idx_);
    }

    /* @brief Increment to the next edge
     * @return An edge iterator pointing to the next edge
     * @post The iterator points to the next edge 
     */
    EdgeIterator &operator++()
    {
      edge_idx_++;
      return *this;
    }

    /* @brief Defines equality between two iterators 
     * @return True if two iterators points to the same edge
     */
    bool operator==(const EdgeIterator &ei) const
    {
      return ((g_ == ei.g_) && (edge_idx_ == ei.edge_idx_));
    }

    /* @brief Defines inequality between two iterators 
     * @return True if two iterators don't point to the same edge 
     */
    bool operator!=(const EdgeIterator &ei) const
    {
      return ((g_ != ei.g_) || (edge_idx_ != ei.edge_idx_));
    }

  private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph *g_;
    size_type edge_idx_;
    EdgeIterator(const Graph *ei_graph, size_type edge_idx)
        : g_(const_cast<Graph *>(ei_graph)), edge_idx_(edge_idx)
    {
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /* @return An iterator pointing at the start of the edges */
  edge_iterator edge_begin() const
  {
    return edge_iterator(this, 0);
  }

  /* @return An iterator pointing at the end of the edges */
  edge_iterator edge_end() const
  {
    return edge_iterator(this, num_edges());
  }

  /** Remove an Edge from the graph, and return the number of edge removed (0 or 1).
   * @pre edge is currently inside the graph
   * @return 1 if node removed, 0 if not.
   * @post has_edge(@a n1 , @a n2) == true
   * @post num.edges()-= return value
   *       
   *
   * Order of edge content data is changed.
   * 
   *
   * Complexity: O(1)
   */
  size_type remove_edge(const Node &n1, const Node &n2)
  {
    if (has_edge(n1, n2))
    {
      auto n1_uid = node_i2u_[n1.index()];
      auto n2_uid = node_i2u_[n2.index()];
      auto rm = EdgeMap[n1_uid][n2_uid];
      edge_u2i_[edge_i2u_.back()] = edge_u2i_[rm];
      std::swap(edge_i2u_[edge_u2i_[rm]], edge_i2u_.back());
      edge_i2u_.pop_back();
      EdgeMap[n1_uid].erase(n2_uid);
      EdgeMap[n2_uid].erase(n1_uid);
      return 1;
    }
    else
      return 0;
  }

  /** Remove an Edge from the graph, and return the number of edge removed (0 or 1).
   * @pre edge is currently inside the graph
   * @return 1 if node removed, 0 if not.
   * @post has_edge(@a n1 , @a n2) == true
   * @post num.edges()-= return value
   *       
   *
   * Order of edge content data is changed.
   * 
   *
   * Complexity: O(1)
   */
  size_type remove_edge(const Edge &e)
  {
    Node n1 = e.node1();
    Node n2 = e.node2();
    remove_edge(n1, n2);
    return 1;
  }

  /** Remove an Edge from the graph, and return an iterator to the edge removed.
   * @pre edge is currently inside the graph
   * @return 1 if node removed, 0 if not.
   * @post has_edge(@a n1 , @a n2) == true
   * @post num.edges()-= return value
   *       
   *
   * Order of edge content data is changed.
   * 
   *
   * Complexity: O(1)
   */
  edge_iterator remove_edge(edge_iterator it)
  {
    Node n1 = (*it).node1();
    Node n2 = (*it).node2();
    remove_edge(n1, n2);
    return it;
  }

  /** Remove a Node from the graph, and return the number of node removed (0 or 1).
   * @pre @a n is currently inside the graph
   * @return 1 if node removed, 0 if not.
   * @post has_node(@a n) == true
   * @post num.nodes()-= return value
   *       
   *
   * Order of node content data is changed. Removes all edges connected to the node.
   * 
   *
   * Complexity: O(num_nodes() + num_edges())
   */
  size_type remove_node(const Node &n)
  {
    if (has_node(n))
    {
      while (true)
      {
        auto it = n.edge_begin();
        if (it == n.edge_end())
        {
          break;
        }
        else
        {
          remove_edge(*it);
        }
      }
      EdgeMap.erase(node_i2u_[n.index()]);
      node_u2i_[node_i2u_.back()] = n.index();
      std::swap(node_i2u_.back(), node_i2u_[n.index()]);
      node_i2u_.pop_back();
      return 1;
    }
    else
    {
      return 0;
    }
  }

  /** Remove a Node from the graph, and return the iterator to the node.
   * @pre @a n is currently inside the graph
   * @return an iterator to n.
   * @post has_node(@a n) == true
   * @post num.nodes()-= return value
   *       
   *
   * Order of node content data is changed. Removes all edges connected to the node.
   * 
   *
   * Complexity: O(num_nodes() + num_edges())
   */
  node_iterator remove_node(node_iterator n_it)
  {
    //--functionality_0
    //--small error
    node_iterator temp(this, (*n_it).index());
    //node_iterator temp(this, *n_it.index());
    //--END
    remove_node(*n_it);
    return temp;
  }

private:
  std::vector<size_type> node_i2u_;
  std::vector<size_type> node_u2i_;
  std::vector<size_type> edge_i2u_;
  std::vector<size_type> edge_u2i_;
  std::vector<NodeContents> NodeVector;
  std::vector<EdgeContents> EdgeVector;
  std::map<size_type, std::map<size_type, size_type>> EdgeMap;
};

#endif // CME212_g_HPP
