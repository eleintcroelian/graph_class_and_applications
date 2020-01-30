#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V>
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

  /** Type of this graph. */
  using node_value_type = V;
  using graph_type = Graph;

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

  struct NodeContents
  {
    Point _Point;
    V _Value;
    node_type _Node;
  };

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
  {
    std::vector<NodeContents> _Node_Struct_Vector;
    std::map<size_type, std::map<size_type, size_type>> _n1_n2_edge;
    std::vector<edge_type> _Edge_Vector;
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
    size_type degree() const
    {
      return _graph_pointer->_n1_n2_edge[_node_id].size();
    }
    incident_iterator edge_begin()
    {
      incident_iterator newiterator(const_cast<Graph *>(_graph_pointer), _node_id);
      return newiterator;
    }
    incident_iterator edge_end()
    {
      incident_iterator newiteratorend(const_cast<Graph *>(_graph_pointer), _node_id);
      newiteratorend._map_iter = (newiteratorend._map_pointer->end());
      newiteratorend._edge_id = newiteratorend._map_iter->second;
      newiteratorend._edge_pointer = &(newiteratorend._graph_pointer->_Edge_Vector[newiteratorend._edge_id]);

      return newiteratorend;
    }
    /** Return this node's position. */
    const Point &position() const
    {
      return _graph_pointer->_Node_Struct_Vector[_node_id]._Point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const
    {
      return _node_id;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;

    node_value_type &value()
    {
      return const_cast<Graph *>(_graph_pointer)->_Node_Struct_Vector[_node_id]._Value;
    };

    const node_value_type &value() const
    {
      const node_value_type constval = _graph_pointer->_Node_Struct_Vector[_node_id]._Value;
      return constval;
    };
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node &n) const
    {
      if ((_node_id == n._node_id) && (_graph_pointer == n._graph_pointer))
      {
        return true;
      }
      else
      {
        return false;
      }
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
      if ((_node_id < n._node_id))
      {
        return true;
      }
      else
      {
        return false;
      }
    }

  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    Node(const Graph *pointer, size_type id)
        : _graph_pointer(const_cast<Graph *>(pointer)), _node_id(id)
    {
    }
    const Graph *_graph_pointer;
    size_type _node_id;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const
  {
    return _Node_Struct_Vector.size();
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
  Node add_node(const Point &position, const node_value_type &Val = node_value_type())
  {
    // if (!has_node(new_node))       // Uncommented since instructions doesn't specify whether
    // {                              // it is required to check if there exists a node with same position
    Node new_node = Node(this, _Node_Struct_Vector.size());
    NodeContents nodecontents;
    nodecontents._Point = position;
    nodecontents._Node = new_node;
    nodecontents._Value = Val;
    _Node_Struct_Vector.push_back(nodecontents);

    return new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node &n) const
  {
    if ((n._node_id < this->_Node_Struct_Vector.size()) && (n._graph_pointer == this))
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  void set_value(size_type index, V val)
  {
    _Node_Struct_Vector[index]._Value = val;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const
  {
    Node temp_node = Node(this, i);
    return temp_node;
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

    /** Return a node of this Edge */
    Node node1() const
    {
      return _edge_graph_pointer->_Node_Struct_Vector[_n1_id]._Node;
    }

    /** Return the other node of this Edge */
    Node node2() const
    {
      return _edge_graph_pointer->_Node_Struct_Vector[_n2_id]._Node;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge &e) const
    {
      if ((_n1_id == e._n1_id) && (_n2_id == e._n2_id))
      {
        return true;
      }
      else
      {
        return false;
      }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge &e) const
    {
      if (_edge_id < e._edge_id)
      {
        return true;
      }
      else
      {
        return false;
      }
    }

  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    const Graph *_edge_graph_pointer; // Pointer to the graph
    size_type _edge_id;               // edge id
    size_type _n1_id;                 // id of node 1
    size_type _n2_id;                 // id of node 2
    Edge(const Graph *pointer, size_type edgeid, size_type n1id, size_type n2id)
        : _edge_graph_pointer(const_cast<Graph *>(pointer)), _edge_id(edgeid), _n1_id(n1id), _n2_id(n2id)
    {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const
  {
    return _Edge_Vector.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const
  {
    return _Edge_Vector[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node &a, const Node &b) const
  {
    if (a < b)
    {
      if (_n1_n2_edge.find(a._node_id) != _n1_n2_edge.end())
      {
        if (_n1_n2_edge.at(a._node_id).find(b._node_id) != _n1_n2_edge.at(a._node_id).end())
        {
          return true;
        }
      }
      return false;
    }
    if (b < a)
    {
      if (_n1_n2_edge.find(b._node_id) != _n1_n2_edge.end())
      {
        if (_n1_n2_edge.at(b._node_id).find(a._node_id) != _n1_n2_edge.at(b._node_id).end()) //_____________
        {
          return true;
        }
      }
      return false;
    }
    return false;
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
  Edge add_edge(const Node &a, const Node &b)
  {
    if (!has_edge(a, b))
    {
      if (a < b)
      {
        edge_type newEdge = Edge(this, _Edge_Vector.size(), a._node_id, b._node_id);

        _n1_n2_edge.insert(std::make_pair(a._node_id, std::map<size_type, size_type>()));
        _n1_n2_edge[a._node_id].insert(std::make_pair(b._node_id, newEdge._edge_id));
        _Edge_Vector.push_back(newEdge);
        return newEdge;
      }
      if (!(a < b))
      {
        edge_type newEdge = Edge(this, _Edge_Vector.size(), b._node_id, a._node_id);
        _n1_n2_edge.insert(std::make_pair(b._node_id, std::map<size_type, size_type>()));
        _n1_n2_edge[b._node_id].insert(std::make_pair(a._node_id, newEdge._edge_id));
        _Edge_Vector.push_back(newEdge);
        return newEdge;
      }
    }
    size_type existingedgeid = _n1_n2_edge[a._node_id][b._node_id]; // If edge already exists
    return _Edge_Vector[existingedgeid];                            // return it
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear()
  {
    _Node_Struct_Vector.clear();
    _n1_n2_edge.clear();
    _Edge_Vector.clear();
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
    using value_type = Node;                             // Element type
    using pointer = Node *;                              // Pointers to elements
    using reference = Node &;                            // Reference to elements
    using difference_type = std::ptrdiff_t;              // Signed difference
    using iterator_category = std::forward_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator()
    {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Node operator*() const
    {
      // std::cout<<(*_node_pointer).index()<<std::endl;
      return *_node_pointer;
    }

    NodeIterator &operator++()
    {
      _node_pointer = &(_graph_pointer->_Node_Struct_Vector[++_id]._Node);
      return *this;
    }
    bool operator==(const NodeIterator &iter) const
    {
      if (_id == iter._id)
      {
        return true;
      }
      else
      {
        return false;
      }
    }
    bool operator!=(const NodeIterator &iter) const
    {
      if (_id != iter._id)
      {
        return true;
      }
      else
      {
        return false;
      }
    }

  private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    NodeIterator(Graph *graph_pointer, size_type idinput)
    {
      //const_cast<Graph *>
      _id = idinput;
      _graph_pointer = graph_pointer;
      _node_pointer = &(_graph_pointer->_Node_Struct_Vector[_id]._Node);
    }
    node_type *_node_pointer;
    Graph *_graph_pointer;
    size_type _id;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const
  {
    NodeIterator newiter(const_cast<Graph *>(this), 0);
    return newiter;
  }
  node_iterator node_end() const
  {
    NodeIterator enditer(const_cast<Graph *>(this), _Node_Struct_Vector.size() - 1);

    // newiter._id = _Node_Struct_Vector.size();//const_cast<Node *>
    // newiter._node_pointer = &((_Node_Struct_Vector[newiter._id]._Node));
    // std::cout << newiter._id << std::endl;
    return enditer;
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
    edge_type *_edge_pointer;
    size_type _edge_id;
    std::map<size_type, size_type> *_map_pointer;
    std::map<size_type, size_type>::iterator _map_iter;
    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const
    bool operator==(const IncidentIterator &iter) const
    {
      if ((_edge_id == iter._edge_id))
      {
        return true;
      }
      return false;
    }

    bool operator!=(const IncidentIterator &iter) const
    {
      // std::cout << "!= called" << std::endl;

      if ((_edge_id != iter._edge_id))
      {
        return true;
      }
        return false;
    }
    Edge operator*() const
    {
      // std::cout << "* called" << std::endl;
      return *_edge_pointer;
    }

    IncidentIterator &operator++()
    {
      // std::cout << "+ called" << std::endl;

      ++_map_iter;
      _edge_id = _map_iter->second;
      _edge_pointer = &(_graph_pointer->_Edge_Vector[_edge_id]);
      return *this;
    }

  private:
    friend class Graph;
    IncidentIterator(Graph *graph_pointer, size_type node_id)
    {
      _graph_pointer = graph_pointer;
      _node_id = node_id;
      _map_pointer = &(_graph_pointer->_n1_n2_edge[_node_id]);
      _map_iter = _map_pointer->begin();
      _edge_id = _map_iter->second;
      _edge_pointer = &(_graph_pointer->_Edge_Vector[_edge_id]);
    }
    // HW1 #3: YOUR CODE HERE
    Graph *_graph_pointer;
    size_type _node_id;
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
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const
    Edge operator*() const
    {
      // std::cout<<(*_node_pointer).index()<<std::endl;
      return *_edge_pointer;
    }

    EdgeIterator &operator++()
    {
      _edge_pointer = &(_graph_pointer->_Edge_Vector[++_id]);
      return *this;
    }
    bool operator==(const EdgeIterator &iter) const
    {
      if (_id == iter._id)
      {
        return true;
      }
      else
      {
        return false;
      }
    }
    bool operator!=(const EdgeIterator &iter) const
    {
      if (_id != iter._id)
      {
        return true;
      }
      else
      {
        return false;
      }
    }

  private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    EdgeIterator(Graph *graph_pointer, size_type idinput)
    {
      //const_cast<Graph *>
      _id = idinput;
      _graph_pointer = graph_pointer;
      _edge_pointer = &(_graph_pointer->_Edge_Vector[_id]);
    }
    edge_type *_edge_pointer;
    Graph *_graph_pointer;
    size_type _id;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  edge_iterator edge_begin() const
  {
    EdgeIterator newiter(const_cast<Graph *>(this), 0);
    return newiter;
  }
  edge_iterator edge_end() const
  {
    EdgeIterator enditer(const_cast<Graph *>(this), _Edge_Vector.size() - 1);
    return enditer;
  }

private:
  // Use this space for your Graph class's internals:
  // helper functions, data members, and so forth.

  std::vector<NodeContents> _Node_Struct_Vector;
  std::map<size_type, std::map<size_type, size_type>> _n1_n2_edge;
  // Nested map with node1_id - [node2_id - edge_id] structure
  std::vector<edge_type> _Edge_Vector; // Vector holding Edge objects (ordered)
};

#endif // CME212_GRAPH_HPP
