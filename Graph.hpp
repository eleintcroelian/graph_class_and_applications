#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP
#define BUGPRINT(x)                                   \
  do                                                  \
  {                                                   \
    std::cerr << __FILE__ << ":" << __LINE__ << ": "; \
    std::cerr << #x << " -> " << (x) << std::endl;    \
  } while (0)
/** @file Graph.hpp
 * @brief An undirected graph type
 */

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

  /** Type of this graph. */
  // using node_value_type = V;
  typedef V node_value_type;
  typedef E edge_value_type;

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
  struct EdgeContents
  {
    E _Value;
    edge_type _Edge;
  };

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
  {
    std::vector<NodeContents> _Node_Struct_Vector;
    std::map<size_type, std::map<size_type, size_type>> _n1_n2_edge;
    std::vector<EdgeContents> _Edge_Struct_Vector;
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

    Point &position()
    {
      return const_cast<Graph *>(_graph_pointer)->_Node_Struct_Vector[_node_id]._Point;
    };

    size_type degree() const
    {
      if (_graph_pointer->_n1_n2_edge.find(_node_id) != _graph_pointer->_n1_n2_edge.end())
      {
        return _graph_pointer->_n1_n2_edge.at(_node_id).size();
      }
      return 0;
    }
    incident_iterator edge_begin() const
    {
      incident_iterator newiterator(const_cast<Graph *>(_graph_pointer), _node_id, true);
      return newiterator;
    }
    incident_iterator edge_end() const
    {
      incident_iterator newiteratorend(const_cast<Graph *>(_graph_pointer), _node_id, false);
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
      return const_cast<Graph *>(_graph_pointer)->_Node_Struct_Vector[_graph_pointer->e2i_node[_node_id]]._Value;
    };

    const node_value_type &value() const
    {
      const node_value_type constval = _graph_pointer->_Node_Struct_Vector[_graph_pointer->e2i_node[_node_id]]._Value;
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
    return e2i_node.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const
  {
    return e2i_node.size();
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
    Node new_node = Node(this, _Node_Struct_Vector.size());
    NodeContents nodecontents;
    nodecontents._Point = position;
    nodecontents._Node = new_node;
    nodecontents._Value = Val;
    e2i_node.push_back(_Node_Struct_Vector.size());
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
    if ((n._node_id < this->e2i_node.size()) && (n._graph_pointer == this))
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
    _Node_Struct_Vector[e2i_node[index]]._Value = val;
  }
  void set_edge_value(size_type index, E val)
  {
    _Edge_Struct_Vector[e2i[index]]._Value = val;
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
    edge_value_type &value()
    {
      return const_cast<Graph *>(_edge_graph_pointer)->_Edge_Struct_Vector[_edge_id]._Value;
    };

    const edge_value_type &value() const
    {
      const edge_value_type constval = _edge_graph_pointer->_Edge_Struct_Vector[_edge_id]._Value;
      return constval;
    };
    double length() const
    {
      return norm((*this).node1().position() - (*this).node2().position());
    }
    /** Return a node of this Edge : EXTERNAL USE*/
    Node node1() const
    {
      return _edge_graph_pointer->_Node_Struct_Vector[_edge_graph_pointer->e2i_node[_n1_id]]._Node;
    }
    /** Return a node of this Edge : INTERNAL USE*/

    Node Node1_int() const
    {
      return _edge_graph_pointer->_Node_Struct_Vector[_n1_id]._Node;
    }

    /** Return the other node of this Edge */
    Node node2() const
    {
      return _edge_graph_pointer->_Node_Struct_Vector[_edge_graph_pointer->e2i_node[_n2_id]]._Node;
    }
    /** Return a node of this Edge : INTERNAL USE*/

    Node Node2_int() const
    {
      return _edge_graph_pointer->_Node_Struct_Vector[_n2_id]._Node;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge &e) const
    {
      if ((_n1_id == e._n1_id) && (_n2_id == e._n2_id) && (_edge_graph_pointer == e._edge_graph_pointer))
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
      if (_edge_graph_pointer < e._edge_graph_pointer)
      {
        return true;
      }
      if (_edge_graph_pointer == e._edge_graph_pointer)
      { //Reference: Graph_1666.hpp
        if (_edge_id < e._edge_id)
        {
          return true;
        }
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
    return e2i.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const
  {
    return _Edge_Struct_Vector[e2i[i]]._Edge;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node &a, const Node &b) const
  {
    Node A(this, e2i_node[a.index()]);
    Node B(this, e2i_node[b.index()]);
    if (_n1_n2_edge.find(A._node_id) != _n1_n2_edge.end())
    {
      if (_n1_n2_edge.at(A._node_id).find(B._node_id) != _n1_n2_edge.at(A._node_id).end())
      {
        return true;
      }
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
    if (!Has_edge(a, b))
    {
      if (a < b)
      {
        edge_type newEdge = Edge(this, _Edge_Struct_Vector.size(), a._node_id, b._node_id);

        _n1_n2_edge.insert(std::make_pair(a._node_id, std::map<size_type, size_type>()));
        _n1_n2_edge[a._node_id].insert(std::make_pair(b._node_id, newEdge._edge_id));
        EdgeContents newstruct;
        newstruct._Edge = newEdge;
        _Edge_Struct_Vector.push_back(newstruct);

        if (_n1_n2_edge.find(b._node_id) == _n1_n2_edge.end())
        {
          _n1_n2_edge.insert(std::make_pair(b._node_id, std::map<size_type, size_type>()));
        }
        _n1_n2_edge[b._node_id].insert(std::make_pair(a._node_id, newEdge._edge_id));
        e2i.push_back(_Edge_Struct_Vector.size() - 1);
        return newEdge;
      }
      if (!(a < b))
      {
        edge_type newEdge = Edge(this, _Edge_Struct_Vector.size(), b._node_id, a._node_id);

        _n1_n2_edge.insert(std::make_pair(b._node_id, std::map<size_type, size_type>()));
        _n1_n2_edge[b._node_id].insert(std::make_pair(a._node_id, newEdge._edge_id));
        EdgeContents newstruct;
        newstruct._Edge = newEdge;
        _Edge_Struct_Vector.push_back(newstruct);
        if (_n1_n2_edge.find(a._node_id) == _n1_n2_edge.end())
        {
          _n1_n2_edge.insert(std::make_pair(a._node_id, std::map<size_type, size_type>()));
        }
        _n1_n2_edge[a._node_id].insert(std::make_pair(b._node_id, newEdge._edge_id));
        e2i.push_back(_Edge_Struct_Vector.size() - 1);
        return newEdge;
      }
    }
    size_type existingedgeid = _n1_n2_edge[a._node_id][b._node_id]; // If edge already exists
    return _Edge_Struct_Vector[existingedgeid]._Edge;               // return it
  }

  size_type remove_node(const Node &n)
  {
    Node InternalVersion(this, e2i_node[n.index()]);
    if (has_node(n))
    {
      if (InternalVersion.degree() != 0)
      {
        std::vector<Edge> eraselist;
        for (auto it = InternalVersion.edge_begin(); it != InternalVersion.edge_end(); ++it)
        {
          eraselist.push_back(*it);
        }
        size_type counter = 0;
        for (auto it = eraselist.begin(); it != eraselist.end(); ++it)
        {
          auto temp = *it;
          auto intnode1 = temp.Node1_int();
          auto intnode2 = temp.Node2_int();
          if (Has_edge(intnode1, intnode2))
          {
            auto edge_id = _n1_n2_edge.at(intnode1.index()).at(intnode2.index());
            _n1_n2_edge.at(intnode1.index()).erase(intnode2.index());
            _n1_n2_edge.at(intnode2.index()).erase(intnode1.index());
            auto external_edge_id = find(e2i.begin(), e2i.end(), edge_id) - e2i.begin();
            std::swap(e2i[external_edge_id], e2i.back());
            e2i.pop_back();
            counter++;
          }
        }
        counter = 0;
      }
      std::swap(e2i_node[n.index()], e2i_node.back());
      e2i_node.pop_back();
      return 1;
    }
    return 0;
  };
  node_iterator remove_node(node_iterator n_it)
  {
    auto a = remove_node(*n_it);
    return n_it;
  };

  size_type remove_edge(const Node &n1, const Node &n2)
  {
    Node intnode1(this, e2i_node[n1.index()]);
    Node intnode2(this, e2i_node[n2.index()]);
    if (Has_edge(intnode1, intnode2))
    {
      auto edge_id = _n1_n2_edge.at(intnode1.index()).at(intnode2.index());
      _n1_n2_edge.at(intnode1.index()).erase(intnode2.index());
      _n1_n2_edge.at(intnode2.index()).erase(intnode1.index());

      auto external_edge_id = find(e2i.begin(), e2i.end(), edge_id) - e2i.begin();

      std::swap(e2i[external_edge_id], e2i.back());
      e2i.pop_back();
      return 1;
    }
    return 0;
  };
  size_type remove_edge(const Edge &edge)
  {
    return remove_edge(edge.node1(), edge.node2());
  };

  edge_iterator remove_edge(edge_iterator e_it)
  {
    remove_edge(*e_it);
    return e_it;
  };

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear()
  {
    _Node_Struct_Vector.clear();
    _n1_n2_edge.clear();
    _Edge_Struct_Vector.clear();
    e2i.clear();
    e2i_node.clear();
  }
  node_type find_external(const node_type& a)
  {
  auto external_node_id = find(e2i_node.begin(), e2i_node.end(), a.index()) - e2i_node.begin();
  Node returnnode(this,external_node_id);
  return returnnode;
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
    std::vector<size_type>::iterator _node_pointer;

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Node operator*() const
    {
      // std::cout<<(*_node_pointer).index()<<std::endl;
      return _graph_pointer->_Node_Struct_Vector[(*_node_pointer)]._Node;
      // return returnnode;
    }

    NodeIterator &operator++()
    {
      ++_node_pointer;
      return *this;
    }
    bool operator==(const NodeIterator &iter) const
    {
      if (_node_pointer == iter._node_pointer)
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
      if (_node_pointer != iter._node_pointer)
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
    NodeIterator(Graph *graph_pointer)
    {
      //const_cast<Graph *>
      _graph_pointer = graph_pointer;
      _node_pointer = _graph_pointer->e2i_node.begin();
    }
    Graph *_graph_pointer;
    size_type _id;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const
  {
    NodeIterator newiter(const_cast<Graph *>(this));
    return newiter;
  }
  node_iterator node_end() const
  {
    NodeIterator enditer(const_cast<Graph *>(this));
    enditer._node_pointer = enditer._graph_pointer->e2i_node.end();
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
    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const
    bool operator==(const IncidentIterator &iter) const
    {
      if ((_map_iter == iter._map_iter))
      {
        return true;
      }
      return false;
    }

    bool operator!=(const IncidentIterator &iter) const
    {
      if ((_map_iter != iter._map_iter))
      {
        return true;
      }
      return false;
    }
    Edge operator*() const
    {
      Edge tempedge(_graph_pointer, _map_iter->second, _node_id, _map_iter->first);
      return tempedge;
    }

    IncidentIterator &operator++()
    {
      ++_map_iter;
      return *this;
    }

  private:
    friend class Graph;
    IncidentIterator(Graph *graph_pointer, size_type node_id, bool end_or_begin)
    {
      _graph_pointer = graph_pointer;
      _node_id = node_id;

      if (graph_pointer->node(_node_id).degree() != 0)
      {
        if (end_or_begin == true)
        {
          _map_iter = _graph_pointer->_n1_n2_edge.at(_node_id).begin();
        }
        else
        {
          _map_iter = _graph_pointer->_n1_n2_edge.at(_node_id).end();
        }
      }
    }
    Graph *_graph_pointer;
    size_type _node_id;
    std::map<size_type, size_type>::iterator _map_iter;
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

    Edge operator*() const
    {
      // auto external_edge_id = find(e2i.begin(), e2i.end(), *_edge_it) - e2i.begin();

      return _graph_pointer->_Edge_Struct_Vector[(*_edge_it)]._Edge;
    }

    EdgeIterator &operator++()
    {
      ++_edge_it;
      return *this;
    }
    bool operator==(const EdgeIterator &iter) const
    {
      if (_edge_it == iter._edge_it)
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
      if (_edge_it != iter._edge_it)
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
    EdgeIterator(Graph *graph_pointer)
    {
      _graph_pointer = graph_pointer;
      _edge_it = _graph_pointer->e2i.begin();
    }
    std::vector<size_type>::iterator _edge_it;
    Graph *_graph_pointer;
  };

  edge_iterator edge_begin() const
  {
    EdgeIterator newiter(const_cast<Graph *>(this));
    return newiter;
  }
  edge_iterator edge_end() const
  {
    EdgeIterator enditer(const_cast<Graph *>(this));
    enditer._edge_it = enditer._graph_pointer->e2i.end();
    return enditer;
  }

private:
  std::vector<NodeContents> _Node_Struct_Vector;
  std::map<size_type, std::map<size_type, size_type>> _n1_n2_edge;
  std::vector<EdgeContents> _Edge_Struct_Vector;
  std::vector<size_type> e2i;
  std::vector<size_type> e2i_node;
  bool Has_edge(const Node &a, const Node &b) const // Internal Use:
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
};

#endif // CME212_GRAPH_HPP
