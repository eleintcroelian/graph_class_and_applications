/**
 * @file mtl_test.cpp
 * Test script for interfacing with MTL4 and it's linear solvers.
 */

// HW3: Need to install/include Boost and MTL in Makefile
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include <cassert>
#include <iostream>

// HW3: YOUR CODE HERE
// Define a IdentityMatrix that interfaces with MTL
class IdentityMatrix
{
  /** Helper function to perform multiplication . Allows for delayed
* evaluation of results .
* Assign :: apply (a, b) resolves to an assignment operation such as
* a += b, a -= b, or a = b.
* @pre @a size (v) == size (w) */
public:
  IdentityMatrix(int m) : m_(m) {}
  template <typename VectorIn, typename VectorOut, typename Assign>
  void mult(const VectorIn &v, VectorOut &w, Assign) const
  {
    assert(size(v) == size(w));
    assert(int(size(v)) == m_);
    for (size_t i = 0; i < size(v); i++)
    {
      w[i] = v[i];
    }
  }
  /** Matvec forwards to MTL 's lazy mat_cvec_multiplier operator */
  template <typename Vector>
  mtl::vec::mat_cvec_multiplier<IdentityMatrix, Vector>
  operator*(const Vector &v) const
  {
    return mtl::vec::mat_cvec_multiplier<IdentityMatrix, Vector>(*this, v);
  }
  unsigned int m_;
};

inline std::size_t size(const IdentityMatrix &A) { return A.m_ * A.m_; }
inline std::size_t num_rows(const IdentityMatrix &A) { return A.m_; }
inline std::size_t num_cols(const IdentityMatrix &A) { return A.m_; }

/** Traits that MTL uses to determine properties of our IdentityMatrix . */
namespace mtl
{
namespace ashape
{
/** Define IdentityMatrix to be a non - scalar type . */
template <>
struct ashape_aux<IdentityMatrix>
{
  typedef nonscal type;
};
} // end namespace ashape

/** IdentityMatrix implements the Collection concept
 * with value_type and size_type */
template <>
struct Collection<IdentityMatrix>
{
  typedef double value_type;
  typedef unsigned size_type;
};
} // end namespace mtl

int main()
{
  using namespace std;
  typedef mtl::dense_vector<double> vt;

  const int size = 10, N = size * size;

  typedef IdentityMatrix matrix_type;
  matrix_type A(size);

  vt v(size);
  iota(v);
  cout << "v is " << v << endl;

  vt w2(size);

  w2 = A * v;
  cout << "A * v is " << w2 << endl;

  w2 += A * v;
  cout << "w2+= A * v is " << w2 << endl;

  w2 -= A * v;
  cout << "w2-= A * v is " << w2 << endl;


  itl::pc::identity<matrix_type> P(A);

  mtl::dense_vector<double> x(size, 1.0), b(size);

  b = A * x;
  x = 0;
  cout << "b is " << b << endl;
  cout << "x is " << x << endl;
  itl::cyclic_iteration<double> iter(b, 100, 1.e-11, 0.0, 5);
  cg(A, x, b, P, iter);
  cout << "b is " << b << endl;
  cout << "x is " << x << endl;


  return 0;
}
