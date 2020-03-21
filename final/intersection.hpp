#include <vector>
#include <algorithm>
/** This function returns the intersection of 3 vectors with type T.
   * Resultant vector is sorted in descending order.
   * Complexity: O(nlogn)
   */
template <typename T>
std::vector<T> intersection(const std::vector<T> &a, const std::vector<T> &b,
                            const std::vector<T> &c)
{

    std::vector<T> resultant_vector;
    std::vector<T> intermediate_vector;
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    std::sort(c.begin(), c.end());

    std::set_intersection(a.begin(), a.end(),
                          b.begin(), b.end(),
                          back_inserter(intermediate_vector));

    std::set_intersection(intermediate_vector.begin(), intermediate_vector.end(),
                          c.begin(), c.end(),
                          back_inserter(resultant_vector));

    std::sort(resultant_vector.begin(), resultant_vector.end(),std::greater<T>());

    return resultant_vector;
};