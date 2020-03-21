#include <vector>
/** This function checks whether a given a valid matrix, represented by a vector<vector<T>>
   * is symmetric. If matrix is non-square result is false.
   * Complexity: Worst case scenario, O((m-1)*(n/2))
   */
template <typename T>
bool isSymmetric(const std::vector<std::vector<T>> &matrix)
{
    if (matrix.size()!=matrix[0].size())
    {
        return false;
    }
    for (int i = 0; i < matrix.size(); i++)
        for (int j = i+1; j < matrix[0].size(); j++)
            if (matrix[i][j] != matrix[j][i])
                return false;
    return true;
};