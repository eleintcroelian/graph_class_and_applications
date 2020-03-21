#include <vector>
/** This function returns the transpose of a given a valid matrix, 
   * represented by a vector<vector<T>>
   * Complexity: O(m*n)
   */
template <typename T>
std::vector<std::vector<T>> transpose(const std::vector<std::vector<T>> &mat)
{
    int n = mat.size();
    int m = mat[0].size();
    std::vector<std::vector<T>> transpose_matrix(m, std::vector<T>(n, 0));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            transpose_matrix[i][j] = mat[j][i];
};