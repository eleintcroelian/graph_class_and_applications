#include <vector>
/** This function reshapes a given a valid matrix, represented by a vector<vector<T>>
   * with dimensions r,c. If the given dimensions are not valid, original matrix is returned.
   * Complexity: O(i*j)
   */
template <typename T>
std::vector<std::vector<T>> matrixReshape(const std::vector<std::vector<T>> &mat,
                                          std::size_t r, std::size_t c)
{
    unsigned int n = mat.size();
    unsigned int m = mat[0].size();
    if (r * c != n * m)
    {
        return mat;
    }
    std::vector<std::vector<T>> new_matrix(r, std::vector<T>(c, 0));
    unsigned int current_row = 0;
    unsigned int current_col = 0;
    for (unsigned int i = 0; i < n; i++)
    {
        for (unsigned int j = 0; j < m; j++)
        {
            new_matrix[current_row][current_col] = mat[i][j];
            current_col++;
            if (current_col == c) // if we hit the new matrix shape column number
            {
                current_row++; // increment row number and zero the column number
                current_col = 0; // meaning get to the next row
            }
        }
    }
    return new_matrix;
};
