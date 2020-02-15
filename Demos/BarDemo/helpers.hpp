#include <Eigen/Dense>
#include "Common/Common.h"
#include <assert.h>

namespace Eigen
{
// https://stackoverflow.com/questions/25389480/how-to-write-read-an-eigen-matrix-from-binary-file
template <class Matrix>
void write_binary(const std::string filename, const Matrix &matrix)
{
    std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
    typename Matrix::Index rows = matrix.rows(), cols = matrix.cols();
    out.write((char *)(&rows), sizeof(typename Matrix::Index));
    out.write((char *)(&cols), sizeof(typename Matrix::Index));
    out.write((char *)matrix.data(), rows * cols * sizeof(typename Matrix::Scalar));
    out.close();
}
template <class Matrix>
void read_binary(const std::string filename, Matrix &matrix)
{
    std::ifstream in(filename, std::ios::in | std::ios::binary);
    typename Matrix::Index rows = 0, cols = 0;
    in.read((char *)(&rows), sizeof(typename Matrix::Index));
    in.read((char *)(&cols), sizeof(typename Matrix::Index));
    matrix.resize(rows, cols);
    in.read((char *)matrix.data(), rows * cols * sizeof(typename Matrix::Scalar));
    in.close();
}
} // namespace Eigen

class TrajectoryData
{
public:
    TrajectoryData(std::vector<Eigen::Isometry3d> _transformations, Eigen::VectorXd _dt_inv)
    {
        transformations = _transformations;
        dt_inv = _dt_inv;
    }

    TrajectoryData(std::string path_prefix)
    {
        // Read trajectory from directory
        const int num_rows = 30;

        transformations.resize(num_rows);
        for (int i = 0; i < num_rows; i++)
        {
            // std::cout << i << std::endl;
            Eigen::read_binary(path_prefix + std::to_string(i) + ".dat", transformations[i].matrix());
            // std::cout<<transformations[i].translation()<<std::endl;
        }

        Eigen::read_binary(path_prefix + "dt.dat", dt_inv);
        // std::cout << dt_inv << std::endl;
        assert(dt_inv.cols() == transformations.size());
    }

    void step(Real currentTime)
    {
        idx++;
        Real dt = (idx >= dt_inv.size()) ? 1e8 : 1 / dt_inv[idx];
        nextUpdateTime += dt;

        std::cout << "STEP:" << idx << "/" << dt_inv.size() << std::endl;
        std::cout << "CUR:" << currentTime << " NEXT:" << nextUpdateTime << std::endl;
        std::cout << transformations[idx].translation() << std::endl
                  << std::endl;
    }

    bool needUpdate(Real currentTime)
    {
        return currentTime >= nextUpdateTime;
    }

    void transform(Vector3r &pos)
    {
        assert(idx < dt_inv.size());
        pos = transformations[idx] * pos;
    }

    void reset()
    {
        idx = 0;
        nextUpdateTime = 0;
    }

    std::vector<Eigen::Isometry3d> transformations;
    Eigen::VectorXd dt_inv;

private:
    int idx = 0;
    Real nextUpdateTime = 0;
};