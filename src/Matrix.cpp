#include "Matrix.hpp"
#include <algorithm>
namespace
{
    void printMatrix(const Matrix& m)
    {
        for(int rows = 0; rows < m.rows(); rows++)
        {
            for(int cols = 0; cols < m.cols(); cols++)
            {
                std::cout << m.coeffRef(rows, cols) << " "; 
            }
            std::cout << "\n";
        }
    }
    const int PRECISION = 99;
    const double EPSILON = 0.000001;

    double _add(double fmatrix, double smatrix)
    {
        return fmatrix + smatrix;
    }

    double _dif(double fmatrix, double smatrix)
    {
        return fmatrix - smatrix;
    }

    Matrix _matrix_math(
        size_t mrows_size, 
        size_t mcols_size, 
        double* matrix,
        const Matrix& mat,
        double (*operation)(double fmatrix, double smatrix))
    {
        Matrix res(mrows_size, mcols_size);
        std::transform(matrix, matrix + mrows_size*mcols_size, mat.data(), res.data(), operation);
        return res;
    }

    void _check(const Matrix& mat, size_t rowIdx, size_t colIdx)
    {
        if(mat.rows() < rowIdx+1 || mat.cols() < colIdx+1)
        {
            throw std::out_of_range("out of range");
        }
    }

    double* mcopy(const Matrix& mat)
    {
        if(mat.isValid() == true)
        {
            double* matrix = new double[mat.rows()*mat.cols()];
            std::copy(mat.data(), mat.data() + mat.rows()*mat.cols(), matrix);
            return matrix;
        }
        return NULL;
    }
    bool operator_check(const Matrix& obj, const Matrix& mat)
    {
        return !(obj.isValid() == false || mat.cols() != obj.cols() || mat.rows() != obj.rows() || mat.isValid() == false); 
    }
}

Matrix::Matrix()
    :mcols_size(0), mrows_size(0), matrix(NULL), valid_status(false){}

Matrix::Matrix(size_t cols)
    :mrows_size(1), mcols_size(cols), valid_status(true)
{
    if(cols <= 0)
    {
        valid_status = false;
    }
    matrix = new double[cols];
}

Matrix::Matrix(size_t rows, size_t cols)
    :mcols_size(cols), mrows_size(rows), valid_status(true)
{
    if(cols <= 0 || rows <= 0)
    {
        valid_status = false;
    }
    matrix = new double[rows*cols];
}

Matrix& Matrix::operator=(const Matrix& mat)
{
    (*this).resize(mat.rows(), mat.cols());
    if(valid_status != false && mat.isValid() != false) 
    {
        double* buf = mcopy(mat);
        delete[] matrix;
        matrix = buf;
    }
    else
    {
        valid_status = false;
    }
    return *this;
}

Matrix::~Matrix()
{
    delete[] matrix;
}

Matrix::Matrix(const Matrix& mat)
: mcols_size(mat.cols()), mrows_size(mat.rows()), valid_status(mat.isValid())
{
    matrix = NULL;
    if(valid_status != false && mat.isValid() != false) 
    {
        matrix = mcopy(mat);
    }
    else
    {
        valid_status = false;
    }
}

const double& Matrix::coeffRef(size_t rowIdx, size_t colIdx) const
{
    _check(*this, rowIdx, colIdx);
    return matrix[colIdx+mcols_size*(rowIdx)];
}

double& Matrix::coeffRef(size_t rowIdx, size_t colIdx)
{
    _check(*this, rowIdx, colIdx);
    return matrix[colIdx+mcols_size*(rowIdx)];
}

Matrix Matrix::operator*(const Matrix& mat) const
{
    if(mcols_size != mat.rows() || valid_status == false || mat.isValid() == false)
    {
        return Matrix();
    }
    Matrix mul;
    mul.resize((*this).rows(), mat.cols());
    for (size_t rows = 0; rows < mrows_size; rows++)
    {
        for(size_t scols = 0; scols < mat.cols(); scols++)
        {
            mul.coeffRef(rows, scols) = 0;
            for(size_t cols = 0; cols < mcols_size; cols++)
            {
                mul.coeffRef(rows, scols) += coeffRef(rows, cols) * mat.coeffRef(cols, scols);
            }
        }
    }
    return mul;
}


Matrix Matrix::operator+(const Matrix& mat) const
{
    if (!operator_check(*this, mat))
    {
        return Matrix();
    }
    return _matrix_math(mrows_size, mcols_size, matrix, mat, _add);
}

Matrix Matrix::operator-(const Matrix& mat) const
{
    if (!operator_check(*this, mat))
    {
        return Matrix();
    }
    return _matrix_math(mrows_size, mcols_size, matrix, mat, _dif);
}

Matrix Matrix::operator*(double value) const
{
    if(valid_status == false)
    {
        return Matrix();
    }
    Matrix mul = Matrix(mrows_size, mcols_size);
    for(int cur = 0; cur < mcols_size*mrows_size; cur++)
    {
        mul.data()[cur] = matrix[cur] * value;
    }
    return mul;
}

Matrix Matrix::operator/(const double value) const
{
    if(value == 0 || valid_status == false)
    {
        return Matrix();
    }
    Matrix del(mrows_size, mcols_size);
    for(int cur = 0; cur < mcols_size*mrows_size; cur++)
    {
        del.data()[cur] = matrix[cur] / value;
    }
    return del;
}

Matrix& Matrix::operator*=(const Matrix& mat)
{
    if(mcols_size != mat.rows() || valid_status == false || mat.isValid() == false)
    {
        valid_status = false;
        return *this;
    }
    Matrix mul(mat.rows(), mat.cols());
    for (size_t rows = 0; rows < mrows_size; rows++)
    {
        for(size_t scols = 0; scols < mat.cols(); scols++)
        {
            mul.coeffRef(rows, scols) = 0;
            for(size_t cols = 0; cols < mcols_size; cols++)
            {
                mul.coeffRef(rows, scols) += coeffRef(rows, cols) * mat.coeffRef(cols, scols);
            }
        }
    }
    std::copy(mul.data(), mul.data() + mul.rows()*mul.cols(), matrix);
    return *this;
}

Matrix& Matrix::operator*=(double value)
{
    if(valid_status == false)
    {
        return *this;
    }
    for(int cur = 0; cur < mcols_size*mrows_size; cur++)
    {
        matrix[cur] *= value;
    }
    return *this;
}

Matrix& Matrix::operator/=(double value)
{
    if(value == 0 || valid_status == false)
    {
        valid_status = false;
        return *this;
    }
    Matrix mul = Matrix(mrows_size, mcols_size);
    for(int cur = 0; cur < mcols_size*mrows_size; cur++)
    {
        matrix[cur] /= value;
    }
    return *this;
}

Matrix& Matrix::operator+=(const Matrix& mat)
{
    if (!operator_check(*this, mat))
    {
        valid_status = false;
        return *this;
    }
    std::transform(matrix, matrix + mrows_size*mcols_size, mat.data(), matrix, _add);
    return *this;
}

Matrix& Matrix::operator-=(const Matrix& mat)
{
    if (!operator_check(*this, mat))
    {
        valid_status = false;
        return *this;
    }
    std::transform(matrix, matrix + mrows_size*mcols_size, mat.data(), matrix, _dif);
    return *this;
}

const size_t Matrix::cols() const
{
    return mcols_size;
}

const size_t Matrix::rows() const
{
    return mrows_size;
}

const double* Matrix::data() const
{
    return matrix;
}

double* Matrix::data()
{
    return matrix;
}

bool Matrix::isValid() const
{
    return valid_status;
}

void Matrix::resize(size_t rows, size_t cols)
{
    if(rows > 0 && cols > 0 && mcols_size == 0 && mrows_size == 0) 
    {
        valid_status = true;
    }
    if(rows == 0 || cols == 0)
    {
        valid_status = false;
        return;
    }
    double* newMatrix = new double[rows*cols];
    if(cols*rows < mcols_size*mrows_size)
    {
        std::copy(matrix, matrix + cols*rows, newMatrix);
    }
    else
    {
        std::copy(matrix, matrix + mcols_size*mrows_size, newMatrix);
    }
    delete[] matrix;
    matrix = newMatrix;
    mrows_size = rows;
    mcols_size = cols;
    return;
}

Matrix& Matrix::setZero()
{
    if (valid_status == false)
    {
        return *this;
    }
    std::fill(matrix, matrix + mcols_size*mrows_size, 0);
    return *this;
}

Matrix& Matrix::setIdentity()
{
    if(valid_status == false)
    {
        return *this;
    }
    setZero();
    for(size_t cur = 0; cur < std::min(mrows_size, mcols_size); cur++)
    {
        coeffRef(cur, cur) = 1;
    }
    return *this;
}

Matrix& Matrix::setConstants(double value)
{
    if (valid_status == false)
    {
        return *this;
    }
    std::fill(matrix, matrix + mcols_size*mrows_size, value);
    return *this;
}

static void swaprows(Matrix& matrix, int first_row, int second_row)
{
    double* row = new double[matrix.rows()];
    std::copy(matrix.data()+matrix.cols()*first_row, matrix.data()+matrix.cols()*(first_row+1), row);
    std::copy(matrix.data()+matrix.cols()*second_row, matrix.data()+matrix.cols()*(second_row+1), matrix.data()+matrix.cols()*first_row);
    std::copy(row, row+matrix.cols(), matrix.data()+matrix.cols()*second_row);
    delete[] row;
    return;
}

static void difrows(Matrix& matrix, int minuend, int subtrahend)
{
    double divider = matrix.coeffRef(subtrahend, subtrahend) / matrix.coeffRef(minuend, subtrahend);
    for(int cursor = 0; cursor < matrix.cols(); cursor++)
    {
        matrix.coeffRef(minuend, cursor) -= matrix.coeffRef(subtrahend, cursor)/divider;
    }
    return;
}

//метод Гаусса
double Matrix::det() const
{
    const double epsilon = 0.0001;
    if (mrows_size != mcols_size || valid_status == false) 
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    Matrix triangle(*this);
    double det = 1;
    for (int cols = 0; cols < mcols_size; cols++)
    {   
        int cursor = cols;
        while(fabs(triangle.coeffRef(cursor, cols)) < epsilon)
        {
            cursor++;
            if(cursor == mcols_size)
            {
                return 0;
            }
        }
        if (cursor != cols)
        {
            swaprows(triangle, cols, cursor);
            det *= -1;
        }
        for(int scols = cols+1; scols < mcols_size; scols++)
        {
            difrows(triangle, scols, cols);
        }
        //printMatrix(triangle);
        det *= triangle.coeffRef(cols,cols);
    }
    return det;
}

Matrix Matrix::transpose() const
{
    if(valid_status == false)
    {
        return Matrix();
    }
    Matrix tpos(mcols_size, mrows_size);
    for(size_t rows = 0; rows < mrows_size; rows++)
    {
        for(size_t cols = 0; cols < mcols_size; cols++)
        {
            tpos.coeffRef(cols, rows) = coeffRef(rows, cols);
        }
    }
    return tpos;
}

Matrix Matrix::inverse() const
{
    if(det() <= EPSILON && det() >= -EPSILON || valid_status == false)
    {
        return Matrix();
    }
    Matrix alAdd(mrows_size, mcols_size);
    alAdd.setZero();
    for(size_t cols = 0; cols < mcols_size; cols++)
    {
        for(size_t rows = 0; rows < mrows_size; rows++)
        {
            Matrix minor(mrows_size-1, mrows_size-1);
            size_t scur = 0;
            for(size_t scols = 0; scols < mcols_size; scols++)
            {
                if(scols != cols)
                {
                    size_t rcur = 0;
                    for(size_t srows = 0; srows < mrows_size; srows++)
                    {
                        if(srows!=rows)
                        {
                            minor.coeffRef(rcur, scur) = coeffRef(srows, scols);
                            rcur++;
                        }
                    }
                    scur++;
                }
            }
            double min = 0;
            if((rows == 0 && cols % 2 == 0) || (cols == 0 && rows % 2 == 0) || (cols % 2 == rows % 2) || (cols == rows))
            {
                min = minor.det();
                alAdd.coeffRef(rows, cols) += minor.det();
            }
            else
            {
                min = minor.det();
                alAdd.coeffRef(rows, cols) -= minor.det();
            }
        }
    }
    double dets = det();
    return (alAdd.transpose()) / det();
}
namespace
{
    Matrix _proj_u_m(const Matrix& u, const Matrix& m)
    {
        double buf = (u.transpose()*u).coeffRef(0,0);
        if(buf == 0)
        {
            return Matrix();
        }
        return u*((u.transpose()*m).coeffRef(0,0)/buf);
    }

    Matrix gram_schmidt(Matrix m)
    {
        if (m.cols() == 0) return m;
        Matrix* ort = new Matrix[m.cols()];
        for(size_t cols = 0; cols < m.cols(); cols++)
        {
            Matrix col(m.rows(), 1);
            for(size_t rows = 0; rows < m.rows(); rows++)
            {
                col.coeffRef(rows, 0) = m.coeffRef(rows, cols);
            }
            ort[cols] = col;
            if(cols>0)
            {
                for(size_t cur = 0; cur < cols; cur++)
                {
                    ort[cols] -= _proj_u_m(ort[cur], col);
                    if(ort[cols].isValid() == false)
                    {
                        delete[] ort;
                        return Matrix();
                    }
                }
            }
            double buf = ((ort[cols]).transpose() * ort[cols]).coeffRef(0,0);
            if(buf == 0)
            {
                delete[] ort;
                return Matrix();
            }
            ort[cols] = ort[cols] / sqrt(buf);
        }
        Matrix out(m.rows(), m.cols());
        for(size_t cols = 0; cols < m.cols(); cols++)
        {
            for(size_t rows = 0; rows < m.rows(); rows++)
            {
                out.coeffRef(rows, cols) = (ort[cols]).coeffRef(rows, 0);
            }
        }
        delete[] ort;
        return out;
    }
}

Matrix Matrix::diag() const
{
    //Алгоритм QR разложения Веры Кублановской
    if (mrows_size != mcols_size || valid_status == false) 
    {
        return Matrix();
    }
    Matrix A(*this);
    if(mrows_size == 1)
    {
        return A;
    }
    if(mrows_size == 2)
    {
        if(coeffRef(1,0) > coeffRef(0,0))
        {
            double buf[2] = {coeffRef(0,0), coeffRef(0,1)};
            A.coeffRef(0,0) = -coeffRef(1,0);
            A.coeffRef(0,1) = -coeffRef(1,1);
            A.coeffRef(1,0) = buf[0];
            A.coeffRef(1,1) = buf[1];
        }
        if(A.coeffRef(0,0) == 0)
        {
            return Matrix();
        }
        A.coeffRef(1,1) -= A.coeffRef(0,1) * A.coeffRef(1,0) / A.coeffRef(0,0);
        A.coeffRef(1,0) = 0;
        return A;
    }
    for(int i = 0; i < PRECISION; i++)
    {
        Matrix Q = gram_schmidt(A);
        if(Q.isValid() == false)
        {
            return Matrix();
        }
        Matrix R = Q.transpose() * A;
        A = R * Q;
    }
    for(size_t cols = 0; cols < A.cols(); cols++)
    {
        for(size_t rows = 0; rows < A.rows(); rows++)
        {
            if(cols != rows)
            {
                A.coeffRef(rows, cols) = 0;
            }
        }
    }
    return A;
}

Matrix Matrix::zeros(size_t rows, size_t cols)
{
    if (rows == 0 || cols == 0)
    {
        return Matrix();
    }
    Matrix out(rows, cols);
    std::fill(out.data(), out.data() + out.cols()*out.rows(), 0);
    return out;
}

Matrix Matrix::constants(size_t rows, size_t cols, double value)
{
    if (rows == 0 || cols == 0)
    {
        return Matrix();
    }
    Matrix out(rows, cols);
    std::fill(out.data(), out.data() + out.cols()*out.rows(), value);
    return out;
}

Matrix Matrix::identity(size_t rows, size_t cols)
{
    if (rows == 0 || cols == 0)
    {
        return Matrix();
    }
    Matrix out = zeros(rows, cols);
    for(size_t cur = 0; cur < std::min(rows, cols); cur++)
    {
        out.coeffRef(cur, cur) = 1;
    }
    return out;
}

Matrix operator*(double value, const Matrix& mat)
{
    if(mat.isValid() == false)
    {
        return Matrix();
    }
    Matrix mul(mat.rows(), mat.cols());
    for(int cur = 0; cur < mat.cols()*mat.rows(); cur++)
    {
        mul.data()[cur] = mat.data()[cur] * value;
    }
    return mul;
}

