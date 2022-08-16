#include <iostream>
#include <math.h>

class Matrix final {
public:
  Matrix();
  Matrix(size_t cols);
  Matrix(size_t rows, size_t cols);
  ~Matrix();

  Matrix(const Matrix& mat);

  Matrix operator*(const Matrix& mat) const;
  Matrix operator-(const Matrix& mat) const;
  Matrix operator+(const Matrix& mat) const;

  Matrix operator*(double value) const;
  Matrix operator/(double value) const;

  Matrix& operator=(const Matrix& mat);
  Matrix& operator*=(const Matrix& mat);
  Matrix& operator+=(const Matrix& mat);
  Matrix& operator-=(const Matrix& mat);

  Matrix& operator*=(double value);
  Matrix& operator/=(double value);

  bool isValid() const;

  void resize(size_t rows, size_t cols);

  const double& coeffRef(size_t rowIdx, size_t colIdx) const;
  double& coeffRef(size_t rowIdx, size_t colIdx);

  const double* data() const;
  double* data();

  const size_t rows() const;
  const size_t cols() const;

  Matrix& setIdentity();
  Matrix& setZero();
  Matrix& setConstants(double value);

  Matrix& setIdentity(size_t rows, size_t cols);
  Matrix& setZero(size_t rows, size_t cols);
  Matrix& setConstants(size_t rows, size_t cols, double value);

  Matrix diag() const;
  Matrix transpose() const;
  Matrix inverse() const;
  double det() const;

  static Matrix identity(size_t rows, size_t cols);
  static Matrix zeros(size_t rows, size_t cols);
  static Matrix constants(size_t rows, size_t cols, double value);

  friend Matrix operator*(double value, const Matrix& mat);
private:
  size_t mrows_size;
  size_t mcols_size;
  double* matrix;
  bool valid_status;
};