#include "Matrix.hpp"
#include <gtest/gtest.h>
#include "valguard-memleak.h"

const double epsilon = 0.01;
Matrix matrixSetter(size_t rows_size, size_t cols_size, double* data)
{
    Matrix m(rows_size, cols_size);
    int cur = 0;
    for(size_t rows = 0; rows < rows_size; rows++)
    {
        for(size_t cols = 0; cols < cols_size; cols++)
        {
            m.coeffRef(rows, cols) = data[cur];
            cur++;
        }
    }
    return m;
}

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

TEST(Base, coeftest)
{
    Matrix a(2,2);
    double data[4] = {1,1,10,1};
    a = matrixSetter(2,2, data);
    ASSERT_TRUE(a.coeffRef(1,0) == 10);
    std::cout << valgrind_leaks().str();
}

TEST(Base, invalid)
{
    Matrix a(0,0);
    ASSERT_TRUE(a.isValid() == false);
    std::cout << valgrind_leaks().str();
}

TEST(Base, valid)
{
    Matrix a(2,2);
    ASSERT_TRUE(a.isValid() == true);
    std::cout << valgrind_leaks().str();
}

TEST(Base, copy)
{
    Matrix ma(2,2);
    double data[4] = {1,2,3,4};
    ma = matrixSetter(2,2, data);
    Matrix mb(ma);
    ASSERT_TRUE(mb.coeffRef(0,0) == 1);
    ASSERT_TRUE(mb.coeffRef(1,0) == 3);
    Matrix mc(3,4);
    mc = mb;
    ASSERT_TRUE(mc.coeffRef(0,0) == 1);
    ASSERT_TRUE(mc.coeffRef(1,0) == 3);
    Matrix md(1,1);
    md = mc;
    ASSERT_TRUE(md.coeffRef(0,0) == 1);
    Matrix ms;
    ms = ma;
    ASSERT_TRUE(ms.coeffRef(0,0) == 1);
    ASSERT_TRUE(ms.coeffRef(1,0) == 3);
    std::cout << valgrind_leaks().str();
}

TEST(Base, copyInvalid)
{
    Matrix ma;
    Matrix mb;
    mb = ma;
    ASSERT_TRUE(mb.isValid() == false);
    std::cout << valgrind_leaks().str();
}

TEST(Base, ThrowTest)
{
    //use google test try castch
    //coeffref tests
    Matrix ma(2,2);
    try
    {
        ma.coeffRef(8,8) = 5;
    }
    catch(const std::out_of_range& error)
    {
        std::cout << "out of range positive" << std::endl;
    }
    try
    {
        ma.coeffRef(-2,-3) = 5;
    }
    catch(const std::out_of_range& error)
    {
        std::cout << "out of range negative" << std::endl;
    }
    std::cout << valgrind_leaks().str();
}

TEST(ArifmeticSimple, SumAndIdentity)
{
    Matrix sum1(2,2);
    Matrix sum2(2,2);
    sum1 = sum1.setIdentity();
    sum2 = sum2.setIdentity();
    Matrix sum3(2,2);
    sum3 = sum1 + sum2;
    ASSERT_TRUE(sum3.coeffRef(0,0) == 2);
    ASSERT_TRUE(sum3.coeffRef(1,1) == 2);
    std::cout << valgrind_leaks().str();
}

TEST(ArifmeticSimple, SumAndIdentityInvalid)
{
    Matrix sumInv1(1,2);
    Matrix sumInv2(2,1);
    sumInv2.coeffRef(1,0) = 5;
    Matrix sumInv3 = sumInv1 + sumInv2;
    ASSERT_TRUE(sumInv3.isValid() == false);
    std::cout << valgrind_leaks().str();
}

TEST(ArifmeticSimple, def)
{
    Matrix sum1(2,2);
    Matrix sum2(2,2);
    sum1 = sum1.setIdentity();
    sum2 = sum2.setIdentity();
    Matrix sum3(2,2);
    sum3 = sum1 - sum2;
    ASSERT_TRUE(sum3.coeffRef(0,0) == 0);
    std::cout << valgrind_leaks().str();
}

TEST(ArifmeticSimple, defInvalid)
{
    Matrix sum1;
    Matrix sum2;
    Matrix sum3;
    sum3 = sum1 - sum2;
    ASSERT_TRUE(sum3.isValid() == false);
    std::cout << valgrind_leaks().str();
}

TEST(ArifmeticSimple, mul)
{
    double mul1_data[4] = {5, 7, 9, 2};
    Matrix mul1 = matrixSetter(2,2, mul1_data);
    double mul2_data[4] = {9, 4, 1, 3};
    Matrix mul2 = matrixSetter(2,2, mul2_data);;
    Matrix mul3 = mul1 * mul2;
    ASSERT_TRUE(mul3.coeffRef(0,0) == 52);
    ASSERT_TRUE(mul3.coeffRef(0,1) == 41);
    ASSERT_TRUE(mul3.coeffRef(1,0) == 83);
    ASSERT_TRUE(mul3.coeffRef(1,1) == 42);
    mul1.resize(0,0);
    Matrix mul4 = mul1 * mul2;
    ASSERT_TRUE(mul4.isValid() == false);
    std::cout << valgrind_leaks().str();
}

TEST(ArifmeticSimple, mulNonSquare)
{
    double mul1_data[12] = {1,0,0, 0,1,0, 0,0,1, 1,1,1};
    Matrix mul1 = matrixSetter(4,3, mul1_data);
    printMatrix(mul1);
    double mul2_data[12] = {1,0,0,0, 0,1,0,0, 0,0,1,1};
    Matrix mul2 = matrixSetter(3,4, mul2_data);
    printMatrix(mul2);
    Matrix m = mul1*mul2;
    printMatrix(m);
    ASSERT_TRUE(m.coeffRef(0,0) == 1);
}

TEST(ArifmeticSimple, mulInvalid)
{
    Matrix mul1;
    Matrix mul2;
    Matrix mul3;
    mul3 = mul1 * mul2;
    ASSERT_TRUE(mul3.isValid() == false);
    std::cout << valgrind_leaks().str();
}

TEST(ArifmeticSimple, mulValue)
{
    double mul1_data[4] = {5, 7, 9, 2};
    Matrix mul1 = matrixSetter(2,2, mul1_data);
    Matrix mul2 = 5 * mul1;
    ASSERT_TRUE(mul2.coeffRef(0,0) == 25);
    ASSERT_TRUE(mul2.coeffRef(0,1) == 35);
    ASSERT_TRUE(mul2.coeffRef(1,0) == 45);
    ASSERT_TRUE(mul2.coeffRef(1,1) == 10);
    std::cout << valgrind_leaks().str();
}

TEST(ArifmeticSimple, mulValueInvalid)
{
    Matrix mul1;
    Matrix mul2 = 5 * mul1;
    ASSERT_TRUE(mul2.isValid() == false);
    std::cout << valgrind_leaks().str();
}

TEST(ArifmeticSimple, divValue)
{
    double mul1_data[4] = {5, 10, 15, 20};
    Matrix mul1 = matrixSetter(2,2, mul1_data);
    Matrix mul2 = mul1 / 5;
    ASSERT_TRUE(mul2.coeffRef(0,0) == 1);
    ASSERT_TRUE(mul2.coeffRef(0,1) == 2);
    ASSERT_TRUE(mul2.coeffRef(1,0) == 3);
    ASSERT_TRUE(mul2.coeffRef(1,1) == 4);
    mul2 = mul1 / 0;
    ASSERT_TRUE(mul2.isValid() == false);
    std::cout << valgrind_leaks().str();
}

TEST(ArifmeticSimple, divValueInvalid)
{
    Matrix mul1;
    Matrix mul2 = mul1 / 5;
    ASSERT_TRUE(mul2.isValid() == false);
    std::cout << valgrind_leaks().str();
}

TEST(ArifmeticSimple, mulEq)
{
    double mul1_data[4] = {5, 7, 9, 2};
    Matrix mul1 = matrixSetter(2,2, mul1_data);
    mul1 *= mul1;
    ASSERT_TRUE(mul1.coeffRef(0,0) == 88);
    ASSERT_TRUE(mul1.coeffRef(0,1) == 49);
    ASSERT_TRUE(mul1.coeffRef(1,0) == 63);
    ASSERT_TRUE(mul1.coeffRef(1,1) == 67);
    std::cout << valgrind_leaks().str();
}

TEST(ArifmeticSimple, mulEqInvalid)
{
    Matrix mul1;
    mul1 *= mul1;
    ASSERT_TRUE(mul1.isValid() == false);
    std::cout << valgrind_leaks().str();
}

TEST(ArifmeticSimple, sumEq)
{
    double sum_data[4] = {5, 7, 9, 2};
    Matrix sum = matrixSetter(2,2, sum_data);
    sum += sum;
    ASSERT_TRUE(sum.coeffRef(0,0) == 10);
    ASSERT_TRUE(sum.coeffRef(0,1) == 14);
    ASSERT_TRUE(sum.coeffRef(1,0) == 18);
    ASSERT_TRUE(sum.coeffRef(1,1) == 4);
    std::cout << valgrind_leaks().str();
}

TEST(ArifmeticSimple, sumEqInvalid)
{
    Matrix sum1;
    sum1 += sum1;
    ASSERT_TRUE(sum1.isValid() == false);
    std::cout << valgrind_leaks().str();
}

TEST(ArifmeticSimple, defEq)
{
    double sum_data[4] = {5, 7, 9, 2};
    Matrix sum = matrixSetter(2,2, sum_data);
    sum -= sum;
    ASSERT_TRUE(sum.coeffRef(0,0) == 0);
    ASSERT_TRUE(sum.coeffRef(0,1) == 0);
    ASSERT_TRUE(sum.coeffRef(1,0) == 0);
    ASSERT_TRUE(sum.coeffRef(1,1) == 0);
    std::cout << valgrind_leaks().str();
}

TEST(ArifmeticSimple, defEqInvalid)
{
    Matrix def1;
    def1 -= def1;
    ASSERT_TRUE(def1.isValid() == false);
    std::cout << valgrind_leaks().str();
}

TEST(ArifmeticSimple, mulEqMatrix)
{
    double mul1_data[4] = {5, 7, 9, 2};
    Matrix mul1 = matrixSetter(2,2, mul1_data);
    mul1 *= 5;
    ASSERT_TRUE(mul1.coeffRef(0,0) == 25);
    ASSERT_TRUE(mul1.coeffRef(0,1) == 35);
    ASSERT_TRUE(mul1.coeffRef(1,0) == 45);
    ASSERT_TRUE(mul1.coeffRef(1,1) == 10);
    std::cout << valgrind_leaks().str();
}

TEST(ArifmeticSimple, mulEqMatrixInvalid)
{
    Matrix mul1;
    mul1 *= 5;
    ASSERT_TRUE(mul1.isValid() == false);
    std::cout << valgrind_leaks().str();
}

TEST(ArifmeticSimple, divEqMatrix)
{
    double mul1_data[4] = {5, 10, 15, 20};
    Matrix mul1 = matrixSetter(2,2, mul1_data);
    mul1 /= 5;
    ASSERT_TRUE(mul1.coeffRef(0,0) == 1);
    ASSERT_TRUE(mul1.coeffRef(0,1) == 2);
    ASSERT_TRUE(mul1.coeffRef(1,0) == 3);
    ASSERT_TRUE(mul1.coeffRef(1,1) == 4);
    std::cout << valgrind_leaks().str();
}

TEST(ArifmeticSimple, divEqMatrixInvalid)
{
    Matrix mul1;
    mul1 /= 5;
    ASSERT_TRUE(mul1.isValid() == false);
    std::cout << valgrind_leaks().str();
}

TEST(Algebra, det)
{
    Matrix eone(1,1);
    eone.coeffRef(0,0) = 1;
    ASSERT_TRUE(eone.det() == 1);
    //google test comp double ASSERT_DOUBLE_EQ 
    //const double epsilon = 0.000001;
    double mat_data[4] = {5, 7, 9, 2};
    Matrix mat = matrixSetter(2,2, mat_data);
    ASSERT_TRUE(mat.det() < -53 + epsilon);
    ASSERT_TRUE(mat.det() > -53 - epsilon);
    //ASSERT_TRUE(mat.det() <= -53 + epsilon && mat.det() >= -53 - epsilon);
    double e_data[9] = {1,0,0,0,1,0,0,0,1};
    Matrix e = matrixSetter(3,3, e_data);
    double det = e.det();
    ASSERT_DOUBLE_EQ(det, 1);
    //ASSERT_TRUE(det >= 1 - epsilon && det <= 1 + epsilon);
    double e4_data[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
    Matrix e4 = matrixSetter(4,4, e4_data);
    det = e4.det();
    ASSERT_DOUBLE_EQ(det, 1);
    //ASSERT_TRUE(det >= 1 - epsilon && det <= 1 + epsilon);
    double rand_data[9] = {2,3,3,6,6,4,7,2,5};
    Matrix rand = matrixSetter(3,3, rand_data);
    det = rand.det();
    ASSERT_TRUE(det < -52 + epsilon);
    ASSERT_TRUE(det > -52 - epsilon);

    double rand_dataa[16] = {1,1,1,1, 1,0,1,1, 1,1,0,1, 1,1,1,0};
    rand = matrixSetter(4,4, rand_dataa);
    det = rand.det();
    ASSERT_TRUE(det < -1 + epsilon);
    ASSERT_TRUE(det > -1 - epsilon);
    //ASSERT_TRUE(det >= -52 - epsilon && det <= -52 + epsilon);
    double data[4] = {0,1,1,0};
    Matrix four = matrixSetter(2,2, data);
    std::cout << valgrind_leaks().str();
}

TEST(Algebra, detInvalid)
{
    Matrix mat2(2,1);
    double d = mat2.det();
    ASSERT_TRUE(d!=d); // проверка на NaN или inf
    std::cout << valgrind_leaks().str();
}
TEST(Algebra, transpose)
{
    double mat_data[4] = {55, 6, 11, 23};
    Matrix mat = matrixSetter(2,2, mat_data);
    Matrix tmat = mat.transpose();
    for(size_t row = 0; row < mat.rows(); row++)
    {
        for(size_t cur = 0; cur < mat.cols(); cur++)
        {
            ASSERT_EQ(mat.coeffRef(row, cur), tmat.coeffRef(cur, row));
        }
    }
    std::cout << valgrind_leaks().str();
}

TEST(Algebra, transposeInvalid)
{
    Matrix m;
    Matrix s = m.transpose();
    ASSERT_TRUE(s.isValid() == false);
    std::cout << valgrind_leaks().str();
}

TEST(Algebra, diag)
{
    double mat_data[4] = {7, 3, 0, 5};
    Matrix mat = matrixSetter(2,2, mat_data);
    Matrix d(2,2);
    d = mat.diag();
    ASSERT_LE(d.coeffRef(1,1), 6);
    ASSERT_GT(d.coeffRef(1,1), 4);
    ASSERT_LE(d.coeffRef(0,0), 8);
    ASSERT_GT(d.coeffRef(0,0), 6);
    std::cout << valgrind_leaks().str();
}

TEST(Algebra, diagInvalid)
{
    Matrix d;
    Matrix s;
    s = d.diag();
    ASSERT_TRUE(s.isValid() == false);
    std::cout << valgrind_leaks().str();
    Matrix g(2,2);
    g *= Matrix();
    s = g.diag();
    Matrix t(3,4);
    ASSERT_TRUE(s.isValid() == false);
    s = t.diag();
    ASSERT_TRUE(s.isValid() == false);
    t *= Matrix();
    s = t.diag();
    ASSERT_TRUE(s.isValid() == false);
    Matrix x(14,8);
    x = x.diag();
    ASSERT_TRUE(x.isValid() == false);
    Matrix inval(3,3);
    inval.setZero();
    inval = inval.diag();
    printMatrix(inval);
}

TEST(Algebra, invert)
{
    double mat_data[4] = {1, 0, 1, 1};
    Matrix mat = matrixSetter(2,2, mat_data);
    Matrix invmat = mat.inverse();
    printMatrix(invmat);
    ASSERT_TRUE(invmat.coeffRef(0,0) < 1 + epsilon && invmat.coeffRef(0,0) > 1 - epsilon);
    ASSERT_TRUE(invmat.coeffRef(0,1) < 0 + epsilon && invmat.coeffRef(0,1) > 0 - epsilon);
    ASSERT_TRUE(invmat.coeffRef(1,0) < -1 + epsilon && invmat.coeffRef(1,0) > -1 - epsilon);
    double newmat_data[9] = {0,1,1, 1,0,1, 1,1,0};
    Matrix newmat = matrixSetter(3,3, newmat_data);
    printMatrix(newmat.inverse());
    double newmat_data2[4] = {1,0,0,1};
    newmat = matrixSetter(2,2, newmat_data2);
    Matrix what = newmat.inverse();
    printMatrix(what);
    //ASSERT_TRUE(invmat.coeffRef(1,1) == 1);
    double m_data[9] = {1, -2, 2, 2, 1, -1, 4, -3, 5};
    Matrix m = matrixSetter(3,3, m_data);
    Matrix inv = m.inverse();
    ASSERT_EQ(inv.coeffRef(0,0), 0.2);
    ASSERT_EQ(inv.coeffRef(0,1), 0.4);
    ASSERT_EQ(inv.coeffRef(0,2), 0);
    ASSERT_EQ(inv.coeffRef(1,0), -1.4);
    ASSERT_EQ(inv.coeffRef(1,1), -0.3);
    ASSERT_EQ(inv.coeffRef(1,2), 0.5);
    ASSERT_EQ(inv.coeffRef(2,0), -1);
    ASSERT_EQ(inv.coeffRef(2,1), -0.5);
    ASSERT_EQ(inv.coeffRef(2,2), 0.5);

    double mat_dataa[25] = {1,0,0,0,0, 0,1,0,0,0, 0,0,1,0,0, 0,0,0,1,0, 0,0,0,0,1};
    mat = matrixSetter(5,5, mat_dataa);
    invmat = mat.inverse();
    printMatrix(invmat);
    double mat_dataaa[25] = {0,1,1,1,1, 1,0,1,1,1, 1,1,0,1,1, 1,1,1,0,1, 1,1,1,1,0};
    mat = matrixSetter(5,5, mat_dataaa);
    printMatrix(mat);
    invmat = mat.inverse();
    printMatrix(invmat);
    std::cout << valgrind_leaks().str();
}

TEST(Algebra, invertInvalid)
{
    Matrix mat;
    Matrix invmat = mat.inverse();
    ASSERT_TRUE(invmat.isValid() == false);
    std::cout << valgrind_leaks().str();
}

Matrix diagonal(size_t dim, double diag = 1.0, double val = 0.0) {
Matrix matrix{ dim, dim };
for (int i = 0; i < matrix.rows(); ++i) {
for (int j = 0; j < matrix.cols(); ++j) {
matrix.coeffRef(i, j) = i == j ? diag: val;
}
}
return matrix;
} 

TEST(Algebra, Determinant) {
auto dims = { 1u,2u,3u,4u,5u,6u,11u,20u,31u,100u };
for (auto dim : dims) {
auto m = diagonal(dim, 0, 1);
double result = ((double)(dim & 1) != 0 ? 1 : -1) * ((double)dim - 1);
EXPECT_NEAR(m.det(), result, -1);
}
{
Matrix m(3, 3);
double m_data[]{
-1, 7, 4,
3, 5, -10,
2, -3, 9 };
std::copy(std::begin(m_data), std::end(m_data), m.data());
EXPECT_DOUBLE_EQ(m.det(), -420.0);
}
{
Matrix m(4, 4);
double m_data[]{
6, -3, 33, 18,
-8, 9, -2, -101,
-21, -1, 11, -8,
0, 0, -16, 1 };
std::copy(std::begin(m_data), std::end(m_data), m.data());
EXPECT_DOUBLE_EQ(m.det(), -51915.0);
}

}


TEST(Algebra, identity)
{
    Matrix e = Matrix::identity(3,3);
    ASSERT_TRUE(e.coeffRef(0,0) == 1);
    ASSERT_TRUE(e.coeffRef(1,1) == 1);
    ASSERT_TRUE(e.coeffRef(2,2) == 1);
    e = Matrix::identity(2,4);
    ASSERT_TRUE(e.coeffRef(0,0) == 1);
    ASSERT_TRUE(e.coeffRef(1,1) == 1);
    e = Matrix::identity(4,2);
    ASSERT_TRUE(e.coeffRef(0,0) == 1);
    ASSERT_TRUE(e.coeffRef(1,1) == 1);
    std::cout << valgrind_leaks().str();
}

TEST(Algebra, identityInvalid)
{
    Matrix e = Matrix::identity(0,0);
    ASSERT_TRUE(e.isValid() == false);
}

TEST(Alebra, constants)
{
    Matrix mat = Matrix::constants(3, 3, 5);
    for(size_t rows = 0; rows < mat.rows(); rows++)
    {
        for(size_t cols = 0; cols < mat.cols(); cols++)
        {
            ASSERT_EQ(mat.coeffRef(rows, cols), 5);
        }
    }
    Matrix nul = Matrix::constants(0,0, 4);
    ASSERT_TRUE(nul.isValid() == false);
    std::cout << valgrind_leaks().str();
}

TEST(Alebra, constantsInvalid)
{
    Matrix mat = Matrix::constants(0, 3, 5);
    ASSERT_TRUE(mat.isValid() == false);
    std::cout << valgrind_leaks().str();
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    #if defined(_MSC_VER)
	    testing::UnitTest::GetInstance()->listeners().Append(new testing::MSVCMemoryLeakListener());
    #endif // defined(_MSC_VER)
	return RUN_ALL_TESTS(); 
}