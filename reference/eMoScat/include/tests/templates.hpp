using std::string;
using std::cout;
using std::endl;

string red("\e[0;31m");
string green("\e[0;32m");
string yellow("\e[0;33m");
string ecolor("\e[0m");

size_t ERROR_COUNT = 0;
size_t TEST_COUNT = 0;

/// This is the main testin call macro. Wraps the testing function to the output environment.

#define TEST(func, ... ) TEST_COUNT+=1; ERROR_COUNT += printer(TESTS::func( __VA_ARGS__ ) );

namespace TESTS
{

string SECTION = "none";

struct TestStatus
{
    int errors;
    string name;
    string msg;

    TestStatus(const std::string& tname) { errors = 0; name=tname; }
    TestStatus& add_msg(const std::string& message) { msg += message; return *this; }
};

int printer(const TestStatus& status)
{
    if (status.errors > 0) {
        cout << red << "[--== FAILED ==--]" << ecolor;
        cout << " " << SECTION << " : " << status.name << " : " << status.msg <<  endl;
        return 1;
    } else {
        cout << green << "[--== PASSED ==--]" << ecolor;
        cout << " " << SECTION << " : " << status.name << endl;
    }
    return 0;
}

template<class P>
TestStatus swap(const P& X, const P& Y)
{
    assert(X.init());
    assert(Y.init());
    assert(X.get_size() == Y.get_size());
  //
    TestStatus status("swap operation");
    P Z1(X);
    P Z2(Y);
    Z1.swap(Z2);

    for (int i=0; i<Y.get_size(); ++i)
        if (Z1[i] != Y[i] && Z2[i] != X[i]) status.errors++;

    if (status.errors > 0)
        return status.add_msg("failed!");

    return status.add_msg("passed");
}

template<class P, typename T>
TestStatus fill(const P& X, const T& val)
{
    assert(X.init());
  //
    TestStatus status("fill operation");
    P Z(X);
    Z.fill(val);

    for (int i=0; i<X.get_size(); ++i)
        if (Z[i] != val) status.errors++;

    if (status.errors > 0)
        return status.add_msg("failed!");

    return status.add_msg("passed");
}

template<class P>
TestStatus copy(const P& X)
{
    assert(X.init());
  //
    TestStatus status("copy/assignement");

    P Y(X);

    for (int i=0; i<Y.get_size(); ++i)
        if (Y[i] != X[i]) status.errors++;

    if (status.errors > 0)
        return status.add_msg("Copy constructor failed!");

    P Z;
    Z = X;

    for (int i=0; i<Y.get_size(); ++i)
        if (Y[i] != X[i]) status.errors++;

    if (status.errors > 0)
        return status.add_msg("Assignement operator failed!");

    Z = X.copy();

    for (int i=0; i<Y.get_size(); ++i)
        if (Y[i] != X[i]) status.errors++;

    if (status.errors > 0)
        return status.add_msg("Copy function failed!");

    return status;
}

template<class P>
TestStatus inplace_add(const P& X, const P& Y)
{
    assert(X.init());
    assert(Y.init());
    assert(X.get_size() == Y.get_size());
  //
    TestStatus status("inplace addition");
    P Z(Y);
    Z += X;

    for (int i=0; i<Y.get_size(); ++i)
        if (Z[i] != X[i] + Y[i]) status.errors++;

    if (status.errors > 0)
        return status.add_msg("failed!");

    return status.add_msg("passed");
}

template<class P>
TestStatus inplace_sub(const P& X, const P& Y)
{
    assert(X.init());
    assert(Y.init());
    assert(X.get_size() == Y.get_size());
  //
    TestStatus status("inplace subtraction");
    // inplace subtraction
    P Z(Y);
    Z -= X;

    for (int i=0; i<Y.get_size(); ++i)
        if (Z[i] != Y[i] - X[i]) status.errors++;

    if (status.errors > 0)
        return status.add_msg("Inplace subtraction failed!");

    return status.add_msg("passed");
}

template<class P, typename T>
TestStatus inplace_scaling(const P& X, const T& scal)
{
    assert(X.init());
  //
    TestStatus status("inplace scaling");
    // inplace subtraction
    P Z(X);
    Z *= scal;

    for (int i=0; i<X.get_size(); ++i)
        if (Z[i] != X[i] * scal) status.errors++;

    if (status.errors > 0)
        return status.add_msg("Inplace scaling failed!");

    return status.add_msg("passed");

}

template<class H, class P, typename T>
TestStatus general_matrix_vector(T alpha, const H& A, const P& x, T beta, P& y)
{
    TestStatus status("gemv operation");

    P s = y.copy();
    A.gemv(alpha, x, beta, s);

    blas_int rows = A.rows();
    blas_int cols = A.columns();

    blas_int nnz = A.num_nonzeros();

    blas_int pos = 0;

    for (blas_int i=0; i<rows; ++i) {
        T sum = beta * y[i];
        for (blas_int n=A.row_index(i); n<A.row_index(i+1); ++n) {
            sum += alpha*x[A.columns(n)] * A.nonzeros(n);
        }
        if (abs(sum - s[i]) > 1e-12*abs(sum)){
            status.errors++;
            cout << i << ": " << sum << " <> " << s[i] << endl;
        }
    }

    if (status.errors > 0)
        return status.add_msg("failed!");

    return status.add_msg("passed");
}

} // namespace TESTS
