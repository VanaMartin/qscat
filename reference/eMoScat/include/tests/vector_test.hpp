/// This file contains case specific tests for vector only class. The tests are templated so
/// all possible specifications may be tested

int vector_size = 65535;    /// Standard Vector length

template<class V>
TestStatus vector_construction()
{
    TestStatus("Default constructor");
    V(vector_size);
}


template<class V>
void general_vector_test(const std::string& secname)
{
    int errors = 0;

    cout << endl;
    V X(vector_size);
    V Y(vector_size);
    for (int i=0; i<vector_size; ++i) {
        X[i] = i;
        Y[i] = i + vector_size;
    }

    TESTS::SECTION = secname;
    TEST(copy, X);
    TEST(fill, X, 0);
    TEST(swap, X, Y);
    TEST(inplace_add, X, Y);
    TEST(inplace_sub, X, Y);

    return errors;
}
