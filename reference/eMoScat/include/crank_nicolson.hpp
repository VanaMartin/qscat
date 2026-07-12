# ifndef __CRANK_NICOLSON__
    # define __CRANK_NICOLSON__

namespace QSCAT
{
        /// General Crank-Nicolson expansion of exponential operator

        ///
        /// Template specifications allow combination of vectors and operators such as:
        ///      -            GridVector & OperatorFull, OperatorRowCompressed
        ///      -          GridVector2D & Operator2D
        ///      -    DoubleGridVector2D & double_operator_2D\
        ///      -          MaskVector2D & MaskOperator2D
        ///
        ///  specifications of template arguments:
        ///      - T: default floating point variable
        ///      - Z: default complex floating point variable
        ///      - H: operator class (associatied with T, Z)
        ///      - P: state vector class (associated with T, Z)
        ///
        /// Requirements:
        ///
        ///    H:
        ///        - operation *= Z                (scale with scalar)
        ///        - operation += Z                (add scalar to diagonal)
        ///        - method    gemv(P)             (matrix-vector multiplication)
        ///        - method    factorizeLU()       (LU decomposition)
        ///        - method    backSubstitution(P) (matrix-vector back substitution)
        ///
        ///    P:
        ///        - copy constructor
        ///        - method    swap(P)             (swap internal members)
        ///

        template<typename T, typename Z, class H, class P>
        class CrankNicolson
        {
            private:
                bool init;
                int order;
                T dt;
                H *numerators;
                H *denominators;
                Vector<Z> roots;
            private:
                void initialize(int inorder, T indt, H& inH)
                {
                    // assertions
                    assert(inH.init());
                    assert(indt!=T(0));
                    assert(inorder>0);
                    assert(inorder<=20);
                    // code
                    if (inorder > 19) {
                        order = 20;
                    } else if (inorder > 14) {
                        order = 15;
                    } else if (inorder > 9) {
                        order = 10;
                    } else {
                        order = inorder;
                    }
                    dt = indt;
                    roots = Vector<Z>(order);
                    Pade_Roots(roots, order);
                    // auxiliary
                    Z idtoverroot;

                    numerators = new H[order];
                    denominators = new H[order];

                    for (int i=0; i<order; ++i){
                        // numerators
                        H& p = numerators[i];
                        idtoverroot = imu * dt / roots[i];
                        p = inH.copy();
                        p *= idtoverroot;
                        p += Z(1.0);
                        assert(p.init());
                        // denominators
                        H& q = denominators[i];
                        idtoverroot = -imu * dt / roots[i];
                        q = inH.copy();
                        q *= idtoverroot;
                        q += Z(1.0);
                        q.LU_factorize();
                        assert(q.init());
                    }
                    init = true;
                }
                void clean()
                {
                    if (init) {
                       delete[] numerators;
                       delete[] denominators;
                       numerators = NULL;
                       denominators = NULL;
                       order = 0;
                       dt = 0;
                    }
                }
            public:
                CrankNicolson()
                {
                    init=false;
                    order = 0;
                    dt = 0;
                    numerators = 0;
                    denominators = 0;
                }
                CrankNicolson(int inorder, T indt, H& inH)
                {
                    initialize(inorder, indt, inH);
                }
                void set(int inorder, T indt, H& inH)
                {
                    if (init) clean();
                    initialize(inorder, indt, inH);
                }
                void clear()
                {
                    clean();
                }
                ~CrankNicolson()
                {
                    clean();
                }
                void one_step(P& x)
                {
                    // assertions
                    assert(init);
                    // auxiliary memory for gemv
                    P y = x.copy();
                    for (int i=0; i<order; ++i){
                       // numerators
                       H& p = numerators[i];
                       p.gemv(Z(1), x, Z(0), y);
                       x.swap(y);
                       H& q = denominators[i];
                       q.LU_back_substitution(x);
                    }
                }
            private:
                CrankNicolson(const CrankNicolson& old);
                CrankNicolson& operator= (const CrankNicolson& rhs);
        };
} // namespace QSCAT
# endif
