# ifndef __CHEBYSHEV__
    # define __CHEBYSHEV__

namespace QSCAT
{
        /// General Chebyshev expansion of exponential operator

        ///
        /// Template specifications allow combination of vectors and operators such as:
        ///      -           grid_vector & Hamiltonian
        ///      -        grid_vector_2D & Hamiltonian 2D
        ///      - double_grid_vector_2D & double_operator_2D\
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
        ///
        ///    P:
        ///        - copy constructor
        ///        - method    swap(P)             (swap internal members)
        ///

        #include "bessel.h"

        template<typename T, typename Z, class H, class P>
        class Chebyshev
        {
            private:
                bool init_;     //!< Initialization controller
                int order_;     //!< order of the approximation
                T dt_;          //!< time step
                Vector<Z> coefficients_;     //!</ Decomposition coefficients
                H op_;          //!< generating Operator
                Z max_;         //!< Spectral maximum
                Z min_;         //!< Spectral minimum

                P s1_;          //!< Auxiliary state
                P s2_;          //!< Auxiliary state
                P s3_;          //!< Auxiliary state
            private:
                void initialize(int inorder, T indt, H& inH)
                {
                    int max_iteration = 10000;
                    op_ = inH.copy();
                    dt_ = indt;
                    order_ = inorder;

                  // First determine the spectral radius by repeated multiplication
                    // TODO set random
                    P test(inH.get_grid());
                    P test2 = test.copy();
                    for (int i=0; i<max_iteration; ++i) {
                        op_.gemv(1.0, test, 0.0, test2);
                        max_ = sqrt(test2*test2);
                        test2 *= 1.0 / max_;
                        test.swap(test2);
                    }
                    // TODO set random again
                    op_ += -max_;
                    for (int i=0; i<max_iteration; ++i) {
                        op_.gemv(1.0, test, 0.0, test2);
                        min_ = sqrt(test2*test2);
                        test2 *= 1.0 / min_;
                        test.swap(test2);
                    }
                    min_ += max_;
                    // TODO check values
                    // reset to initial state
                    // H = ( 2 * H - (max+min) Id ) / (max-min)
                    op_ = inH.copy();
                    Z diagonal = - (max_ + min_) / abs(max_ - min_);
                    Z factor   = 2.0 / abs(max_ - min_);

                    op_ *= factor;  // scale with factor
                    op_ += diagonal;

                    Z R = -imu*dt_*(max_ - min_)/2.0;
                    Z G = -imu*dt_*min_;
                    T r = abs(R);
                    T q;

                    coefficients_ = Vector<Z>(order_);
                    for (int i=0; i<order_; ++i){
                        Z j = BesselJ(r, i);
                        (i == 0)? q=1.0 : q=0.0;
                        coefficients_[i] = (pow(-imu,i)) * exp(R + G)*(2.0 - q) * j;
                    }

                    init_ = true;
                }
                void clean()
                {
                    if (init_) {
                       order_ = 0;
                       dt_ = 0;
                       min_ = 0;
                       max_ = 0;
                       init_ = false;
                    }
                }
            public:
                Chebyshev()
                {
                    init_=false;
                    order_ = 0;
                    dt_ = 0;
                    min_ = 0;
                    max_ = 0;
                }
                Chebyshev(int inorder, T indt, H& inH)
                {
                    initialize(inorder, indt, inH);
                }
                void set(int inorder, T indt, H& inH)
                {
                    if (init_) clean();
                    initialize(inorder, indt, inH);
                }
                void clear()
                {
                    clean();
                }
                ~Chebyshev()
                {
                    clean();
                }
                void one_step(P& x)
                {
                    s1_ = x.copy();       // n-th element
                    s2_ = x.copy();       // n-1 th element
                    s3_ = x.copy();       // n-2 th element
                    Z a;

                    for (int k=0;k<order_;++k){
                        a = coefficients_[k];
                        if (k==0) {         // Zeroth polynom is identity
                            x.ax(a,s1_);
                        } else if (k==1) {
                            op_.gemv(Z(1.0), s2_,Z(0.0), s1_); // Extended Matrix multiplication
                            x.axpy(a,s1_);
                            s2_ = s1_.copy();
                        } else {
                        // The recursive relation for Chebyshev polynomials: T_(n+1)(X) = 2XT_n - T_(n-1)
                            op_.gemv(Z(2.0), s1_, Z(-1.0), s3_); // 2 A_n - A_(n-2)
                            s2_.swap(s3_);        // Swap with the new value (n-1) <-> (n+1)
                            s1_ = s2_.copy();        // n+1 copied for next step
                            x.axpy(a,s1_);
                        }                       // repeat cycle
                    }
                }
        };
}   // namespace QSCAT
# endif
