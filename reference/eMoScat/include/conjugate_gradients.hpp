#ifndef CC_GRADIENTS_H_
    #define CC_GRADIENTS_H_
    // PRECONDITIONERS SECTION

    template<typename Z, class Grid, class Operator, class State>
    class DiagonalPreconditioner
    {
        Grid grid_;
        Vector<Z> diagonal_, inv_diagonal_;
     public:
        DiagonalPreconditioner(const Grid& grid, const Operator& H) : grid_(grid), diagonal_(grid.get_size()), inv_diagonal_(grid.get_size())
        {
            Z h;
            for (blas_int i=0; i<grid_.get_size(); ++i) {
                h = H(i,i);
                diagonal_[i] = h;
                inv_diagonal_[i] = Z(1) / h;
            }
        }
        void gemv(Z alpha, const State& x, Z beta, State& y)
        {
            if (beta == 0){
                y.ax(alpha, x);
                y.body().element_wise_multiplication(diagonal_);
            } else {
                State a(grid_);
                a.ax(alpha, x);
                a.body().element_wise_multiplication(diagonal_);
                a.axpy(beta, y);
                y.swap(a);
            }
        }
        void igemv(Z alpha, const State& x, Z beta, State& y)
        {
            if (beta == 0){
                y.ax(alpha, x);
                y.body().element_wise_multiplication(inv_diagonal_);
            } else {
                State a(grid_);
                a.ax(alpha, x);
                a.body().element_wise_multiplication(inv_diagonal_);
                a.axpy(beta, y);
                y.swap(a);
            }
        }
    };

    // Template for bi-conjugate gradients method
    template<typename T, typename Z, class Operator, class State>
    State compute_biconjugate_gradients(Operator& A, State b, T prec, State x=State())
    {
        //State x(b);
        State p(b), r(b);
        State q(b), s(b);

        if (x.init()) {
            A.gemv(-1.0, x, 1.0, r);

            A.conjugate();
            State y(x);
            y.complex_conjugate();
            A.gemv(-1.0, y, 1.0, s);
            A.conjugate();

            p = r;
            q = s;
        } else {
            x = b;
            A.gemv(-1.0, x, 1.0, r);
            p=r;
            //x.fill(0.0);
        }

        // auxiliary
        State ap(x), aq(x);
        State ar(r), r0(r);
        //r0.fill(1.0);
        Z alpha=1, beta=1, delta=s*r, omega=1, rho=r*r0;

        blas_int k=0;
        while( abs(r*r) > prec && k < 100000){
            A.gemv(1.0, p, 0.0, ap);
            A.conjugate();
            A.gemv(1.0, q, 0.0, aq);
            A.conjugate();

            alpha = (s*r) / (q * ap);

            x.axpy(alpha, p);

            r.axpy(-alpha, ap);
            s.axpy(-conj(alpha), aq);

            beta = 1.0 / delta;
            delta = s * r;
            beta *= delta;

            p *= beta;
            p += r;

            q *= conj(beta);
            q += s;
          //
//            A.gemv(1.0, p, 0.0, ap);
//            alpha = rho / (ap * r0);
//
//            x.axpy(alpha, p);
//            r.axpy(-alpha, ap);
//
//            A.gemv(1.0, r, 0.0, ar);
//            omega = (ar * r) / (ar*ar);
//
//            x.axpy(omega, r);
//            r.axpy(-omega,ar);
//            
//            delta = r * r0;
//            beta = (delta/rho) * (alpha / omega);
//            rho = delta;
//
//            p *= beta;
//            p += r;
//            p.axpy(-beta*omega, ap);
          //

            k++;
            //if (k%1000==0) {
            //    cout << "iteration " << k << " : " <<  r*r << endl;
            //    cout << "TEST:\t" << ar*ar << "\ta:"<< alpha << "\tb:" << beta << "\td:" <<  delta <<"\tw:" << omega << "\ta/w:" << alpha/omega << endl;
            //}
        }

        cout << "prec: " << r*r / (x*x) << " done in k = " << k << " iterations" << endl;
        return x;
    }

    // FIXME this function is broken
    template<typename T, typename Z, class Operator, class Preconditioner, class State>
    State compute_preconditioned_biconjugate_gradients(Operator& A, Preconditioner& P, State b, T prec, State x=State())
    {

        //State x(b);
        State p(b), r(b), mr(b);
        State q(b), s(b), sm(b);

        if (x.init()) {
            A.gemv(-1.0, x, 1.0, r);

            State y(x);
            //y.complex_conjugate();
            A.conjugate();
            A.gemv(-1.0, x, 1.0, s);
            A.conjugate();

            p = r;
            q = s;

        } else {
            x = b;
            x.fill(0.0);
        }

        P.igemv(1.0, r, 0.0, p);
        P.conjugate();
        P.igemv(1.0, s, 0.0, q);
        P.conjugate();
        mr = p;
        sm = q;

        // auxiliary
        Z alpha, beta, delta=s*mr;
        State ap(x), aq(x), z(x);

        blas_int k=0;

        while( abs(r*r) > prec && k < 100000){
            A.gemv(1.0, p, 0.0, ap);
            A.conjugate();
            A.gemv(1.0, q, 0.0, aq);
            A.conjugate();

            alpha = delta / (q * ap);

            x.axpy(alpha, p);

            r.axpy(-alpha, ap);
            s.axpy(-conj(alpha), aq);

            P.igemv(1.0, r, 0.0, mr);
            P.conjugate();
            P.igemv(1.0, s, 0.0, sm);
            P.conjugate();

            beta = 1.0 / delta;
            delta = s * mr;
            beta *= delta;

            p *= beta;
            p += mr;

            q *= conj(beta);
            q += sm;

            k++;
            State f(b);
            A.gemv(-1.0, x, 1.0, f);

            if (k%1000==0)
                cout << "iteration " << k << " : " <<  r*r << endl;
        }

        cout << "prec: " << r*r / (x*x) << " done in k = " << k << " iterations" << endl;
        return x;
    }

    /// \brief Conjugate Orthoghonal Conjugate Gradients method
    /// \detail Iterative solution of complex symmetrical linear system $Ax=b$.
    template<typename T, typename Z, class Op, class State>
    State COCG(Op& A, State b, T prec, State x=State())
    {
        FILE * file;
        fopen_s(&file,"cocg_residual.txt","w");

        State p(b), r(b);

        if (!x.init())
            x = b;

        A.gemv(-1.0, x, 1.0, r);
        p=r;

        // auxiliary
        State ap(x);
        Z alpha=1, beta=1, delta=r.reduction(r);

        blas_int k=0;
        while( abs(r*r) > prec && k < 2000000){
            A.gemv(1.0, p, 0.0, ap);

            alpha = delta / p.reduction(ap);
            x.axpy(alpha, p);

            r.axpy(-alpha, ap);

            beta = 1.0 / delta;
            delta = r.reduction(r);
            beta *= delta;

            p *= beta;
            p += r;

            k++;
            fprintf(file, "%.12e\n", abs(r*r));
            if (k%100==0) {
                cout << "iteration " << k << " : " <<  r*r << endl;
                //cout << "TEST:\t" << r*r << "\ta:"<< alpha << "\tb:" << beta << "\td:" <<  delta << endl;
            }
        }

        fclose(file);
        cout << "prec: " << r*r / (x*x) << " done in k = " << k << " iterations" << endl;
        return x;

    }

    /// \brief Preconditioned Conjugate Orthoghonal Conjugate Gradients method
    /// \detail Iterative solution of complex symmetrical linear system $Ax=b$.
    template<typename T, typename Z, class Op, class Preconditioner, class State>
    State PCOCG(Op& A, Preconditioner& P, State b, T prec, State x=State())
    {
        FILE * file;
        fopen_s(&file,"pcocg_residual.txt","w");

        State p(b), r(b), mr(b), mp(b);

        if (!x.init())
            x = b;
        //P.igemv(1.0, x, 0.0, mp);
        //mp = x;
        A.gemv(-1.0, x, 1.0, r);
        //p=r;
        P.igemv(1.0, r, 0.0, mr);
        //P.gemv(1.0, r, 0.0, mr);
        //mr = r;
        p = mr;

        // auxiliary
        State q(x);
        Z beta=1, delta=r.reduction(mr);

        A.gemv(1.0, p, 0.0, q);
        Z alpha = delta / p.reduction(q);
        x.axpy(alpha, p);
        r.axpy(-alpha, q);

        blas_int k=0;
        while( abs(r*r) > prec && k < 2000000){
            beta = 1.0 / delta;
            P.igemv(1.0, r, 0.0, mr);
            delta = r.reduction(mr);     // r <=> mr

            beta *= delta;
            p *= beta;
            p += mr;     // r <=> mr

            //P.igemv(1.0, p, 0.0, mp );
            //mp = p;
            A.gemv(1.0, p, 0.0, q);

            alpha = delta / p.reduction(q);
            x.axpy(alpha, p);
            r.axpy(-alpha, q);
            //P.gemv(1.0, r, 0.0, mr);
            //mr = r;


            k++;
            fprintf(file, "%.12e\n", abs(r*r));
            fflush(file);

            if (k%100==0) {
                cout << "iteration " << k << " : " <<  r*r << endl;
                //cout << "TEST:\t" << r*r << "\ta:"<< alpha << "\tb:" << beta << "\td:" <<  delta << endl;
            }
        }

        fclose(file);
        cout << "prec: " << r*r / (x*x) << " done in k = " << k << " iterations" << endl;
        //r = x;
        //P.gemv(1.0, r, 0.0, x);
        return x;

    }
#endif
