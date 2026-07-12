#ifndef _MOELLER_H
    #define _MOELLER_H

    /// General definition of Moeller-like oprator 

    /// template arguments:
    ///  - T .. default floating point type
    ///  - Z .. default complex floating point type
    ///  - U .. unitary evolution operator class
    ///  - H .. Hamiltonian representation class
    ///  - P .. evolved state container class

    template<typename T, class U, class H, class P>
    class MoellerOperator 
    {
        U *forward_, *backward_;
        T dt_;
        int order_;
     public:
        MoellerOperator(int order, T dt, H& h_0, H& h_f)
            : forward_(0), backward_(0), dt_(dt), order_(order)
        {
            assert(dt!=0.0);
            assert(h_0.init());
            assert(h_f.init());
          //
            forward_ = new U(order_, dt_, h_f);
            H aux(h_0);

            aux.complex_conjugate();    // necessary step : ( outgoing * -1 -> incoming )
            backward_ = new U(order_, -dt_, aux);
        }
        P project(const P& state, T range)
        {
            assert(state.init());
          //
            P out(state);
            int taskSize = range/abs(dt_); // TODO change to upper bound
           
            while(taskSize--) 
                backward_->one_step(out);
            
            taskSize = range/abs(dt_); // TODO change to upper bound
            while(taskSize--)
                forward_->one_step(out);

            return out;
        }
        ~MoellerOperator() 
        {
            if (forward_) delete forward_;
            if (backward_) delete backward_;
        }
        private:
            MoellerOperator();
    };
#endif // _MOELLER_H
