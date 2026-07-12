#ifndef INCLUDE_INPUT_H_
#define INCLUDE_INPUT_H_

#include <iostream>
#include <fstream>
#include <string>

/// Input parametrisation

/// This code is about to be replaced with .json format 
///
namespace QSCAT
{
namespace parameters{ 
/// Class containing the parameres of the grid and a procedure for reading them from a file. 
    template<typename T>
    class grid 
    {
    public:
        bool init;
        int tnel;
        T theta;
        int nq;
        
        int nel[3];
        int elmt[3];
        
        Vector<T> * aa;
    public:
        grid()
        {
            aa = NULL;
        }
        ~grid()
        {
            if (aa) delete aa;  
            aa = NULL;
        }
        grid(const grid<T>& old):
            init(old.init),
            tnel(old.tnel),
            theta(old.theta),
            nq(old.nq)
        {
            
            nel[0] = old.nel[0]; nel[1] = old.nel[1]; nel[2] = old.nel[2];
            elmt[0] = old.elmt[0]; elmt[1] = old.elmt[1]; elmt[2] = old.elmt[2];
            
            aa = new Vector<T>(*(old.aa));
        }
        grid & operator= (grid<T> tmp)
        {
            std::swap(init,tmp.init);
            std::swap(tnel, tmp.tnel);
            std::swap(theta, tmp.theta);
            std::swap(nq, tmp.nq);
            
            std::swap(nel[0], tmp.nel[0]);
            std::swap(nel[1], tmp.nel[1]);
            std::swap(nel[2], tmp.nel[2]);
            std::swap(elmt[0], tmp.elmt[0]);
            std::swap(elmt[1], tmp.elmt[1]);
            std::swap(elmt[2], tmp.elmt[2]);
            
            std::swap(aa, tmp.aa);
            return *this;
        }
        void read_grid(std::ifstream & file)
        {
        // Auxiliary variables
            char ch[2];
            int * lc_n, *r_n, *rc_n; 
            T * lc_x, *r_x, *rc_x;
            tnel = 0;
        // Reading the parameters
            file>>elmt[0]>>ch>>elmt[1]>>ch>>elmt[2];
            file.ignore(1000,'\n');
        // Left complex scaled part of the grid
            if (elmt[0]!=0) {
                lc_n = new int[elmt[0]];
                lc_x = new T[elmt[0]];
                for (int j=0;j<elmt[0];++j){
                    file>>lc_n[j]>>ch;
                    tnel += lc_n[j];
                    if (ch[0]!=','){ break; }
                }
                file.ignore(1000,'\n');
                for (int j=0;j<elmt[0];++j){
                    file>>lc_x[j]>>ch;
                    if (ch[0]!=','){ break; }
                }
                file.ignore(1000,'\n');
            } else {
                lc_n = new int[1];
                *lc_n = 0;
                lc_x = new T[1];
                *lc_x = 0.0;
                skipline(file,2);
            }
        // Center real region of the grid
            if (elmt[1]!=0) {
                r_n = new int[elmt[1]];
                r_x = new T[elmt[1]];
                for (int j=0;j<elmt[1];++j){
                    file>>r_n[j]>>ch;
                    tnel += r_n[j];
                    if (ch[0]!=','){ break; }
                }
                file.ignore(1000,'\n');
                for (int j=0;j<elmt[1];++j){
                    file>>r_x[j]>>ch;
                    if (ch[0]!=','){ break; }
                }
                file.ignore(1000,'\n');
            } else {
                r_n = new int[1];
                *r_n = 0;
                skipline(file,2);
            }
        // Right complex scaled part of the grid
            if (elmt[2]!=0) {
                rc_n = new int[elmt[2]];
                rc_x = new T[elmt[2]];
                for (int j=0;j<elmt[2];++j){
                    file>>rc_n[j]>>ch;
                    tnel += rc_n[j];
                    if (ch[0]!=','){ break; }
                }
                file.ignore(1000,'\n');
                for (int j=0;j<elmt[2];++j){
                    file>>rc_x[j]>>ch;
                    if (ch[0]!=','){ break; }
                }
                file.ignore(1000,'\n');
            } else {
                rc_n = new int[1];
                *rc_n = 0;
                rc_x = new T[1];
                *rc_x = 0.0;
                skipline(file,2);
            }
        // Generating the ordered lengths vector for the element lengths
            int pos = 1;
            aa = new Vector<T>(tnel+1);
            T cut = 0;
            (*aa)[0] = 0.0;
            nel[0] = 0;
            for (int i=0; i<elmt[0]; ++i){
                nel[0] += lc_n[i];
                for (int k=0; k<lc_n[i]; ++k){
                    (*aa)[pos] = (lc_x[i]-cut)/lc_n[i]; // The desired vector stores only the length of the current elements. The "cut" variable subtracts the previous length
                    ++pos;
                }
                cut = lc_x[i];
            }
            nel[1] = 0;
            for (int i=0; i<elmt[1]; ++i){
                nel[1] += r_n[i];
                for (int k=0; k<r_n[i]; ++k){
                    (*aa)[pos] = (r_x[i]-cut)/r_n[i];
                    ++pos;
                }
                cut = r_x[i];
            }
            nel[2] = 0;
            for (int i=0; i<elmt[2]; ++i){
                nel[2] += rc_n[i];
                for (int k=0; k<rc_n[i]; ++k){
                    (*aa)[pos] = (rc_x[i]-cut)/rc_n[i];
                    ++pos;
                }
                cut = rc_x[i];
            }
            file>>theta;
            file.ignore(1000,'\n');
            file>>nq;
            file.ignore(1000,'\n');

            delete[] lc_n;
            delete[] r_n;
            delete[] rc_n;
            delete[] lc_x;
            delete[] r_x;
            delete[] rc_x;
            init = true;
        }
    };
/// Class containing multiple grid parametrisations
    template<typename T>
    class multi_grid
    {
    public:
        int n;                      //!< Number of defined grids
        grid<T> * gp;               //!< pointer to grid parametrs array of size n
    private:
        void ferror(const char *filename) {
            std::cout << "File '" << filename << "' is corrupted!" << std::endl;
            exit(3333);
        }
    public:
        multi_grid(){
            gp = NULL;
        }
        ~multi_grid(){
            if (gp) delete[] gp;
        }
        multi_grid(const char *filename) {
            char dumm[150]; 
            std::ifstream file(filename);
            if (file.is_open()){
                skipline(file, 4);
                file>>n>>dumm;
                gp = new grid<T>[n];
                skipline(file, 3);
                for (int i=0; i<n;++i){
                    gp[i].read_grid(file);
                    skipline(file, 2);
                }
                file.close();
            } else {
                std::cout << "Error! Could not open the grids file. " << std::endl;
                exit(44556);
            }

        }
        multi_grid(const multi_grid &old){
            n = old.n;
            gp = new grid<T>[n];
            for (int i=0;i<n; ++i){
                gp[i] = (old.gp)[i];
            }
        }
        multi_grid& operator= (multi_grid tmp){
            std::swap(n, tmp.n);
            std::swap(gp, tmp.gp);
            return *this;
        }
        bool check(const int & i){
            bool out;
            (i>=0 && i<n)? out=true : out=false;
            return out;
        }
    };

/// Testfunction parameters
    template<typename T>
    class testfunction 
    {
    public:
        char usage;         //!< Usage controller: y - used, n - unused
        char axis;          //!< The outgoing channel axis: x - electronic, y - nuclear
        int io;             //!< In/Out-going controller: (1) - in, (-1) - out
        char method;        //!< Method of S-matrix element derivation: t - tannor&weeks, d - delta-t&w, f - flux, k - t+d, l - t+f, m - d+f, a: all  
        T position;         //!< Position of the center of the testfunction or of the integral surface
        T impulse;          //!< Mean impulse of the testfunction (T&W only)
        T sigma;            //!< Wave packet mean width (T&W only)
        int channels;       //!< Number of channels, i.e. transversal states to be determined
    public:
        testfunction(){}
        void read(std::ifstream & file) {
            file>>usage;
            skipline(file,1);
            file>>axis;
            skipline(file,1);
            file>>io;
            skipline(file,1);
            file>>method;
            skipline(file,1);
            file>>position;
            skipline(file,1);
            file>>impulse;
            skipline(file,1);
            file>>sigma;
            skipline(file,1);
            file>>channels;
            skipline(file,2);
        }
        testfunction(testfunction<T> &old):
            usage(old.usage),
            axis(old.axis),
            io(old.io),
            method(old.method),
            position(old.position),
            impulse(old.impulse),
            sigma(old.sigma),
            channels(old.channels){}
        testfunction& operator= (testfunction<T> tmp){
            std::swap(usage,tmp.usage);
            std::swap(axis,tmp.axis);
            std::swap(io,tmp.io);
            std::swap(method,tmp.method);
            std::swap(position,tmp.position);
            std::swap(impulse,tmp.impulse);
            std::swap(sigma,tmp.sigma);
            std::swap(channels,tmp.channels);
            return *this;
        }
    };
/// Initial state 
    template<typename T>
    class initial_state 
    {
    public: 
        char axis;          //!< The incident coordinate: x - electronic, y - nuclear
        T position;         //!< Initial wave packet center position
        T impulse;          //!< Mean impulse of the wave packet
        T sigma;            //!< Wave packet mean width
        int channel;        //!< The incident channel, i.e. transversal state
    public: 
        void read(std::ifstream & file) {
            file>>axis;
            skipline(file,1);
            file>>position;
            skipline(file,1);
            file>>sigma;
            skipline(file,1);
            file>>impulse;
            skipline(file,1);
            file>>channel;
            skipline(file,2);
        }
        initial_state(){}
        initial_state(initial_state<T>& old):
            axis(old.axis),
            position(old.position),
            impulse(old.impulse),
            sigma(old.sigma),
            channel(old.channel){}
        initial_state& operator= (initial_state<T> tmp){
            std::swap(axis,tmp.axis);
            std::swap(position,tmp.position);
            std::swap(impulse,tmp.impulse);
            std::swap(sigma,tmp.sigma);
            std::swap(channel,tmp.channel);
            return *this;
        }
    };
/// Evolution
    template<typename T>
    class evolution 
    {
    public:
        T dt;               //!< Evolution time step (equidistant discretisation of the continuous variable
        char evolution_c;       //!< Evolution type switch (n: Crank-Nicolson, c: Chebyshev)
        int pade;           //!< Order of the pade appoximation   
        int cheb;           //!< Order of the Chebyshev polynomial
        int Q;              //!< Integration method order
        T ncutoff;          //!< The cutoff value of the normalisation
        T tcutoff;          //!< The cutoff value of the evolution time
        int steps;          //!< Number of the enegry discretisation steps
        T e_max;            //!< Upper bound of the energy
        T e_min;            //!< Lower bound of the energy
        int loop;           //!< Evolution loop steps
    public:
        evolution(){}
        void read(std::ifstream & file){
            file>>dt;
            skipline(file,1);
            file>>evolution_c;
            skipline(file,1);
            file>>pade;
            skipline(file,1);
            file>>cheb;
            skipline(file,1);
            file>>Q;
            skipline(file,1);
            file>>ncutoff;
            skipline(file,1);
            file>>tcutoff;
            skipline(file,1);
            file>>steps;
            skipline(file,1);
            file>>e_min;
            skipline(file,1);
            file>>e_max;
            skipline(file,1);
            file>>loop;
            skipline(file,2);
        }
        evolution(evolution<T>& old):
            dt(old.dt),
            evolution_c(old.evolution_c),
            pade(old.pade),
            cheb(old.cheb),
            Q(old.Q),
            ncutoff(old.ncutoff),
            tcutoff(old.tcutoff),
            steps(old.steps),
            e_max(old.e_max),
            e_min(old.e_min),
            loop(old.loop){}
        evolution& operator= (evolution<T> tmp){
            std::swap(dt,tmp.dt),
            std::swap(evolution_c,tmp.evolution_c),
            std::swap(pade,tmp.pade),
            std::swap(cheb,tmp.cheb),
            std::swap(Q,tmp.Q),
            std::swap(ncutoff,tmp.ncutoff),
            std::swap(tcutoff,tmp.tcutoff),
            std::swap(steps,tmp.steps),
            std::swap(e_max,tmp.e_max),
            std::swap(e_min,tmp.e_min),
            std::swap(loop,tmp.loop);
            return *this;
        }
    };

/// LCP Approximation parameters
    template<typename T>
    class LCP 
    {
        bool init;
    public:
        int nuclear_grid;           //!< The nuclear grid index (zero based)
        int nel_grids;              //!< The total of electronic grids
        int * electronic_grids;     //!< Indices of the electronic grids (zero based)

        int ve_channels;            //!< The total of vibrational excitation channels to be computed
        int da_channels;            //!< Dissociative attachment controller: 0 -false, 1 -true
        char method;                //!< The evolution operator approximation method
        int cn_order;               //!< Crank-Nicolson approximation order
        int cheb_order;             //!< Chebyshev approximation order

        int order;                  //!< Time quadrature order
        int steps;                  //!< Evolution loop steps
        int e_steps;                //!< Energy discretization total

        T e_min;                    //!< Energy range lower bound
        T e_max;                    //!< Energy range upper bound
    private:
        void read(const char *filename){
            std::ifstream file(filename);
            char ch[2];
            if (file.is_open()){
                skipline(file, 5);      // Skips the head
                file>>nuclear_grid;
                nuclear_grid--;         // Sets the proper value in the zero based indexing
                skipline(file, 1);
                file>>nel_grids;
                skipline(file, 1);
                electronic_grids = new int[nel_grids];
                for (int i=0; i<nel_grids; ++i){
                    file>>electronic_grids[i]>>ch;
                    electronic_grids[i]--;      // Sets the proper value in the zero based indexing
                    if (ch[0]!=','){ break; }
                }

                skipline(file, 3);
                file>>ve_channels;
                skipline(file, 1);
                file>>da_channels;
                skipline(file, 1);
                file>>method;
                skipline(file, 1);
                file>>cn_order;
                skipline(file, 1);
                file>>cheb_order;

                skipline(file, 3);
                file>>order;
                skipline(file, 1);
                file>>steps;
                skipline(file, 1);
                file>>e_steps;
                skipline(file, 1);
                file>>e_min;
                skipline(file, 1);
                file>>e_max;
            } else {
                std::cout << "Error! Could not open the 2D model parameters file." << std::endl;
                exit(44556);
            }
        }
    public:
        LCP(){
            electronic_grids = NULL;
        }
        ~LCP(){
            if (electronic_grids) delete[] electronic_grids;
        }
        LCP(const char *filename){
            read(filename); 
        }
        LCP(LCP& old):
            nuclear_grid(old.nuclear_grid), 
            nel_grids(old.nel_grids),
            ve_channels(old.ve_channels),
            da_channels(old.da_channels),
            method(old.method),
            cn_order(old.cn_order),
            cheb_order(old.cheb_order),
            order(old.order),
            steps(old.steps),
            e_steps(old.e_steps),
            e_min(old.e_min),
            e_max(old.e_max){
            
            electronic_grids = new int[nel_grids];
            for (int i=0; i<nel_grids; ++i){
                electronic_grids[i] = old.electronic_grids[i];
            }
        }
        LCP& operator= (LCP tmp){
            std::swap(nuclear_grid,tmp.nuclear_grid);   
            std::swap(nel_grids,tmp.nel_grids);
            std::swap(ve_channels,tmp.ve_channels);
            std::swap(da_channels,tmp.da_channels);
            std::swap(method,tmp.method);
            std::swap(cn_order,tmp.cn_order);
            std::swap(cheb_order,tmp.cheb_order);
            std::swap(order,tmp.order);
            std::swap(steps,tmp.steps);
            std::swap(e_steps,tmp.e_steps);
            std::swap(e_min,tmp.e_min);
            std::swap(e_max,tmp.e_max);
            std::swap(electronic_grids,tmp.electronic_grids);
        }
    };
/// Nonlocal resonance model
    template<typename T>
    class NRM
    {
    public:
        int nuclear_grid;           //!< The nuclear grid index (zero based)
        int electronic_grid;        //!< Indices of the electronic grids (zero based)

        int ve_channels;            //!< The total of vibrational excitation channels to be computed
        int da_channels;            //!< Dissociative attachment controller: 0 -false, 1 -true
        char method;                //!< The evolution operator approximation method
        int cn_order;               //!< Crank-Nicolson approximation order
        int cheb_order;             //!< Chebyshev approximation order

        int order;                  //!< Time quadrature order
        int steps;                  //!< Evolution loop steps
        int e_steps;                //!< Energy discretization total

        T e_min;                    //!< Energy range lower bound
        T e_max;                    //!< Energy range upper bound

    private:
        void read(const char *filename){
            std::ifstream file(filename);
            //char ch[2];
            if (file.is_open()){
                skipline(file, 5);      // Skips the head
                file>>nuclear_grid;
                nuclear_grid--;         // Sets the proper value in the zero based indexing
                skipline(file, 1);
                file>>electronic_grid;
                electronic_grid--;      // Sets the proper value in the zero based indexing
                
                skipline(file, 3);
                file>>ve_channels;
                skipline(file, 1);
                file>>da_channels;
                skipline(file, 1);
                file>>method;
                skipline(file, 1);
                file>>cn_order;
                skipline(file, 1);
                file>>cheb_order;

                skipline(file, 3);
                file>>order;
                skipline(file, 1);
                file>>steps;
                skipline(file, 1);
                file>>e_steps;
                skipline(file, 1);
                file>>e_min;
                skipline(file, 1);
                file>>e_max;
            } else {
                std::cout << "Error! Could not open the 2D model parameters file." << std::endl;
                exit(44556);
            }
        }
    public:
        NRM(const char *filename){
            read(filename);
        }
        NRM(){}
        NRM(NRM& old):
            nuclear_grid(old.nuclear_grid),
            electronic_grid(old.electronic_grid),       
            ve_channels(old.ve_channels),
            da_channels(old.da_channels),
            method(old.method),
            cn_order(old.cn_order),
            cheb_order(old.cheb_order),
            order(old.order),
            steps(old.steps),
            e_steps(old.e_steps),
            e_min(old.e_min),
            e_max(old.e_max){}
        NRM& operator= (NRM tmp){
            std::swap(nuclear_grid,tmp.nuclear_grid);
            std::swap(electronic_grid,tmp.electronic_grid);     
            std::swap(ve_channels,tmp.ve_channels);
            std::swap(da_channels,tmp.da_channels);
            std::swap(method,tmp.method);
            std::swap(cn_order,tmp.cn_order);
            std::swap(cheb_order,tmp.cheb_order);
            std::swap(order,tmp.order);
            std::swap(steps,tmp.steps);
            std::swap(e_steps,tmp.e_steps);
            std::swap(e_min,tmp.e_min);
            std::swap(e_max,tmp.e_max);
            return *this;
        }
    };
/// 2D model Main parameters class
    template<typename T>
    class model_2D 
    {
    public:
    // Model parameters
        bool init;
        std::string folder;
        T mu;       //!< Reduced mass
        int l;      //!< Partial wave
        T alpha_0;
        T alpha_c;
        T D_0;
        T lambda_1;
        T lambda_c;
        T lambda_inf;
        T R_0;
        T R_c;          //!< approximate crossing point of potential curves
        T R_lambda;
        
        T xRange[2];    //!< Saving range in the electronic coordinate
        T yRange[2];    //!< Saving range in the nuclear coordinate
    
        int ngx;        //!< The number of the X_coordinate (electronic) grid
        int ngy;        //!< The number of the Y_coordinate (nuclear) grid

        int xSamples;   //!< Saving samples in the electronic coordinate
        int ySamples;   //!< Saving samples in the nuclear coordinate

        testfunction<T> test_par[2];    //!< Test functions -> see definitions above 
        initial_state<T> init_par;  //!< Initial state -> see definitions above
        evolution<T> evol_par;      //!< Evolution -> see definitions above 
        
    public:
        model_2D(){}
        model_2D(const char *filename) 
        {
            read(filename);
            init=true;
        }
        model_2D(model_2D& old):
            init(old.init),
            folder(old.folder),
            mu(old.mu),
            l(old.l),
            alpha_0(old.alpha_0),
            alpha_c(old.alpha_c),
            D_0(old.D_0),
            lambda_1(old.lambda_1),
            lambda_c(old.lambda_c),
            lambda_inf(old.lambda_inf),
            R_0(old.R_0),
            R_c(old.R_c),
            R_lambda(old.R_lambda),
            ngx(old.ngx),
            ngy(old.ngy),
            init_par(old.init_par),
            evol_par(old.evol_par),
            xSamples(old.xSamples),
            ySamples(old.ySamples)
        {
            
            xRange[0] = old.xRange[0];
            xRange[1] = old.xRange[1];
            yRange[0] = old.yRange[0];
            yRange[1] = old.yRange[1];

            test_par[0] = old.test_par[0];
            test_par[1] = old.test_par[1];

        }
        model_2D& operator= (model_2D tmp)
        {
            std::swap(init,tmp.init);
            std::swap(folder,tmp.folder);
            std::swap(mu,tmp.mu);
            std::swap(l,tmp.l);
            std::swap(alpha_0,tmp.alpha_0);
            std::swap(alpha_c,tmp.alpha_c);
            std::swap(D_0,tmp.D_0);
            std::swap(lambda_1,tmp.lambda_1);
            std::swap(lambda_c, tmp.lambda_c);
            std::swap(lambda_inf,tmp.lambda_inf);
            std::swap(R_0,tmp.R_0);
            std::swap(R_c,tmp.R_c);
            std::swap(R_lambda,tmp.R_lambda);
            std::swap(ngx,tmp.ngx);
            std::swap(ngy,tmp.ngy);
            std::swap(init_par,tmp.init_par);
            std::swap(evol_par,tmp.evol_par);
            std::swap(xSamples, tmp.xSamples);
            std::swap(ySamples, tmp.ySamples);
            
            std::swap(xRange[0], tmp.xRange[0]);
            std::swap(xRange[1], tmp.xRange[1]);
            std::swap(yRange[0], tmp.yRange[0]);
            std::swap(yRange[1], tmp.yRange[1]);
            
            std::swap(test_par[0], tmp.test_par[0]);
            std::swap(test_par[1], tmp.test_par[1]);
            return *this;
        }
        void read(const char *filename) 
        {
            std::ifstream file(filename);
            if (file.is_open()){
                char d[2];
                skipline(file, 5);
                file>>folder;
                skipline(file,1);
                file>>mu>>d>>l>>d>>D_0>>d>>alpha_0>>d>>R_0>>d>>lambda_inf>>d>>lambda_1>>d>>R_lambda>>d>>lambda_c>>d>>R_c>>d>>alpha_c;
                skipline(file,2);
                file>>ngx>>d>>ngy;
                skipline(file,4);
                test_par[0].read(file);
                test_par[1].read(file);
                skipline(file,1);
                init_par.read(file);
                skipline(file,2);
                evol_par.read(file);
                skipline(file,2);
                file >> xSamples >> d >> ySamples;
                skipline(file,1);
                file >> xRange[0] >> d >> xRange[1];
                skipline(file,1);
                file >> yRange[0] >> d >> yRange[1];
                file.close();
                
            } else {
                std::cout << "Error! Could not open the 2D model parameters file." << std::endl;
                exit(44556);
            }

        }
    };
}
} // namspace QSCAT
#endif
