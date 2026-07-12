#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <complex>
#include <cassert>
#include <math.h>
#include <stdlib.h>

#ifdef linux
    #include <sec_stream.h>
#endif

#include "FemDvrEcs/DiscreteStates.h"

namespace QSCAT
{
    // internals
    bool DiscreteStates1D::read(const char *name, const FemDvrEcsGrid& grid)
    {
        grid_ = grid;
        std::ifstream file(name);
        char ch[2];
        bool stat;
        blas_int basis_size;
        dfloat work;
        if (file.is_open()){
          // Success - File opened
            skipline(file, 1);
            file>>ch>>num_discrete_states_;
            skipline(file, 1);
            file>>ch>>basis_size;
            skipline(file, 1);
            states_ = new GridVector[num_discrete_states_];
            spectrum_ = new dcomp[num_discrete_states_];
            for (blas_int j=0; j<num_discrete_states_; ++j){
                states_[j] = GridVector(grid);
                spectrum_[j] = 0;
                file>>ch>>work;
                spectrum_[j] += work;
                file>>ch>>work;
                spectrum_[j] += imu*work;
            }
            skipline(file, 1);
            for (blas_int i=0; i<basis_size; ++i){
                for (blas_int j=0; j<num_discrete_states_; ++j){
                    states_[j][i] = 0;
                    file>>work;
                    states_[j][i] += work;
                    file>>ch[0]>>work;
                    states_[j][i] += imu*work;
                }
                skipline(file, 1);
            }
            stat = true;
            *init_ = true;
            file.close();
        } else {
        // Fail - File not found
            stat = false;
            *init_ = false;
        }
        return stat;
    }
    void DiscreteStates1D::clean()
    {
        if (*init_ && num_discrete_states_!=0){
            delete[] spectrum_;
            delete[] states_;
        }
        *init_ = false;
        num_discrete_states_ = 0;
    }
    bool DiscreteStates1D::save_bin_body(std::ofstream & file) const
    {
        assert(init());
      //
        bool stat = file.is_open();
        if (stat) {
            stat = grid_.save_binary(file);
            file.write((char*) &num_discrete_states_, sizeof(blas_int));
            file.write((char*) &max_energy_, sizeof(dfloat));
            file.write((char*) &machine_epsilon_, sizeof(dfloat));
            file.write((char*) &precision_, sizeof(dfloat));
            file.write((char*) spectrum_, num_discrete_states_*sizeof(dcomp));
            for (blas_int i=0; i<num_discrete_states_; ++i){
                if(stat) stat = states_[i].save_binary(file);
            }
        }
        return stat;
    }
    bool DiscreteStates1D::read_bin_body(std::ifstream & file)
    {
        bool stat = file.is_open();
        if (stat){
            stat = grid_.read_binary(file);
            if (init()) {
                delete[] states_;
                delete[] spectrum_;
            }
            file.read((char*) &num_discrete_states_, sizeof(blas_int));
            file.read((char*) &max_energy_, sizeof(dfloat));
            file.read((char*) &machine_epsilon_, sizeof(dfloat));
            file.read((char*) &precision_, sizeof(dfloat));
            spectrum_ = new dcomp[num_discrete_states_];
            file.read((char*) spectrum_, num_discrete_states_*sizeof(dcomp));
            states_ = new GridVector[num_discrete_states_];
            for (blas_int i=0;i<num_discrete_states_;++i){
                if(stat) stat = states_[i].read_binary(file);
            }
            *init_ = true;
        }
        return stat;
    }

    // constructors
    DiscreteStates1D::DiscreteStates1D() : Object()
    {
        *init_ = false;
        num_discrete_states_ = 0;
        max_energy_ = 0;
        machine_epsilon_ = 0;
        precision_ = 0;
        spectrum_ = 0;
        states_ = 0;
    }
    DiscreteStates1D::DiscreteStates1D(blas_int num_grids, FemDvrEcsGrid* grids, GridVector* potentials, dfloat mass) : Object()
    {
        assert(grids);
        assert(potentials);
        assert(num_grids>1);
      //
        grid_ = grids[0];
        max_energy_ = 1.0;
        machine_epsilon_ = 1.0e-16;
        precision_ = 1.0e-4;
        num_discrete_states_ = 0;
        EigenSystem<dcomp> * eSys = new EigenSystem<dcomp>[num_grids];

      // Building the hamiltonians and its spectra
        for (blas_int i=0; i<num_grids; ++i) {
            assert(grids[i] == potentials[i].get_grid());
          //
            OperatorFull H(grids[i]);
            H.add_kinetic_term(mass);
            H += potentials[i];
            eSys[i] = H.eigen_system();
        }

    // Comparison of the eigenvalues
        blas_int grid_index[2];
        blas_int upper_bound = grids[0].nb();

        for (blas_int i=1; i<num_grids; ++i){
            if (grids[i].nb() < upper_bound) {
                upper_bound = grids[i].nb();
            }
        }

        blas_int *index = new blas_int[upper_bound];

    // Change the ordering?

        for (blas_int i=0; i<num_grids-1; ++i){
            for (blas_int j=i+1; j<num_grids; ++j){
                for (blas_int k=0; k<upper_bound-1; ++k){
                    for (blas_int l=0; l<upper_bound; ++l){
                        //if (abs(Ham[i]->eigs->energies[k] - Ham[j]->eigs->energies[l]) < precision_*abs(Ham[j]->eigs->energies[l])) {
                        if (abs(eSys[i].eigen_value(k) - eSys[j].eigen_value(l)) < precision_*abs(eSys[j].eigen_value(l))) {
                            index[num_discrete_states_] = k;
                            num_discrete_states_++;
                            grid_index[0] = i;
                            grid_index[1] = j;
                        }
                    }
                }
                if (num_discrete_states_ != 0) {
                    break;
                }
            }

            if (num_discrete_states_ != 0) {
                break;
            }
        }

        if (num_discrete_states_ != 0) {
            std::cout << "... success. " << num_discrete_states_ << " discrete states were found!" << std::endl;
            spectrum_ = new dcomp[num_discrete_states_];
            states_ = new GridVector [num_discrete_states_];
            for (blas_int i=0; i<num_discrete_states_; ++i){
                spectrum_[i] = eSys[grid_index[0]].eigen_value(index[i]);
                states_[i] = GridVector(grids[grid_index[0]]);
                eSys[grid_index[0]].eigen_vector((states_[i].body()),index[i]);
            }

        } else {
            std::cout << "... failure. No discrete states were found!" << std::endl;
        }

        delete[] eSys;
        *init_ = true;
    }
    DiscreteStates1D::DiscreteStates1D(blas_int num_grids, FemDvrEcsGrid * grids, GridVector * potentials, dfloat mass, dfloat a_prec, dfloat min_eig, dfloat max_eig) : Object()
    {
        grid_ = grids[0];
        max_energy_ = 1.0;
        machine_epsilon_ = 1.0e-16;
        precision_ = a_prec;
        num_discrete_states_ = 0;
        EigenSystem<dcomp> *eSys = new EigenSystem<dcomp>[num_grids];
        zVector *auxSpec = new zVector[num_grids];

    // Building the hamiltonians and its spectra
        for (blas_int i=0; i<num_grids; ++i){
            OperatorFull H(grids[i]);
            H.add_kinetic_term(mass);
            H += potentials[i];
            eSys[i] = H.eigen_system();
        }

    // Comparison of the eigenvalues
        blas_int grid_index[2];
        blas_int upper_bound = grids[0].nb();

        for (blas_int i=1; i<num_grids; ++i){
            if (grids[i].nb() < upper_bound) {
                upper_bound = grids[i].nb();
            }
        }

        blas_int * index = new blas_int[upper_bound];

        for (blas_int i=0; i<1; ++i){		// #MODIFIED FOR FIRST GRID ONLY -- assures the phi_d compatibility
            for (blas_int j=i+1; j<num_grids; ++j){
                for (blas_int k=0; k<upper_bound-1; ++k){
                    for (blas_int l=0; l<upper_bound; ++l){
                        //if (real(Ham[i]->eigs->energies[k]) < max_eig && real(Ham[i]->eigs->energies[k]) > min_eig) {
                        if (real(eSys[i].eigen_value(k)) < max_eig && real(eSys[i].eigen_value(k)) > min_eig) {
                            dfloat distance = abs(eSys[i].eigen_value(k) - eSys[j].eigen_value(l));
                            dfloat relative_precision = precision_*abs(eSys[i].eigen_value(k));

                            if (distance < ((relative_precision > 1e-8)? relative_precision:precision_)) {
                                if (num_discrete_states_ < upper_bound) {
                                    index[num_discrete_states_] = k;
                                    num_discrete_states_++;
                                    grid_index[0] = i;
                                    grid_index[1] = j;
                                } else {
                                    break;
                                }
                            }
                        }
                    }
                }
                if (num_discrete_states_ != 0) {
                    break;
                }
            }

            if (num_discrete_states_ != 0) {
                break;
            }
        }
    // Testing area: saving the spectrum_
//        //T test1, test2, test3;
//        zVector ** out = new zVector* [3] ;
//        out[0] = new zVector(eSys[0].get_size());
//        out[1] = new zVector(eSys[0].get_size());
//        out[2] = new zVector(eSys[0].get_size());
//        for (blas_int i=0; i<upper_bound; ++i){ // eSys[0].get_size()
//            (*(out[0]))[i] = eSys[0].eigen_value(i);
//            (*(out[1]))[i] = eSys[1].eigen_value(i);
//            (*(out[2]))[i] = eSys[2].eigen_value(i);
//        }
//        SaveMultipleVectors(3, upper_bound, out, "output/eigen.erg");  // eSys[0].get_size()
//        delete out[0];
//        delete out[1];
//        delete out[2];
//        delete out;
    // end of the testing area

        if (num_discrete_states_ != 0) {
            std::cout << "... success. " << num_discrete_states_ << " discrete state(s) were found";
            spectrum_ = new dcomp[num_discrete_states_];
            states_ = new GridVector [num_discrete_states_];
            for (blas_int i=0; i<num_discrete_states_; ++i){
                spectrum_[i] = eSys[grid_index[0]].eigen_value(index[i]);
                states_[i] = GridVector(grids[grid_index[0]]);
                eSys[grid_index[0]].eigen_vector((states_[i].body()),index[i]);
            }

        } else {
            std::cout << "... failure. No discrete state(s) were found";
        }

        delete[] auxSpec;
        delete[] index;

    // Cleaning the Hamiltonians
        //for (int i=0; i<num_grids; ++i){
        //	if (Ham[i] != num_gridsULL){
        //		delete Ham[i];
        //	}
        //}
        delete[] eSys;
        *init_ = true;
        incref();
    }
    DiscreteStates1D::DiscreteStates1D(char *name, const FemDvrEcsGrid& grid) : Object()
    {
        read(name, grid);
    }
    DiscreteStates1D::~DiscreteStates1D()
    {
        if(decref()==0)
            clean();
    }

    // accessors
    void DiscreteStates1D::retrieve(zVector*& E, zVector ** V)
    {
        // TODO reimpelment
        assert(0);
    //
    //    E = new vector<Z>(*spectrum_);
    //    for (int i=0; i<num_discrete_states_; ++i){
    //        V[i] = new GridVector(*(states_[i]));
    //    }
    }
    // TODO : change all occurences to exact copy
    void DiscreteStates1D::get_state(GridVector& destination, blas_int i)
    {
        assert(init());
        assert(i < num_discrete_states_);
      //
        destination = states_[i].copy();
    }
    const FemDvrEcsGrid& DiscreteStates1D::get_grid() const
    {
        return grid_;
    }
    dcomp DiscreteStates1D::get_energy(blas_int i)
    {
        assert(i < num_discrete_states_);
      //
        return spectrum_[i];
    }
    blas_int DiscreteStates1D::number_of_states()
    {
        assert(init());
      //
        return num_discrete_states_;
    }

    // storage
    bool DiscreteStates1D::from_file(const char *name, FemDvrEcsGrid& grid)
    {
        if (init()) { clean(); }
        return read(name,grid);
    }
    void DiscreteStates1D::save_states(const char *name)
    {
        assert(init());
      //
        FILE * file;
        fopen_s(&file,name,"w");
        blas_int nb = states_[0].get_size();
        fprintf(file,"#Multiple vectors stored in one file.\n");
        fprintf(file,"# %lld - number of states_.\n", num_discrete_states_);
        fprintf(file,"# %lld - number of basis functions of the grid.\n", nb);
        fprintf(file,"# %.12E, %.12E", real(spectrum_[0]), imag(spectrum_[0]));
        for (blas_int i=1; i<num_discrete_states_; ++i){
            fprintf(file,"#%.12E, %.12E", real(spectrum_[i]), imag(spectrum_[i]));
        }
        fprintf(file,"\n");
        for (blas_int j=0; j<nb; ++j){
            fprintf(file,"%.12E\t%.12E", real(states_[0][j]), imag(states_[0][j]));
            for (blas_int i=1; i<num_discrete_states_; ++i){
                fprintf(file,"%.12E\t%.12E", real(states_[i][j]), imag(states_[i][j]));
            }
            fprintf(file,"\n");
        }
        fclose(file);
    }
} // namespcece QSCAT
