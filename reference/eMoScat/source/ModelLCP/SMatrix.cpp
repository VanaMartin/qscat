#include "ModelLCP/SMatrix.h"

using std::cout;
using std::endl;
using std::exp;

namespace QSCAT
{

void LCP::SMatrix::initialize(const pjvalue& parameters)
{
	mu_ = parameters["model"]["reduced_mass"].asDouble();
	init_channel_ = parameters["initial_state"]["channel"].asInt();
	order_ = 2; // FIXME
	channels_ = parameters["cross_sections"]["testfunctions"]["vibrational_excitation"]["channels"].asInt();
	ve_channels_ = channels_;
	da_channels_ = 0;
    if (parameters["cross_sections"]["testfunctions"].isMember("dissociative_attachment")) {
	    da_channels_ = parameters["cross_sections"]["testfunctions"]["dissociative_attachment"]["channels"].asInt();
        channels_ += da_channels_;
    }
	steps_ = parameters["evolution"]["loop_steps"].asInt();
    def_float emin = parameters["cross_sections"]["range"][0].asDouble();
    def_float emax = parameters["cross_sections"]["range"][1].asDouble();
	size_ = (emax - emin)/parameters["cross_sections"]["dE"].asDouble();
	make_coefficients();
	energies_ = dVector(size_, emin, emax, false);
	buffer_ = new zVector[channels_];
	s_ = new zVector[channels_];
	for (int i=0; i<channels_; ++i){
		buffer_[i] = zVector(steps_+1);
		buffer_[i].fill(0.0);
		s_[i] = zVector(size_);
		s_[i].fill(0.0);
	}
	init_ = true;
}

void LCP::SMatrix::clean()
{
	if (init_) {
		delete[] buffer_;
		delete[] s_;
	}
    init_ = false;
}

void LCP::SMatrix::make_coefficients()
{
	switch (order_){
		case 0:
			coefficients_ = dVector(steps_);
			coefficients_.fill(1.0);
			break;
		case 1:
			coefficients_ = dVector(steps_+1);
			coefficients_.fill(1.0);
			coefficients_[steps_-1] = 0.5;
			coefficients_[steps_] = 0.5;		// The last element has to be added again if the integration continues another time
			break;
		case 2:
			coefficients_ = dVector(steps_+1);
			coefficients_[0] = 4.0/3.0;
			for (int i=1;i<steps_-1;++i){
				if (i%2==1) {
					coefficients_[i]= 4.0/3.0;
				} else {
					coefficients_[i]= 2.0/3.0;
				}
			}
			coefficients_[steps_-1] = 1.0/3.0;
			coefficients_[steps_] = 1.0/3.0;		// The last element has to be added again if the integration continues another time
	}
}

LCP::SMatrix::SMatrix()
{
	init_=false;
}

LCP::SMatrix::SMatrix(const pjvalue& parameters)
{
	initialize(parameters);
}

LCP::SMatrix::~SMatrix()
{
	if (init_) {
		clean();		// Flush previous data
	}
}

//LCP::SMatrix& LCP::SMatrix::set(const pjvalue& parameters)
//{
//	initialize(parameters);
//	return *this;
//}

void LCP::SMatrix::contribution(gVector& psi, gVector* states, int i)
{
	for (int j=0; j<ve_channels_; ++j){
		// Assuming the states to already contain the factor sqrt(Gamma(R)/2pi)
		buffer_[j][i] = states[j]*psi;
	}
	for (int j=ve_channels_; j<channels_; ++j){
		buffer_[j][i] = psi.f(psi.get_grid().nr());
	}
}

void LCP::SMatrix::close_multistep(def_float time, def_float dt, dVector& ergs, const def_float& ierg, const def_float& X)
{
	for (int i=0; i<ve_channels_; ++i){		// Contributions to vibrational excitations
		def_float shift = ierg - ergs[i];
		for (int k=0; k<size_; ++k){
			if (energies_[k] + shift>0) {
				if (order_ > 0) { // Contribution of the previous loop end
					s_[i][k] += (1.0/imu)*exp(imu*(energies_[k] + ierg)*(time))*buffer_[i][steps_]*dt*coefficients_[steps_];
				}
				for (int j=0; j<steps_; ++j){
					s_[i][k] += (1.0/imu)*exp(imu*(energies_[k] + ierg)*(time + (j+1)*dt))*buffer_[i][j]*dt*coefficients_[j];
				}
			}
		}
		buffer_[i][steps_] = buffer_[i][0];
	}
	for (int i=ve_channels_; i<channels_; ++i){	// Contributions to dissociative attachment
		def_float shift = ierg - ergs[i];
		for (int k=0; k<size_; ++k){
			if (energies_[k] + shift>0) {
				def_float K = sqrt(2.0*mu_*(energies_[k]+shift));
				if (order_ > 0) { // Contribution of the previous loop end
					s_[i][k] += sqrt(K/(2.0*pi*mu_)) * exp(-imu*K*X) *
						exp( imu*(energies_[k] + ierg)*(time)) * buffer_[i][steps_] *
						dt * coefficients_[steps_];
				}
				for (int j=0; j<steps_; ++j){
					s_[i][k] += sqrt(K/(2.0*pi*mu_)) * exp(-imu*K*X) *
						exp( imu*(energies_[k] + ierg)*(time + (j+1)*dt)) * buffer_[i][j] *
						dt * coefficients_[j];
				}
			}
		}
		if (order_ > 0) {
			buffer_[i][steps_] = buffer_[i][0];
		}
	}
}

// Derives the cross sections and stores them into a file with appropriate time
void LCP::SMatrix::cross_sections(const def_float& time, std::string& folder)
{
	char name[50];
	def_float val;
	def_float erg;
	if (ve_channels_ != 0){
		sprintf_s(name, "/CS/CSVE_t=%.1f.dat", time);
		FILE * file;
		fopen_s(&file,(folder+name).c_str(),"w");
		fprintf(file,"#Cross Sections of vibrational excitations: channels, values\n");
		fprintf(file,"#\t%d, %d\n", ve_channels_, size_);
		for (int i=0; i<size_; ++i){
			erg = energies_[i];
			fprintf(file, "%.12E", erg);
			for (int j=0; j<ve_channels_; ++j){
				val = pow(abs(s_[j][i]),2)*4.0*pow(pi,3)/(2.0*erg);
				fprintf(file, "\t%.12E", val);
			}
			fprintf(file, "\n");
		}
		fclose(file);
	}
	if (channels_>ve_channels_) {
		sprintf_s(name, "/CS/CSDA_t=%.1f.dat", time);
		FILE * file;
		fopen_s(&file,(folder+name).c_str(),"w");
		fprintf(file,"#Cross Sections of vibrational excitations: channels, values\n");
		fprintf(file,"#\t%d, %d\n", channels_ - ve_channels_, size_);
		for (int i=0; i<size_; ++i){
			erg = energies_[i];
			fprintf(file, "%.12E", erg);
			for (int j=ve_channels_; j<channels_; ++j){
				val = pow(abs(s_[j][i]),2)*4.0*pow(pi,3)/(2.0*erg);
				fprintf(file, "\t%.12E", val);
			}
			fprintf(file, "\n");
		}
		fclose(file);
	}
}

} // namespace QSCAT
