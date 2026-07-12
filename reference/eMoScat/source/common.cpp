#include <iostream>
#include <fstream>
#include <complex>

#include "common.h"
// Auxiliary functions for the eMoScat 

bool functions::check_mem_size(const int elements)
{
	float mem_size;
	mem_size = 8*elements*2;
	if (mem_size > 500*pow(1024.0,2)) {
		std::cout << "Warning! The estimated array memory size " << mem_size/pow(1024.0,2) << " MB is exceeding the upper limit. Check the parametrization for lowering the size." << std::endl;
	}
	return true;
}

// Globaly defined function for skipping the rest of line & multiple others
void skipline(std::ifstream& file, const int & i)
{
	for (int j=0;j<i;++j){
		file.ignore(1000,'\n');
	}
}
