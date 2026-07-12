#include "qscat.h"
#include "potentials.h"

using namespace std;
using namespace QSCAT;

int main(int argc, char **argv) {

    pjvalue cfg;
    if (argc == 1){
        cout << "Configuration File not found, running N_2 model" << endl;
        cfg = read_json_file("input/experimental/N2-model.json");
    } else {
        cout << "Configuration File: " << argv[1] << endl;
        cfg = read_json_file(argv[1]);
    }

	LCP::ModelLCP LCP(cfg);
    gVector2D PhiD = LCP.get_discrete_state(cfg);
	TimeDependentModel2D M2D(cfg, &PhiD);


    const pjvalue ep = cfg["evolution"];
	for (int i=0; i<int(ep["time_cutoff"].asDouble()/(ep["time_step"].asDouble() * ep["loop_steps"].asDouble()));++i){
		M2D.multistep();
        LCP.multistep();
		//M2D.save_binary((m2dp.folder + "frame.M2D").c_str());
	}
	return 0;
}

