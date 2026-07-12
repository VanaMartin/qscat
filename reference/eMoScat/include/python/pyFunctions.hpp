
    def("gVectorGaussian", &getGaussian);

    def("gVectorBesselJ", &getBesselJ);
    def("gVectorHankel1", &getHankel1);
    def("gVectorCoulombF", &getCoulombF);
    def("gVectorCoulombH1", &getCoulombH1);

    def("gVectorMorse", &pyMorse);
    def("gVectorLambdaInf", &pyLambdaInf);
    def("gVectorElectronicLambda", &pyElectronicLambda);
    
    def("gVector2dLambda", &pyLambdaInteraction);
    def("gVector2dPotential", &py2dPotential);
