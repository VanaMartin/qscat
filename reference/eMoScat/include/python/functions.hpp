
    // TODO from dictionary
    gVector getGaussian(const femGrid& g, def_float x0, def_float sigma, def_float p0)
    {
        def_comp core;
        gVector x(g);
        for (int i=0; i<g.nb(); ++i) {
            core = exp(-pow((g.x(i)-x0),2)/(4*pow(sigma,2)))/sqrt(sqrt(pi*2)*sigma);
            x.f(core * exp(imu*p0*g.x(i)), i);

        }
        return x;

        //pjvalue cfg;
        //std::string err = picojson::parse(cfg, src);
        //gVector x(g);
        //fill_grid_vector(x, cfg, Gaussian);
        //return x;
    }

    gVector getBesselJ(const femGrid& g, def_float p, def_float j, def_float mu)
    {
        gVector x(g);
        for (int i=0; i<g.nb(); ++i) {
            x.f(sphBesselJEn(g.x(i), p, mu, j), i);
        }
        return x;
    }

    gVector getHankel1(const femGrid& g, def_float p, def_float j, def_float mu)
    {
        gVector x(g);
        for (int i=0; i<g.nb(); ++i) {
            x.f(sphHankel1En(g.x(i), p, mu, j), i);
        }
        return x;
    }

    gVector getCoulombF(const femGrid& g, def_float p, def_float j, def_float mu, def_float q)
    {
        gVector x(g);
        for (int i=0; i<g.nb(); ++i) {
            x.f(coulomb::sF_en(g.x(i), p, q, mu, j), i);
        }
        return x;
    }

    gVector getCoulombH1(const femGrid& g, def_float p, def_float j, def_float mu, def_float q)
    {
        gVector x(g);
        for (int i=0; i<g.nb(); ++i) {
            x.f(coulomb::sH1_en(g.x(i), p, q, mu, j), i);
        }
        return x;
    }

    gVector pyMorse(const femGrid& g, const std::string& src)
    {
        pjvalue cfg;
        std::string err = picojson::parse(cfg, src);
        gVector x(g);
        fill_grid_vector(x, cfg, MorsePotential);
        return x;
    }

    gVector pyLambdaInf(const femGrid& g, const std::string& src)
    {
        pjvalue cfg;
        std::string err = picojson::parse(cfg, src);
        gVector x(g);
        fill_grid_vector(x, cfg, AsymptoticLambda);
        return x;
    }

    gVector pyElectronicLambda(const femGrid& g, dcomp y, const std::string& src)
    {
        pjvalue cfg;
        std::string err = picojson::parse(cfg, src);
        gVector x(g);
        fill_grid_vector_xaxis(x, y, cfg, ElectronicLambdaInteraction);
        return x;
    }

    // 2d

    gVector2D pyLambdaInteraction(const femGrid2D& g, const std::string &src)
    {
        pjvalue cfg;
        std::string err = picojson::parse(cfg, src);
        gVector2D x(g);
        fill_grid_vector_2d(x, cfg, LambdaInteraction);
        return x;
    }

    gVector2D py2dPotential(const femGrid2D& g, const std::string &src)
    {
        pjvalue cfg;
        std::string err = picojson::parse(cfg, src);
        gVector2D x(g);
        fill_grid_vector_2d(x, cfg, Neutral2dPotential);
        return x;
    }
