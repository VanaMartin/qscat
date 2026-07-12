
    boost::python::def("bessel_jn", &BesselJ);

    class_<parameters::grid<double> >("gridParameters", init<>())
        .def(init<const parameters::grid<double>&>())
    ;

    class_<parameters::multi_grid<double> >("multiGridParameters", init<const char*>())
        .def(init<const parameters::multi_grid<double>&>())
        .def("__getitem__", &multiGridItem)
    ;
