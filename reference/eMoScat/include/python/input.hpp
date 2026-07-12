
parameters::grid<double> multiGridItem(const parameters::multi_grid<double>& object, int i)
{
    assert(i < object.n);
    return object.gp[i];
}

