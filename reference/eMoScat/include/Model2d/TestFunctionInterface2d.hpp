namespace QSCAT
{
template<typename State>
TestFunctionInterface2d<State>::TestFunctionInterface2d()
{
    init_ = false;
    opened_ = false;
}

template<typename State>
TestFunctionInterface2d<State>::TestFunctionInterface2d(const TestFunctionInterface2d<State>& old):
    init_(old.init_),
    axis_(old.axis_),
    channel_(old.channel_),
    in_out_(old.in_out_),
    size_(old.size_),
    xsize_(old.xsize_),
    ysize_(old.ysize_),
    mu_x_(old.mu_x_),
    mu_y_(old.mu_y_),
    impulse_momentum_(old.impulse_momentum_),
    charge_(old.charge_),
    energy_(old.energy_),
    initial_energy_(old.initial_energy_),
    energy_shift_(old.energy_shift_),
    bound_state_(old.bound_state_),
    fourier_coefficients_(old.fourier_coefficients_),
    energies_(old.energies_),
    quad_order_(old.quad_order_),
    buffer_(old.buffer_),
    coefficients_(old.coefficients_),
    grid_(old.grid_)
{
    outfile_ = NULL;
    opened_ = false;
}

template<typename State>
TestFunctionInterface2d<State>& TestFunctionInterface2d<State>::swap(TestFunctionInterface2d<State>& rhs)
{
    std::swap(init_, rhs.init_);
    std::swap(axis_, rhs.axis_);
    std::swap(channel_, rhs.channel_);
    std::swap(in_out_, rhs.in_out_);
    std::swap(size_, rhs.size_);
    std::swap(xsize_, rhs.xsize_);
    std::swap(ysize_, rhs.ysize_);
    std::swap(mu_x_, rhs.mu_x_);
    std::swap(mu_y_, rhs.mu_y_);
    std::swap(impulse_momentum_, rhs.impulse_momentum_);
    std::swap(charge_, rhs.charge_);
    std::swap(energy_, rhs.energy_);
    std::swap(initial_energy_, rhs.initial_energy_);
    std::swap(energy_shift_, rhs.energy_shift_);
    bound_state_.swap(rhs.bound_state_);
    fourier_coefficients_.swap(rhs.fourier_coefficients_);
    energies_.swap(rhs.energies_);
    std::swap(quad_order_, rhs.quad_order_);
    //buffer_.swap(rhs.buffer_);                // NOTE buffer information lost : will be fixed by shallow copy semantics
    coefficients_.swap(rhs.coefficients_);
    std::swap(outfile_, rhs.outfile_);
    std::swap(opened_, rhs.opened_);
    grid_.swap(rhs.grid_);
    return *this;
}

template<typename State>
bool TestFunctionInterface2d<State>::init() const
{
    return init_;
}
template<typename State>
def_float TestFunctionInterface2d<State>::energy(int i) const
{
    assert(init_);
    assert(i >= 0);
    assert(i<energies_.get_size());
  //
    return energies_[i];
}
template<typename State>
const def_comp& TestFunctionInterface2d<State>::initial_energy() const
{
    assert(init_);
  //
    return initial_energy_;
}
template<typename State>
const def_float& TestFunctionInterface2d<State>::energy_shift() const
{
    assert(init_);
  //
    return energy_shift_;
}
template<typename State>
const def_float& TestFunctionInterface2d<State>::reduced_mass() const
{
    assert(init_);
  //
    return (axis_ == 'x')? mu_x_: mu_y_;
}

template<typename State>
void TestFunctionInterface2d<State>::set_output(const std::string& filename)
{
    assert(!opened_);
  //
    outfile_ = new std::ofstream;
    outfile_->open(filename.c_str(), std::ios::out | std::ios::binary);
    opened_ = true;
    return;
}

template<typename State>
TestFunctionInterface2d<State>::~TestFunctionInterface2d()
{
    if (opened_) {
        outfile_->close();
        delete outfile_;
    }
}

}; // namespace QSCAT
