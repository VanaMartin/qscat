# The upstream quantum-chemistry slot is named psi4. Its SA backend reads the
# explicit PySCF MO inventory, preserving core/active/virtual Molden block order.
use JSON::PP;

sub read_state_average_output {
    my ($r_par) = @_;
    open(my $input, '<', "$ENV{'UKRMOL_RUN_DIR'}/target.json") or die $!;
    my $report = decode_json(do { local $/; <$input> });
    close($input);
    die "Unconverged state-averaged target" unless $report->{'converged'};
    my $data = $r_par->{'data'};
    my $orbs = $data->{'orbitals'};
    $orbs->{'target_all'} = [(0) x $data->{'nir'}];
    my @ordered;
    foreach my $orb (@{$report->{'orbitals'}}) {
        my $ir = 1;
        $ir++ until $irred_repr{'C2v'}->[$ir - 1] eq $orb->{'irrep'};
        my $key = "$ir.$orb->{'index_in_irrep'}";
        push(@ordered, $key);
        $orbs->{'energies'}->{$key} = $orb->{'energy_hartree'};
        $orbs->{'target_all'}->[$ir - 1]++;
    }
    $data->{'target'}->{'noccupied'} = 7;
    $data->{'target'}->{'nfull_occ'} = 7;
    $data->{'target'}->{'nhalf_occ'} = 0;
    $orbs->{'target_sum'} = scalar @ordered;
    $data->{'scf_ok'} = 1;
    &check_numbers_of_orbitals($r_par);
    &set_orbitals($r_par, \@ordered);
    return 1;
}
