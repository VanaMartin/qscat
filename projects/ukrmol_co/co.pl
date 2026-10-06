# Defaults are read from the checksum-pinned UKRmol-scripts release.
my $upstream = $ENV{'UKRMOL_SCRIPTS'};
do "$upstream/scripts/model.pl" or die "Cannot load model defaults: $@ $!";
do "$upstream/scripts/config.pl" or die "Cannot load run defaults: $@ $!";

%dirs = (
    'bin_in' => "/opt/ukrmolp/bin.$ENV{'UKRMOL_PRECISION'}",
    'bin_out' => "/opt/ukrmolp/bin.$ENV{'UKRMOL_PRECISION'}",
    'psi4' => '/opt/ukrmolp/bin',
    'molpro' => '', 'molcas' => '',
    'basis' => "$upstream/basis.sets",
    'templates' => "$ENV{'UKRMOL_RUN_DIR'}/templates",
    'libs' => "$ENV{'UKRMOL_RUN_DIR'}/lib",
    'output' => "$ENV{'UKRMOL_RUN_DIR'}/output",
);

$model{'directory'} = 'CO';
$model{'molecule'} = 'CO';
$model{'atoms'} = ['C', 'O'];
$model{'nelectrons'} = 14;
$model{'symmetry'} = 'C2v';
$model{'basis'} = $ENV{'UKRMOL_BASIS'};
$model{'orbitals'} = $ENV{'UKRMOL_ORBITALS'};
$model{'select_orb_by'} = 'energy';
$model{'model'} = $ENV{'UKRMOL_MODEL'};
my $nfrozen = $ENV{'UKRMOL_FROZEN_ORBITALS'};
$model{'nfrozen'} = $nfrozen;
$model{'frozen_orbs'} = [$nfrozen,0,0,0,0,0,0,0];
my @active = split /,/, $ENV{'UKRMOL_ACTIVE_ORBITALS'};
$model{'active_orbs'} = [@active,0,0,0,0];
$model{'nactive'} = 0;
$model{'nactive'} += $_ foreach @active;
# Explicit equal B1/B2 virtual spaces preserve the two Pi components.
my @virtual = split /,/, $ENV{'UKRMOL_VIRTUAL_ORBITALS'};
$model{'virtual_orbs'} = [@virtual,0,0,0,0];
$model{'nvirtual'} = 0;
$model{'nvirtual'} += $_ foreach @virtual;
$model{'reference_orbs'} = [0,0,0,0,0,0,0,0];
$model{'ncasscf_states'} = {'singlet' => [1,0,0,0,0,0,0,0]};
$model{'ntarget_states'} = {'singlet' => [1,0,0,0,0,0,0,0]};
$model{'ntarget_states_used'} = 1;
if ($model{'model'} eq 'CAS-A') {
    my $roots = $ENV{'UKRMOL_TARGET_ROOTS'};
    $model{'ntarget_states'} = {
        'singlet' => [($roots) x 4,0,0,0,0],
        'triplet' => [($roots) x 4,0,0,0,0],
    };
    $model{'ntarget_states_used'} = $ENV{'UKRMOL_TARGET_STATES_USED'};
}
# B1/B2 contain the two 2Pi components and higher odd angular projections.
$model{'scattering_states'} = {'doublet' => [0,1,1,0,0,0,0,0]};
$model{'use_GTO'} = 1;
$model{'radius_GTO'} = $ENV{'UKRMOL_RADIUS'};
$model{'rmatrix_radius'} = $ENV{'UKRMOL_RADIUS'};
$model{'maxl_GTO'} = $ENV{'UKRMOL_MAXL'};
$model{'use_BTO'} = 0;
$model{'delthres'} = [($ENV{'UKRMOL_DELETION_THRESHOLD'}) x 8];
$model{'e_unit'} = 2;  # UKRmol input/output electron energies in eV.
$model{'r_unit'} = 0;  # Geometry and propagation radii in bohr.
$model{'x_unit'} = 1;  # Symmetry-resolved cross sections in bohr^2.
$model{'nescat'} = $ENV{'UKRMOL_ENERGIES'};
$model{'einc'} = "$ENV{'UKRMOL_ENERGY_START'}, $ENV{'UKRMOL_ENERGY_STEP'}";
$model{'raf'} = $ENV{'UKRMOL_PROPAGATION_RADIUS'};
$model{'max_multipole'} = 2;

$run{'molpro'} = 0;
$run{'psi4'} = 1;
$run{'print_info'} = 'both';
$run{'buffer_size'} = 512;
$run{'mpi_integrals'} = "mpirun -np $ENV{'UKRMOL_RANKS'}";
$run{'mpi_scatci'} = "mpirun -np $ENV{'UKRMOL_RANKS'}";
$run{'mpi_rsolve'} = "mpirun -np $ENV{'UKRMOL_RANKS'}";
$run{'parallel_geom'} = 1;
$run{'parallel_symm'} = 1;
$run{'clean'} = 0;
$run{'remove_moints'} = 0;
$run{'save_channels'} = 1;
$run{'save_rmat_amp'} = 1;
$run{'save_Kmatrix'} = 1;
$run{'only'} = '';

my $R = $ENV{'UKRMOL_BOND_LENGTH'};
%geometry = (
    'suffix' => '',
    'geometry_labels' => 'R_bohr',
    'correct_cm' => 1,
    'length_unit' => 0,
    'start_at_geometry' => 1,
    'stop_at_geometry' => 1,
    'geometries' => [{
        'description' => "$R",
        'gnuplot_desc' => "R = $R bohr",
        'atoms' => [['C',0,0,0], ['O',0,0,$R]],
    }],
);
