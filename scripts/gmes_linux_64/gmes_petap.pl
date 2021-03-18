#!/usr/bin/env perl
# ==============================================================
# Alex Lomsadze
#
# GeneMark-ES Suite version 4.*
#
# Last modified: January, 2020
#
# Please report problems to:
#   Alex Lomsadze  alexl@gatech.edu
#   Mark Borodovsky  borodovsky@gatech.edu
#
# Copyright:
#   Georgia Institute of Technology, Atlanta, Georgia, USA
#
# Dr. Mark Borodovsky Bioinformatics Lab
#
# Affiliation:
#   Joint Georgia Tech and Emory Wallace H Coulter Department of Biomedical Engineering
#   Center for Bioinformatics and Computational Genomics at Georgia Tech
#   School of Computational Science and Engineering at Georgia Tech
#
# This eukaryotic gene prediction suite contains:
#   GeneMark.hmm,  GeneMark-ES,  GeneMark-ET,  GeneMark-EP  and  *.PLUS algorithms
#
# Name of this package:  gmes_petap.pl
#   GeneMark.hmm  -> gm
#   Eukaryotic    -> e
#   Self-training -> s
#   Plus          -> p
#   Evidence      -> e
#   Transcripts   -> t
#     and         -> a
#   Proteins      -> p
#
# ==============================================================
# Algorithms included into this package were described in the following publications:
#
# Gene prediction algorithm GeneMark-EP
#    Bruna T., Lomsadze A. Borodovsky M.
#    "GeneMark-EP and -EP+: automatic eukaryotic gene prediction
#     supported by spliced aligned proteins"
#    in preparation
#
# Gene prediction algorithm GeneMark-ET
#    Lomsadze A., Burns P. and Borodovsky M.
#    "Integration of RNA-Seq Data into Eukaryotic Gene Finding Algorithm
#     with Semi-Supervised Training."
#    Nucleic Acids Research, 2014, July 2
#
# Gene prediction algorithms GeneMark.hmm and GeneMark-ES with Branch Point model
#    Ter-Hovhannisyan V., Lomsadze A., Chernoff Y. and Borodovsky M.
#    "Gene prediction in novel fungal genomes using an ab initio
#     algorithm with unsupervised training."
#    Genome Research, 2008, Dec 18(12):1979-90
#
# Gene prediction algorithm GeneMark-ES
#    Lomsadze A., Ter-Hovhannisyan V., Chernoff Y. and Borodovsky M.
#    "Gene identification in novel eukaryotic genomes by
#     self-training algorithm."
#    Nucleic Acids Research, 2005, Vol. 33, No. 20, 6494-6506
# ==============================================================

use strict;
use warnings;

# modules from standard Perl distribution
use Getopt::Long qw( GetOptions );
use FindBin qw( $Bin );
use File::Spec;
use File::Path qw( make_path );
use Cwd qw( abs_path cwd );
use Data::Dumper;

# modules from CPAN
use YAML;
use Hash::Merge qw( merge );
use Parallel::ForkManager;
use MCE::Mutex;

# some modules from CPAN are pre-installed locally
#  ReadBackwards
#  LineFit
use lib './lib';

# ------------------------------------------------
my $v = 0;
my $debug = 0;
my $cfg = {};  # parameters are stored here; reference to hash of hashes
my $log;
my $logMutex;
my $key_bin = 0;
# ------------------------------------------------
my $bin = $Bin;  # code directory
my $work_dir = cwd;  # use current working directory to store temporary and output files
# ------------------------------------------------

ReadParameters();

CreateDirectories()   if $cfg->{'Run'}->{'set_dirs'};
CommitData()          if $cfg->{'Run'}->{'commit_input_data'};
DataReport()          if $cfg->{'Run'}->{'input_data_report'};
CommitTrainingData()  if $cfg->{'Run'}->{'commit_training_data'};
TrainingDataReport()  if $cfg->{'Run'}->{'training_data_report'};

if( $cfg->{'Parameters'}->{'predict_with'} )
{
	RunPredictionWithModel();
}
elsif( $cfg->{'Parameters'}->{'ES'} )
{
	RunES();
}
elsif( $cfg->{'Parameters'}->{'ET'} )
{
	RunET();
}
elsif( $cfg->{'Parameters'}->{'EP'} )
{
	if ( $cfg->{'Parameters'}->{'EP'} =~ /\S/ )
	{
		RunEP();
	}
	else
	{
		RunES();
		RenameOutput("es");
		RunProtHint( "genemark_es.gtf" );
		UpdateSettingAfterProtHint();
		CleanOutputFolder();
		$cfg->{'Parameters'}->{'ini_mod'} = '';
		RunEP();
	}
}
elsif( $cfg->{'Parameters'}->{'ETP'} )
{
	# under development

	RunET();
	RenameOutput("et");
	RunProtHint( "genemark_et.gtf" );
	UpdateSettingAfterProtHint();

	# create T&P hints for training
	# create T&P hints for evidence
	# create T hints for evidence
	# merge above hints for ETP run
	
	CleanOutputFolder();
	$cfg->{'Parameters'}->{'ini_mod'} = ResolvePath( "gmhmm_et.mod", $work_dir );
	RunETP();
}
else
	{ print  "error, expected command was not detected\n"; exit 1; }

PredictionReport() if $cfg->{'Run'}->{'prediction_report'};

close($log);

exit 0;

# ================= subs =========================
sub UpdateSettingAfterProtHint
{
	chdir $work_dir;

	my $ep_input    = "prothint/prothint.gff";
	my $ep_evidence = "prothint/evidence.gff";

	$ep_input = ResolvePath( $ep_input, $work_dir );
	$ep_evidence = ResolvePath( $ep_evidence, $work_dir );

	$cfg->{'Parameters'}->{'EP'} = $ep_input;
	$cfg->{'Parameters'}->{'evidence'} = $ep_evidence;

	CommitEvidence();	
	CommitDataEP();

	CommitTrainingEvidence();
	CommitTrainingEP();
}
# ------------------------------------------------
sub CleanOutputFolder
{
	chdir $work_dir;

	if ( ! -d "output" )
		{ print "error, folder not found: output\n"; exit 1; }

	chdir "output";

	my $dir = "gmhmm";

	opendir( DIR, $dir ) or die "error on open directory $dir: $!\n";
	foreach my $file ( grep{ /dna.fa_\d+$/ } readdir(DIR) )
	{
		unlink "$dir/$file" or die "error on delete file $file: $!\n";
	}
	closedir DIR;

	$dir = "data";

	opendir( DIR, $dir ) or die "error on open directory $dir: $!\n";
	foreach my $file ( grep{ /dna.fa_\d+.out$/ } readdir(DIR) )
	{
		unlink "$dir/$file" or die "error on delete $file: $!\n";
	}
	closedir DIR;

	if ( -e "gmhmm.mod" )
	{
		unlink "gmhmm.mod" or die "error on delete gmhmm.mod: $!\n";
	}

	chdir $work_dir;
}
# ------------------------------------------------
sub RunProtHint
{
	my $gtf = shift;

	$gtf = ResolvePath( $gtf, $work_dir );

	chdir $work_dir;

	mkdir "prothint";
	if ( ! -d "prothint" )
		{ print "error, folder not found: prothint\n"; exit 1; }

	chdir "prothint";

	my $com = "$Bin/ProtHint/bin/prothint.py  --geneMarkGtf $gtf  $cfg->{'Parameters'}->{'sequence'}  $cfg->{'Parameters'}->{'dbep'} ";

	if ( $cfg->{'Parameters'}->{'pbs'} )
	{
		$com .= " --pbs ";
	}
	else
	{
		$com .= " --threads $cfg->{'Parameters'}->{'cores'} ";
	}

	RunCom( $com );

	chdir $work_dir;
}
# ------------------------------------------------
sub RenameOutput
{
	my $label = shift;

	chdir $work_dir;

	if ( -e "genemark.gtf" )
	{
		RunCom( "mv genemark.gtf genemark_". $label .".gtf" );
	}
	else
		{ print "error, file not found: genemark.gtf\n"; exit 1; }

	if ( -e "gmhmm.mod" )
	{
		RunCom( "mv gmhmm.mod gmhmm_". $label .".mod" );
	}
	else
		{ print "error, file not found: gmhmm.mod\n"; exit 1; }
}
# ------------------------------------------------
sub RunPredictionWithModel
{
	PredictGenes( $cfg->{'Parameters'}->{'predict_with'}, $cfg->{'Parameters'}->{'min_gene_in_predict'} );
}
# ------------------------------------------------
sub RunES
{
	PrepareInitialModel( \&BuildInitialModelES ) if $cfg->{'Run'}->{'prepare_ini_mod'};

	if (  $cfg->{'Run'}->{'run_training'} )
	{
		RunIterations('ES_A', \&Training_ES_A ) if $cfg->{'ES_A'}->{'iterations'};
		RunIterations('ES_B', \&Training_ES_B ) if $cfg->{'ES_B'}->{'iterations'};
		RunIterations('ES_C', \&Training_ES_C ) if $cfg->{'ES_C'}->{'iterations'};

		my $final_mod = ResolvePath( $cfg->{'ES_C'}->{ 'out_mod'}, "run" );
		RunCom( "cp $final_mod $work_dir/gmhmm.mod" );
	}

	TrainingReportES() if $cfg->{'Run'}->{'training_report'};
	PredictGenes( "$work_dir/gmhmm.mod", $cfg->{'Parameters'}->{'min_gene_in_predict'} ) if $cfg->{'Run'}->{'run_prediction'};
}
# ------------------------------------------------
sub RunET
{
	PrepareInitialModel( \&BuildInitialModelET ) if $cfg->{'Run'}->{'prepare_ini_mod'};

	if (  $cfg->{'Run'}->{'run_training'} )
	{
		RunIterations('ET_A', \&Training_ET_A ) if $cfg->{'ET_A'}->{'iterations'};
		RunIterations('ET_B', \&Training_ET_B ) if $cfg->{'ET_B'}->{'iterations'};
		RunIterations('ET_C', \&Training_ET_C ) if $cfg->{'ET_C'}->{'iterations'};

		my $final_mod = ResolvePath( $cfg->{'ET_C'}->{ 'out_mod'}, "run" );
		RunCom( "cp $final_mod $work_dir/gmhmm.mod" );
	}

	TrainingReportET() if $cfg->{'Run'}->{'training_report'};
	PredictGenes( "$work_dir/gmhmm.mod", $cfg->{'Parameters'}->{'min_gene_in_predict'} ) if $cfg->{'Run'}->{'run_prediction'};
}
# ------------------------------------------------
sub RunEP
{
	PrepareInitialModel( \&BuildInitialModelEP ) if $cfg->{'Run'}->{'prepare_ini_mod'};

	if (  $cfg->{'Run'}->{'run_training'} )
	{
		RunIterations('EP_A', \&Training_EP_A ) if $cfg->{'EP_A'}->{'iterations'};
		RunIterations('EP_B', \&Training_EP_B ) if $cfg->{'EP_B'}->{'iterations'};
		RunIterations('EP_C', \&Training_EP_C ) if $cfg->{'EP_C'}->{'iterations'};

		my $final_mod = ResolvePath( $cfg->{'EP_C'}->{ 'out_mod'}, "run" );
		RunCom( "cp $final_mod $work_dir/gmhmm.mod" );
	}

	TrainingReportEP() if $cfg->{'Run'}->{'training_report'};
	PredictGenes( "$work_dir/gmhmm.mod", $cfg->{'Parameters'}->{'min_gene_in_predict'} ) if $cfg->{'Run'}->{'run_prediction'};
}
# ------------------------------------------------
sub RunETP
{
	PrepareInitialModel( \&BuildInitialModelETP ) if $cfg->{'Run'}->{'prepare_ini_mod'};

	if (  $cfg->{'Run'}->{'run_training'} )
	{
		RunIterations('ETP_C', \&Training_ETP_C ) if $cfg->{'ETP_C'}->{'iterations'};

		my $final_mod = ResolvePath( $cfg->{'ETP_C'}->{ 'out_mod'}, "run" );
		RunCom( "cp $final_mod $work_dir/gmhmm.mod" );
	}

	TrainingReportETP() if $cfg->{'Run'}->{'training_report'};
	PredictGenes( "$work_dir/gmhmm.mod", $cfg->{'Parameters'}->{'min_gene_in_predict'} ) if $cfg->{'Run'}->{'run_prediction'};
}
# ------------------------------------------------
sub BuildInitialModelEP
{
        print "build initial EP model\n" if $v;

        my $mod = shift;
        chdir $work_dir;

        my $dir_for_build = "run/EP_ini";
        $dir_for_build = SetDir( $dir_for_build );

        chdir "data";
        RunCom("$bin/parse_by_introns.pl  --section EP_ini  --cfg  $cfg->{'Config'}->{'run_cfg'}  --parse_dir $dir_for_build" );

        chdir $dir_for_build;
        RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section donor_GT    --format DONOR" );
        RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section acceptor_AG --format ACCEPTOR" );
        RunCom( "$bin/histogram.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section intron_DUR" );

        if( $cfg->{'Parameters'}->{'fungus'} )
        {
                RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section acceptor_short_AG --format ACC_BP " );

                # run Gibss3 on subset of introns
                RunCom( "$bin/bp_seq_select.pl --seq_in $cfg->{'EP_ini'}->{'bp_region'} --seq_out $cfg->{'EP_ini'}->{'gibbs_seq'}  --max_seq $cfg->{'EP_ini'}->{'gibbs_seq_max'}  --bp_region_length  $cfg->{'EP_ini'}->{'bp_region_length'}  --min_bp_region_length $cfg->{'EP_ini'}->{'min_bp_region_length'} " );

                # 9     motif length
                # -n    Use nucleic acid alphabet
                # -r    turn off reverse complements with DNA
                # -nopt Don't print Near Optimal output
                # -m    Do not maximize after near optimal sampling
                # -w    pseduocount weight
                # -Z    Don't write progress info
                # -s    random number generator seed

                RunCom( "$bin/Gibbs3  gibbs.seq  9 -n -r -o gibbs.out -nopt -m -w 0.001 -Z  -s 1 -P $bin/prior.bp -F" );

                RunCom( "$bin/parse_gibbs.pl --seq gibbs.seq  --gibbs gibbs.out --motif_seq $cfg->{'branch_point'}->{'infile'}  --spacer_len $cfg->{'spacer_DUR'}->{'in'}  --spacer_seq spacer.seq  --tr $cfg->{'EP_ini'}->{'tr_bp'}  "  );
                RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section branch_point --format BRANCH " );
                RunCom( "$bin/histogram.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section spacer_DUR" );

                RunCom( "$bin/scan_for_bp.pl --seq_in $cfg->{'EP_ini'}->{'bp_region'}  --gibbs_in gibbs.out  --pos_out $cfg->{'prespacer_DUR'}->{'in'} ");
                RunCom( "$bin/histogram.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section prespacer_DUR" );

                my $str = " --MKCHAIN_L_MARGING 0  --MKCHAIN_R_MARGING 0  --MKCHAIN_PSEUDOCOUNTS 1  --MKCHAIN_PRECISION 6 --ORDM 1 ";
                RunCom( "$bin/probuild --non spacer.seq --mkmod_non spacer.mkch  $str" );
        }

        RunCom( "$bin/build_mod.pl --cfg $cfg->{'Config'}->{'run_cfg'}  --section EP_ini --def $mod ");

        $mod = ResolvePath( $cfg->{'EP_ini'}->{'mod'} );

        chdir $work_dir;

        return $mod;
}
# ------------------------------------------------
sub BuildInitialModelET
{
	print "build initial ET model\n" if $v;
	
	my $mod = shift;
	chdir $work_dir;

	my $dir_for_build = "run/ET_ini";
	$dir_for_build = SetDir( $dir_for_build );
	
	chdir "data";
	RunCom("$bin/parse_by_introns.pl  --section ET_ini  --cfg  $cfg->{'Config'}->{'run_cfg'}  --parse_dir $dir_for_build" );
	
	chdir $dir_for_build;
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section donor_GT    --format DONOR " ); 
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section acceptor_AG --format ACCEPTOR " );
	RunCom( "$bin/histogram.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section intron_DUR" );

	if( $cfg->{'Parameters'}->{'fungus'} )
	{
		RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section acceptor_short_AG --format ACC_BP " );

		# run Gibss3 on subset of introns
		RunCom( "$bin/bp_seq_select.pl --seq_in $cfg->{'ET_ini'}->{'bp_region'} --seq_out $cfg->{'ET_ini'}->{'gibbs_seq'}  --max_seq $cfg->{'ET_ini'}->{'gibbs_seq_max'}  --bp_region_length  $cfg->{'ET_ini'}->{'bp_region_length'}  --min_bp_region_length $cfg->{'ET_ini'}->{'min_bp_region_length'} " );
		
		# 9     motif length
		# -n    Use nucleic acid alphabet
		# -r    turn off reverse complements with DNA
		# -nopt Don't print Near Optimal output
		# -m    Do not maximize after near optimal sampling
		# -w    pseduocount weight
		# -Z    Don't write progress info
		# -s    random number generator seed
		
		RunCom( "$bin/Gibbs3  gibbs.seq  9 -n -r -o gibbs.out -nopt -m -w 0.001 -Z  -s 1 -P $bin/prior.bp -F" );

		RunCom( "$bin/parse_gibbs.pl --seq gibbs.seq  --gibbs gibbs.out --motif_seq $cfg->{'branch_point'}->{'infile'}  --spacer_len $cfg->{'spacer_DUR'}->{'in'}  --spacer_seq spacer.seq  --tr $cfg->{'ET_ini'}->{'tr_bp'}  "  );
		RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section branch_point --format BRANCH " );
		RunCom( "$bin/histogram.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section spacer_DUR" );
		
		RunCom( "$bin/scan_for_bp.pl --seq_in $cfg->{'ET_ini'}->{'bp_region'}  --gibbs_in gibbs.out  --pos_out $cfg->{'prespacer_DUR'}->{'in'} ");
		RunCom( "$bin/histogram.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section prespacer_DUR" );

		my $str = " --MKCHAIN_L_MARGING 0  --MKCHAIN_R_MARGING 0  --MKCHAIN_PSEUDOCOUNTS 1  --MKCHAIN_PRECISION 6 --ORDM 1 ";
		RunCom( "$bin/probuild --non spacer.seq --mkmod_non spacer.mkch  $str" );
	}

	RunCom( "$bin/build_mod.pl --cfg $cfg->{'Config'}->{'run_cfg'}  --section ET_ini --def $mod ");
	
	$mod = ResolvePath( $cfg->{'ET_ini'}->{'mod'} );
	
	chdir $work_dir;
	
	return $mod;
}
# ------------------------------------------------
sub BuildInitialModelES
{
	print "build initial ES model\n" if $v;
	
	my $mod = shift;
	chdir $work_dir;
	
	my $dir_for_build = "run/ES_ini";
	$dir_for_build = SetDir( $dir_for_build );
	
	chdir $dir_for_build;
	RunCom( "$bin/build_mod.pl --cfg $cfg->{'Config'}->{'run_cfg'}  --section ES_ini --def $mod ");
	
	$mod = ResolvePath( $cfg->{'ES_ini'}->{'mod'} );
	
	chdir $work_dir;
	
	return $mod;
}
# ------------------------------------------------
sub Training_EP_C
{
	my ( $path, $name ) = @_;

	Training_E_anchored_C($path, $name, "EP_C");
}
# ------------------------------------------------
sub Training_ET_C
{
	my ( $path, $name ) = @_;

	Training_E_anchored_C($path, $name, "ET_C");
}
# ------------------------------------------------
sub Training_ETP_C
{
	my ( $path, $name ) = @_;

	Training_E_anchored_C($path, $name, "ETP_C");
}
# ------------------------------------------------
sub Training_E_anchored_C
{
	my ( $path, $name, $section ) = @_;
	print "training level $section: $path\n" if $v;
	
	chdir $path; 
	RunCom("$bin/parse_ET.pl --section $section --cfg  $cfg->{'Config'}->{'run_cfg'}  --v" );
	
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section start_ATG  --format INI" ); 
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section stop_TAA   --format TERM_TAA" );
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section stop_TAG   --format TERM_TAG" );
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section stop_TGA   --format TERM_TGA" );

	if ( $cfg->{'Parameters'}->{'gc3'} < 100 )
	{
		RunCom( "mv cod.seq cod.seq_all"  );
		RunCom( "$bin/get_below_gc.pl  < cod.seq_all > cod.seq" );
	}

	my $str = " --MKCHAIN_L_MARGING 0  --MKCHAIN_R_MARGING 0  --MKCHAIN_PSEUDOCOUNTS 1  --MKCHAIN_PRECISION 8";
	$str .= " --revcomp_non  --ORDM 5 ";
	RunCom( "$bin/probuild --cod cod.seq --non non.seq --mkmod_euk mkch   $str" );
	
	RunCom( "$bin/histogram.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section initial_DUR" );
	RunCom( "$bin/histogram.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section internal_DUR" );
	RunCom( "$bin/histogram.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section terminal_DUR" );
	RunCom( "$bin/histogram.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section single_DUR" );

	if ( $cfg->{'Parameters'}->{'fungus'} and ($cfg->{'Fungi'}->{'intergenic_type'} eq "uniform") )
	{
		;
	}
	else
	{
		RunCom( "$bin/histogram.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section intergenic_DUR" );
	}

	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section donor_GT    --format DONOR_0    --phase 0  --quiet" ); 
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section acceptor_AG --format ACCEPTOR_0 --phase 0 " );
	RunCom( " cat  GT.mat > donor.mat " );
	RunCom( " cat  AG.mat > acceptor.mat " );
	
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section donor_GT    --format DONOR_1    --phase 1  --quiet" ); 
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section acceptor_AG --format ACCEPTOR_1 --phase 1 " );
	RunCom( " cat  GT.mat >> donor.mat " );
	RunCom( " cat  AG.mat >> acceptor.mat " );
	
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section donor_GT    --format DONOR_2    --phase 2  --quiet" ); 
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section acceptor_AG --format ACCEPTOR_2 --phase 2 " );
	RunCom( " cat  GT.mat >> donor.mat " );
	RunCom( " cat  AG.mat >> acceptor.mat " );
	
	RunCom( " mv donor.mat     GT.mat " );
	RunCom( " mv acceptor.mat  AG.mat " );

	if( $cfg->{'Parameters'}->{'gc_donor'} and ($cfg->{'Parameters'}->{'gc_donor'} ne "auto") )
	{
		RunCom( "$bin/make_nt_freq_mat.pl  --cfg $cfg->{'Config'}->{'run_cfg'}  --section donor_GC  --format DONOR_GC_0  --phase 0  --force" );
		RunCom( " cat  GC.mat > donor_GC.mat " );
		RunCom( "$bin/make_nt_freq_mat.pl  --cfg $cfg->{'Config'}->{'run_cfg'}  --section donor_GC  --format DONOR_GC_1  --phase 1  --force" );
		RunCom( " cat  GC.mat >> donor_GC.mat " );
		RunCom( "$bin/make_nt_freq_mat.pl  --cfg $cfg->{'Config'}->{'run_cfg'}  --section donor_GC  --format DONOR_GC_2  --phase 2  --force" );
		RunCom( " cat  GC.mat >> donor_GC.mat " );
		RunCom( " mv donor_GC.mat GC.mat " );
	}
	
	if( $cfg->{'Parameters'}->{'fungus'} )
	{
		RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section acceptor_short_AG --format ACC_BP_0 --phase 0 " );
		RunCom( " cat AG_SHORT.mat >  acceptor_short.mat " );
		RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section acceptor_short_AG --format ACC_BP_1 --phase 1 " );
		RunCom( " cat AG_SHORT.mat >> acceptor_short.mat " );
		RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section acceptor_short_AG --format ACC_BP_2 --phase 2 " );
		RunCom( " cat AG_SHORT.mat >> acceptor_short.mat " );
		RunCom( " mv acceptor_short.mat  AG_SHORT.mat " );
	}
	
	RunCom( "$bin/build_mod.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section $section --def prev.mod  --out $name ");
	
	chdir $work_dir;
}
# ------------------------------------------------
sub Training_EP_B
{
	my ( $path, $name ) = @_;

	Training_E_anchored_B($path, $name, "EP_B");
}
# ------------------------------------------------
sub Training_ET_B
{
	my ( $path, $name ) = @_;

	Training_E_anchored_B($path, $name, "ET_B");
}
# ------------------------------------------------
sub Training_E_anchored_B
{
	my ( $path, $name, $section ) = @_;
	print "training level $section: $path\n" if $v;

	chdir $path; 
	RunCom("$bin/parse_ET.pl --section $section --cfg  $cfg->{'Config'}->{'run_cfg'}  --v" );
	
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section start_ATG  --format INI" ); 
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section stop_TAA   --format TERM_TAA" );
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section stop_TAG   --format TERM_TAG" );
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section stop_TGA   --format TERM_TGA" );	

	if ( $cfg->{'Parameters'}->{'gc3'} < 100 )
	{
		RunCom( "mv cod.seq cod.seq_all"  );
		RunCom( "$bin/get_below_gc.pl  < cod.seq_all > cod.seq" );
	}

	my $str = " --MKCHAIN_L_MARGING 0  --MKCHAIN_R_MARGING 0  --MKCHAIN_PSEUDOCOUNTS 1  --MKCHAIN_PRECISION 8";
	$str .= " --revcomp_non  --ORDM 5 ";
	RunCom( "$bin/probuild --cod cod.seq --non non.seq --mkmod_euk mkch   $str" );
	
	RunCom( "$bin/histogram.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section initial_DUR" );
	RunCom( "$bin/histogram.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section internal_DUR" );
	RunCom( "$bin/histogram.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section terminal_DUR" );
	RunCom( "$bin/histogram.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section single_DUR" );
	
	RunCom( "$bin/build_mod.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section $section --def prev.mod  --out $name ");
	
	chdir $work_dir;
}
# ------------------------------------------------
sub Training_EP_A
{
	my ( $path, $name ) = @_;

	Training_E_anchored_A($path, $name, "EP_A");
}
# ------------------------------------------------
sub Training_ET_A
{
	my ( $path, $name ) = @_;

	Training_E_anchored_A($path, $name, "ET_A");
}
# ------------------------------------------------
sub Training_E_anchored_A
{
	my ( $path, $name, $section ) = @_;
	print "training level $section: $path\n" if $v;
	
	chdir $path; 
	RunCom("$bin/parse_ET.pl --section $section --cfg  $cfg->{'Config'}->{'run_cfg'}  --v" );

	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section start_ATG  --format INI" ); 
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section stop_TAA   --format TERM_TAA" );
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section stop_TAG   --format TERM_TAG" );
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section stop_TGA   --format TERM_TGA" );

	if ( $cfg->{'Parameters'}->{'gc3'} < 100 )
	{
		RunCom( "mv cod.seq cod.seq_all"  );
		RunCom( "$bin/get_below_gc.pl  < cod.seq_all > cod.seq" );
	}

	my $str = " --MKCHAIN_L_MARGING 0  --MKCHAIN_R_MARGING 0  --MKCHAIN_PSEUDOCOUNTS 1  --MKCHAIN_PRECISION 8";
	$str .= " --revcomp_non  --ORDM 5 ";
	RunCom( "$bin/probuild --cod cod.seq --non non.seq --mkmod_euk mkch   $str" );
	
	RunCom( "$bin/build_mod.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section $section --def prev.mod  --out $name ");
	
	chdir $work_dir;
}
# ------------------------------------------------
sub Training_ES_C
{
	my ( $path, $name ) = @_;
	print "training level ES_C: $path\n" if $v;
	chdir $work_dir;

	chdir $path;
	RunCom("$bin/parse_set.pl --section ES_C --cfg  $cfg->{'Config'}->{'run_cfg'}  --v" );

	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section start_ATG   --format INI" ); 
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section stop_TAA    --format TERM_TAA" );
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section stop_TAG    --format TERM_TAG" );
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section stop_TGA    --format TERM_TGA" );

	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section donor_GT    --format DONOR_0    --phase 0  --quiet" ); 
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section acceptor_AG --format ACCEPTOR_0 --phase 0 " );
	RunCom( " cat  GT.mat > donor.mat " );
	RunCom( " cat  AG.mat > acceptor.mat " );
	
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section donor_GT    --format DONOR_1    --phase 1  --quiet" ); 
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section acceptor_AG --format ACCEPTOR_1 --phase 1 " );
	RunCom( " cat  GT.mat >> donor.mat " );
	RunCom( " cat  AG.mat >> acceptor.mat " );
	
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section donor_GT    --format DONOR_2    --phase 2  --quiet" ); 
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section acceptor_AG --format ACCEPTOR_2 --phase 2 " );
	RunCom( " cat  GT.mat >> donor.mat " );
	RunCom( " cat  AG.mat >> acceptor.mat " );
	
	RunCom( " mv donor.mat     GT.mat " );
	RunCom( " mv acceptor.mat  AG.mat " );

	if( $cfg->{'Parameters'}->{'gc_donor'} and ($cfg->{'Parameters'}->{'gc_donor'} ne "auto") )
	{
		RunCom( "$bin/make_nt_freq_mat.pl  --cfg $cfg->{'Config'}->{'run_cfg'}  --section donor_GC  --format DONOR_GC_0  --phase 0  --force" );
		RunCom( " cat  GC.mat > donor_GC.mat " );
		RunCom( "$bin/make_nt_freq_mat.pl  --cfg $cfg->{'Config'}->{'run_cfg'}  --section donor_GC  --format DONOR_GC_1  --phase 1  --force" );
		RunCom( " cat  GC.mat >> donor_GC.mat " );
		RunCom( "$bin/make_nt_freq_mat.pl  --cfg $cfg->{'Config'}->{'run_cfg'}  --section donor_GC  --format DONOR_GC_2  --phase 2  --force" );
		RunCom( " cat  GC.mat >> donor_GC.mat " );
		RunCom( " mv donor_GC.mat GC.mat " );
	}

	if ( $cfg->{'Parameters'}->{'gc3'} < 100 )
	{
		RunCom( "mv cod.seq cod.seq_all"  );
		RunCom( "$bin/get_below_gc.pl  < cod.seq_all > cod.seq" );
	}

	my $str = " --MKCHAIN_L_MARGING 0  --MKCHAIN_R_MARGING 0  --MKCHAIN_PSEUDOCOUNTS 1  --MKCHAIN_PRECISION 8";
	$str .= " --revcomp_non  --ORDM 5 ";
	RunCom( "$bin/probuild --cod cod.seq --non non.seq --mkmod_euk mkch   $str" );
	
	RunCom( "$bin/histogram.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section intron_DUR" );
	RunCom( "$bin/histogram.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section initial_DUR" );
	RunCom( "$bin/histogram.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section internal_DUR" );
	RunCom( "$bin/histogram.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section terminal_DUR" );
	RunCom( "$bin/histogram.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section single_DUR" );

	if ( $cfg->{'Parameters'}->{'fungus'} and ($cfg->{'Fungi'}->{'intergenic_type'} eq "uniform") ) 
	{
		;
	}
	else
	{
		RunCom( "$bin/histogram.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section intergenic_DUR" );
	}

	if( $cfg->{'Parameters'}->{'fungus'} )
	{
		RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section acceptor_short_AG --format ACC_BP_0 --phase 0 " );
		RunCom( " cat AG_SHORT.mat >  acceptor_short.mat " );
		RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section acceptor_short_AG --format ACC_BP_1 --phase 1 " );
		RunCom( " cat AG_SHORT.mat >> acceptor_short.mat " );
		RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section acceptor_short_AG --format ACC_BP_2 --phase 2 " );
		RunCom( " cat AG_SHORT.mat >> acceptor_short.mat " );
		RunCom( " mv acceptor_short.mat  AG_SHORT.mat " );
		
		# run Gibss3 on subset of introns
		RunCom( "$bin/bp_seq_select.pl --seq_in $cfg->{'ES_C'}->{'bp_region'} --seq_out $cfg->{'ES_C'}->{'gibbs_seq'}  --max_seq $cfg->{'ES_C'}->{'gibbs_seq_max'}  --bp_region_length  $cfg->{'ES_C'}->{'bp_region_length'}  --min_bp_region_length $cfg->{'ET_ini'}->{'min_bp_region_length'} ");
		
		# 9     motif length
		# -n    Use nucleic acid alphabet
		# -r    turn off reverse complements with DNA
		# -nopt Don't print Near Optimal output
		# -m    Do not maximize after near optimal sampling
		# -w    pseduocount weight
		# -Z    Don't write progress info
		# -s    random number generator seed
		
		RunCom( "$bin/Gibbs3  gibbs.seq  9 -n -r -o gibbs.out -nopt -m -w 0.001 -Z  -s 1 -P $bin/prior.bp -F" );

		RunCom( "$bin/parse_gibbs.pl --seq gibbs.seq  --gibbs gibbs.out --motif_seq $cfg->{'branch_point'}->{'infile'}  --spacer_len $cfg->{'spacer_DUR'}->{'in'}  --spacer_seq spacer.seq  --tr $cfg->{'ES_C'}->{'tr_bp'}  "  );
		RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section branch_point --format BRANCH " );
		RunCom( "$bin/histogram.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section spacer_DUR" );
		
		RunCom( "$bin/scan_for_bp.pl --seq_in $cfg->{'ES_C'}->{'bp_region'}  --gibbs_in gibbs.out  --pos_out $cfg->{'prespacer_DUR'}->{'in'} ");
		RunCom( "$bin/histogram.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section prespacer_DUR" );

		my $str = " --MKCHAIN_L_MARGING 0  --MKCHAIN_R_MARGING 0  --MKCHAIN_PSEUDOCOUNTS 1  --MKCHAIN_PRECISION 6 --ORDM 1 ";
		RunCom( "$bin/probuild --non spacer.seq --mkmod_non spacer.mkch  $str" );
	}
	
	RunCom( "$bin/build_mod.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section ES_C --def prev.mod  --out $name ");
	
	chdir $work_dir;
}
# ------------------------------------------------
sub Training_ES_B
{
	my ( $path, $name ) = @_;
	print "training level ES_B: $path\n" if $v;
	chdir $work_dir;
	
	chdir $path; 
	RunCom("$bin/parse_set.pl --section ES_B --cfg  $cfg->{'Config'}->{'run_cfg'}  --v" );

	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section start_ATG   --format INI" ); 
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section stop_TAA    --format TERM_TAA" );
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section stop_TAG    --format TERM_TAG" );
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section stop_TGA    --format TERM_TGA" );

	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section donor_GT    --format DONOR_0    --phase 0 " ); 
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section acceptor_AG --format ACCEPTOR_0 --phase 0 " );
	RunCom( " cat  GT.mat > donor.mat " );
	RunCom( " cat  AG.mat > acceptor.mat " );
	
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section donor_GT    --format DONOR_1    --phase 1 " ); 
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section acceptor_AG --format ACCEPTOR_1 --phase 1 " );
	RunCom( " cat  GT.mat >> donor.mat " );
	RunCom( " cat  AG.mat >> acceptor.mat " );
	
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section donor_GT    --format DONOR_2    --phase 2 " ); 
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section acceptor_AG --format ACCEPTOR_2 --phase 2 " );
	RunCom( " cat  GT.mat >> donor.mat " );
	RunCom( " cat  AG.mat >> acceptor.mat " );
	
	RunCom( " mv donor.mat     GT.mat " );
	RunCom( " mv acceptor.mat  AG.mat " );

	if ( $cfg->{'Parameters'}->{'gc3'} < 100 )
	{
		RunCom( "mv cod.seq cod.seq_all"  );
		RunCom( "$bin/get_below_gc.pl  < cod.seq_all > cod.seq" );
	}

	my $str = " --MKCHAIN_L_MARGING 0  --MKCHAIN_R_MARGING 0  --MKCHAIN_PSEUDOCOUNTS 1  --MKCHAIN_PRECISION 8";
	$str .= " --revcomp_non  --ORDM 5 ";
	RunCom( "$bin/probuild --cod cod.seq --non non.seq --mkmod_euk mkch   $str" );
	
	RunCom( "$bin/histogram.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section intron_DUR" );
	RunCom( "$bin/histogram.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section initial_DUR" );
	RunCom( "$bin/histogram.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section internal_DUR" );
	RunCom( "$bin/histogram.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section terminal_DUR" );
	RunCom( "$bin/histogram.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section single_DUR" );
	
	RunCom( "$bin/build_mod.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section ES_B --def prev.mod  --out $name ");
	
	chdir $work_dir;
}
# ------------------------------------------------
sub Training_ES_A
{
	my ( $path, $name ) = @_;
	print "training level ES_A: $path\n" if $v;
	chdir $work_dir;
	
	chdir $path; 
	RunCom("$bin/parse_set.pl --section ES_A --cfg  $cfg->{'Config'}->{'run_cfg'}  --v " );

	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section start_ATG   --format INI" ); 
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section stop_TAA    --format TERM_TAA" );
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section stop_TAG    --format TERM_TAG" );
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section stop_TGA    --format TERM_TGA" );

	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section donor_GT    --format DONOR " ); 
	RunCom( "$bin/make_nt_freq_mat.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section acceptor_AG --format ACCEPTOR " );

	if ( $cfg->{'Parameters'}->{'gc3'} < 100 )
	{
		RunCom( "mv cod.seq cod.seq_all"  );
		RunCom( "$bin/get_below_gc.pl  < cod.seq_all > cod.seq" );
	}

	my $str = " --MKCHAIN_L_MARGING 0  --MKCHAIN_R_MARGING 0  --MKCHAIN_PSEUDOCOUNTS 1  --MKCHAIN_PRECISION 8";
	$str .= " --revcomp_non  --ORDM 5 ";
	RunCom( "$bin/probuild --cod cod.seq --non non.seq --mkmod_euk mkch   $str" );
	
	RunCom( "$bin/build_mod.pl --cfg $cfg->{'Config'}->{'run_cfg'} --section ES_A --def prev.mod  --out $name ");
	
	chdir $work_dir;
}
# ------------------------------------------------
sub ConcatenatePredictions
{
	my ( $path ) = @_;
	print "concatenate predictions: $path\n" if $v;
	
	my $dir = ResolvePath( "hmmout", $path );
	
	unlink "$path/set.out" if ( -e "$path/set.out" );
	
	opendir( DIR, $dir );
	foreach  my $file ( grep{ /dna\.fa_\d+\.out$/ } readdir(DIR) )
	{
		$file = "$dir/$file";
		my $com = "cat $file >> $path/set.out";
		system( $com ) && die"$!\n";
		
		unlink $file;
	}
	closedir DIR;
}
# ------------------------------------------------
sub SetDir
{
	my $dir = shift;
	mkdir $dir;
	if ( ! -d $dir) { print "error, directory not found: $dir\n"; exit 1; }
	return abs_path( $dir );
}
# ------------------------------------------------
sub RunIterations
{
	my ( $label , $function ) = @_;
	print "running step $label\n" if $v;

	my $current_mod = ResolvePath( $cfg->{$label}->{'in_mod'}, 'run' );
	my $new_mod;

	my $end = $cfg->{$label}->{'iterations'};
	foreach my $i (1..$end)
	{
		my $training_dir = SetDir( "run/$label" ."_". $i );
		my $hmmout_dir   = SetDir( "$training_dir/hmmout" );
		
		if ( $cfg->{$label}->{'run_prediction'} )
		{
			my $plus = '';
			if ( $cfg->{'Parameters'}->{'evidence'} )
			{
				if ( -e "$work_dir/data/evidence_training.gff" and -s "$work_dir/data/evidence_training.gff" )
				{
					# Verify evidence for the current model
					RunCom( "$bin/verify_evidence_gmhmm.pl --in $work_dir/data/evidence_training.gff  --out $training_dir/evidence_hmm.gff  --mod $current_mod " );
					$plus = "$training_dir/evidence_hmm.gff";
				}
			}

			RunHmm( $current_mod, $hmmout_dir, '', '', $plus );
			ConcatenatePredictions( $training_dir );
		}
		
		$new_mod = "$label\_$i.mod";
		if (  $cfg->{$label}->{'run_training'} )
		{
			RunCom( "ln -sf $current_mod  $training_dir/prev.mod" );
			$function->( $training_dir, $new_mod );
		}
		$current_mod = ResolvePath( $new_mod, $training_dir );
	}
	
	my $out_mod = $cfg->{$label}->{'out_mod'};
	$out_mod = File::Spec->catfile( 'run', $out_mod );

	unlink $out_mod  if ( -e $out_mod );
	symlink $current_mod, $out_mod;
}
# ------------------------------------------------
sub PredictGenes
{
	# parameter-1: predict genes using this file with parameters
	# parameter-2: filter out genes shorter than minimum
	my( $name, $min ) = @_;

	print "predict final gene set\n" if $v;
	chdir $work_dir;
	
	# copy parameter file to output folder
	if( !$name )     { print "error, hmm parameter file not specified $0\n"; exit 1; }
	if( ! -e $name ) { print "error, hmm parameter file not found $0: $name\n"; exit 1; }
	RunCom( "cp $name output/gmhmm.mod" );
	
	# split sequence for prediction step
	my $dir = "output/data";
	mkdir $dir;
	SplitFasta( abs_path("data/dna.fna"), $dir, $cfg->{'Parameters'}->{'min_contig_in_predict'}, "prediction.trace" );
	chdir $work_dir;

	my $evi = '';

	# evidence
	if( $cfg->{'Parameters'}->{'evidence'} )
	{
		RunCom( "$bin/verify_evidence_gmhmm.pl --in data/evidence.gff  --out data/evidence_hmm.gff  --mod  output/gmhmm.mod " );
		RunCom( "$bin/rescale_gff.pl  --in data/evidence_hmm.gff  --trace info/prediction.trace  --out $dir/evidence_prediction.gff" );
		$evi = abs_path( "$dir/evidence_prediction.gff" );
	}
	
	# predict genes
	$dir = "output/gmhmm";
	mkdir $dir;
	
	RunHmm( abs_path("output/gmhmm.mod"), "output/gmhmm" , "output/data", ' -f lst -n ', $evi );

	# reformat from LST to GFF
	my @files = ReadContigNames( "output/data" );
	my $final =  "genemark.gtf";
	unlink $final;
	foreach my $f (@files)
	{
		$f = File::Spec->splitpath( $f ) .".out";
		$f = ResolvePath( $f, "output/gmhmm" );
		
		RunCom( "$bin/hmm_to_gtf.pl  --in $f  --app  --out $final  --min $min " );
	}

	RunCom( "$bin/reformat_gff.pl --out $final.tmp --trace info/dna.trace --in $final  --back" );
	RunCom( "mv $final.tmp $final" );

	if ( $cfg->{'Parameters'}->{'EP'} && $cfg->{'Parameters'}->{'EP'} =~ /\S/ )
	{
		RunCom( "$bin/flag_anchored_elements.py $final $cfg->{'Parameters'}->{'EP'} > $final.tmp" );
		RunCom( "mv $final.tmp $final" );
	}
	elsif ( $cfg->{'Parameters'}->{'ET'} )
	{
		RunCom( "$bin/flag_anchored_elements.py $final $cfg->{'Parameters'}->{'ET'} > $final.tmp" );
		RunCom( "mv $final.tmp $final" );
	}

}
# ------------------------------------------------
sub RunHmmPBS
{
	my( $mod, $path_output, $path_input, $opt, $evi ) = @_;
	print "running gm.hmm on PBS\n" if $v;
	
	$path_input = "$work_dir/data/training" if !$path_input;
	
	$opt = ' -f tr ' if !$opt;
	
	if( $evi )
	{
		if( ! -e $evi  ) { print "error, evidence file not found\n"; exit 1;}
		
		$opt .= ' -d '. $evi .' ';
	}

	if ( $cfg->{'Parameters'}->{'mask_penalty'} )
	{
		$opt .= " -k $cfg->{'Parameters'}->{'mask_penalty'} ";
	}
	
	my $hmm_par = "\" -m $mod  $opt \"";
	
	my $com = "$bin/run_hmm_pbs.pl --par $hmm_par --out $path_output  --in $path_input  --w ";
	
	writeToLog($com);
	RunComWithWarning( $com );
}
# ------------------------------------------------
sub RunComWithWarning
{
	my ( $com, $mess ) = @_;
	$mess = $com if !defined $mess;
	my $res = system( $com );
	if( $res )
	{
		print "warning on: $mess\n";
		writeToLog("warning on: $mess");
	}

	writeToLog("warning on: $mess") if $debug;
}
# ------------------------------------------------
sub RunHmmLocal
{
	my( $mod, $path_output, $path_input, $opt, $evi ) = @_;
	print "running gm.hmm on local system\n" if $v;
	
	my $gm_hmm =  $cfg->{'Config'}->{'gm_hmm'};
	
	$path_input = "$work_dir/data/training" if !$path_input;
	
	my @contigs = ReadContigNames( abs_path( $path_input ) );
	
	$opt = ' -f tr ' if !$opt;
	$opt .= " -T $Bin " if $key_bin;
	
	if( $evi )
	{
		if( ! -e $evi  ) { print "error, evidence file not found\n"; exit 1;}
		
		$opt .= ' -d '. $evi .' ';
	}

	if ( $cfg->{'Parameters'}->{'mask_penalty'} )
	{
		$opt .= " -k $cfg->{'Parameters'}->{'mask_penalty'} ";
	}

	foreach my $file ( @contigs )
	{
		print "$file\n" if $debug;
		
		my $out = File::Spec->splitpath( $file );
		$out = File::Spec->catfile( $path_output, $out .".out" );
		
		RunComWithWarning( "$gm_hmm  -m $mod  $opt  -o $out  $file");
	}
}
# ------------------------------------------------
sub RunHmmOnCores
{
	my( $mod, $path_output, $cores, $path_input, $opt, $evi ) = @_;
	print "running gm.hmm on local multi-core system\n" if $v;
	
	my $gm_hmm =  $cfg->{'Config'}->{'gm_hmm'};
	
	$path_input = "$work_dir/data/training" if !$path_input ;	
	
	my @contigs = ReadContigNames( abs_path( $path_input ) );

	$opt = ' -f tr ' if !$opt;
	$opt .= " -T $Bin " if $key_bin;
	
	if( $evi )
	{
		if( ! -e $evi  ) { print "error, evidence file not found\n"; exit 1;}
		
		$opt .= ' -d '. $evi .' ';
	}

	if ($cfg->{'Parameters'}->{'mask_penalty'})
	{
		$opt .= " -k $cfg->{'Parameters'}->{'mask_penalty'} ";
	}

	my $manager = new Parallel::ForkManager( $cores );
	
	for my $file (@contigs)
	{
		print "$file\n" if $debug;
		
		my $out = File::Spec->splitpath( $file );
		$out = File::Spec->catfile( $path_output, $out .".out" );

		$manager->start and next;
		RunComWithWarning( "$gm_hmm  -m $mod  $opt  -o $out  $file");
		$manager->finish;
	}
	
	$manager->wait_all_children;
}
# ------------------------------------------------
sub RunHmm
{
	my( $mod, $path_output, $path_input, $options, $evi ) = @_;

	if( $cfg->{'Parameters'}->{'pbs'} )
	{
		RunHmmPBS( $mod, $path_output, $path_input, $options, $evi );
	}
	else
	{
		if( $cfg->{'Parameters'}->{'cores'} > 1 )
		{
			RunHmmOnCores( $mod, $path_output, $cfg->{'Parameters'}->{'cores'}, $path_input, $options, $evi );
		}
		else
		{
			RunHmmLocal( $mod, $path_output, $path_input, $options, $evi );
		}
	}
}
# ------------------------------------------------
sub ReadContigNames
{
	my $dir = shift;
	my @list;

	my @list_with_id = ();

	opendir( DIR, $dir ) or die "error on open directory $0: $dir, $!\n";
	foreach my $file ( grep{ /dna.fa_\d+$/ } readdir(DIR) )
	{
		if ( $file =~ /^dna.fa_(\d+)$/ )
		{
			push @list_with_id, $1;
		}
	}
	closedir DIR;

	foreach my $file_id (sort{$a<=>$b} @list_with_id)
	{
		my $file =  File::Spec->catfile( $dir, "dna.fa_". $file_id );
		if ( -f $file )
		{
			$file = abs_path($file);
			push @list, $file;
		}
		else { print "error, unexpected name found $0: $file\n"; exit 1; }
	}

	my $message = scalar @list ." contigs in list";
	writeToLog($message);
 	print "$message\n" if $v;
 	
 	return @list;
}
# ------------------------------------------------
sub GetHeuristicFileName
{
	my( $GC ) = @_;
	$GC = int $GC;

	my $MIN_HEURISTIC_GC = 32;
	my $MAX_HEURISTIC_GC = 70;
	if( $GC < $MIN_HEURISTIC_GC ) { $GC = $MIN_HEURISTIC_GC; }
	if( $GC > $MAX_HEURISTIC_GC ) { $GC = $MAX_HEURISTIC_GC; }
	
	return ResolvePath( "heu_05_gcode_1_gc_$GC.mod", $cfg->{'Config'}->{'heu_dir'} );
}
# ------------------------------------------------
sub GetGCfromStatfile
{
	print "get GC of sequence\n" if $v;
	
	my $name = shift;
	my $GC = 0;
	
	open( my $IN, $name ) or die "error on open file $0: $name, $!\n";
	while(<$IN>)
	{
		if( /^GC\s+(\S+)\s*$/ )
		{
			$GC = $1;
			last;
		}
	}
	close $IN;
	
	if ( !$GC ) { print "error, GC of sequence is zero $0\n"; exit 1; }
	
	return int( $GC + 0.5 );
}
# ------------------------------------------------
sub PrepareInitialModel
{
	my $function = shift;

	print "prepare initial model\n" if $v;
	chdir $work_dir;
	
	my $ini_mod = $cfg->{'Parameters'}->{'ini_mod'};
	
	if( !$ini_mod )
	{
		# choose one of heuristic models
		# default method in ES and prerequisite for ET/EP
		
		my $GC = GetGCfromStatfile( "info/training.general" );
		$ini_mod = GetHeuristicFileName( $GC );
		
		# save info about heuristic model name here
		$cfg->{'Parameters'}->{'ini_mod'} = $ini_mod;
		
		print "GC $GC\n" if $v;

		$ini_mod = $function->( $ini_mod );
	}
	
	if ( ! $ini_mod ) { print "error, initiation model file not specified $0\n"; exit 1; }
	if ( ! -e $ini_mod ) { print "error, initiation model file not found $0: $ini_mod\n"; exit 1; }
	
	chdir $work_dir;
	RunCom( "ln -sf  $ini_mod  run/ini.mod" );
}
# ------------------------------------------------
sub TrainingDataReport
{
	print "training data report\n" if $v;
	chdir $work_dir;
	if ( ! -e "data/training.fna" ) { print "error, file not found: info/training.fna\n"; exit 1; }
	RunCom( "$bin/probuild --seq data/training.fna --stat info/training.general --allow_x  --GC_PRECISION 0  --details");
}
# ------------------------------------------------
sub SplitFasta
{
	my( $file, $dir, $min_contig, $trace_name ) = @_;

	# this function splits sequence for training or prediction step
	# difference is in the length of "ignored" short contigs
	# in training:   $min_contig from parameter section
	# in prediction: $min_contig_in_predict
	
	chdir $work_dir;
	
	# remove old files if any
	opendir( DIR, $dir ) or die "error on open directory: $dir, $!\n";
	foreach my $file ( grep{ /dna.fa_\d+$/ } readdir(DIR) )
	{
		unlink "$dir/$file" or die "Could not unlink $file: $!";
	}
	closedir DIR;

	my $max_contig = $cfg->{'Parameters'}->{'max_contig'};
	my $max_gap    = $cfg->{'Parameters'}->{'max_gap'};
	my $max_mask   = $cfg->{'Parameters'}->{'max_mask'};
	
	chdir $dir;
	
	RunCom( "$bin/probuild  --seq $file  --split dna.fa  --max_contig $max_contig --min_contig $min_contig --letters_per_line 100 --split_at_n $max_gap --split_at_x $max_mask --allow_x --x_to_n  --trace ../../info/$trace_name " );
	chdir $work_dir;
}
# ------------------------------------------------
sub CommitTrainingData
{
	print "commit training data\n" if $v;
	
	chdir $work_dir;

	CommitTrainingSequence();
	CheckTrainingSetSize();

	CommitTrainingEvidence() if $cfg->{'Parameters'}->{'evidence'};

	CommitTrainingET()  if $cfg->{'Parameters'}->{'ET'};
	CommitTrainingEP()  if $cfg->{'Parameters'}->{'EP'};
	CommitTrainingETP() if $cfg->{'Parameters'}->{'ETP'};
}
# ------------------------------------------------
sub CheckTrainingSetSize
{
	chdir $work_dir;

	if ( ! -e "data/training.fna" )
	{
		die "error, file not found: data/training.fna\n"; 
	}

	my $size = -s "data/training.fna";

	if ( $size < 1000000 )
	{
		die "error, not enough sequnce to run training in file $work_dir/data/training.fna\n";
	}
}
# ------------------------------------------------
sub CommitTrainingSequence
{	
	chdir $work_dir;
	my $dir = "data/training";
	
	# remove old files if any
	if ( -e "data/training.fna" )
	{
		unlink "data/training.fna" or die "Could not unlink: $!";
	}
	
	opendir( DIR, $dir );
	foreach my $file ( grep{ /dna.fa_\d+$/ } readdir(DIR) )
	{
		 unlink "$dir/$file" or die "Could not unlink: $!";
	}
	
	# create new files
	SplitFasta( abs_path( "data/dna.fna" ), $dir, $cfg->{'Parameters'}->{'min_contig'}, "training.trace" );	
	
	chdir $work_dir;
	
	opendir( DIR, $dir );
	foreach  my $file ( grep{ /dna.fa_\d+$/ } readdir(DIR) )
	{
		my $com = "cat $dir/$file >> data/training.fna";
		system( $com ) && die "$!\n";
	}
	closedir DIR;
}
# ------------------------------------------------
sub CommitTrainingEvidence
{
	RunCom( "$bin/rescale_gff.pl  --in data/evidence.gff  --trace info/training.trace  --out data/evidence_training.gff");
}
# ------------------------------------------------
sub CommitTrainingET
{
	RunCom( "$bin/rescale_gff.pl  --in data/et.gff  --trace info/training.trace  --out data/et_training.gff" );
}
# ------------------------------------------------
sub CommitTrainingEP
{
	if ( $cfg->{'Parameters'}->{'EP'} =~ /\S/ )
	{
		RunCom( "$bin/rescale_gff.pl  --in data/ep.gff  --trace info/training.trace  --out data/ep_training.gff" );
	}
}
# ------------------------------------------------
sub CommitTrainingETP
{
	;
}
# ------------------------------------------------
sub ResolveAutoSoftMask
{
	if ( $cfg->{'Parameters'}->{'soft_mask'} eq "auto" )
	{
		my $gsize = -s $cfg->{'Parameters'}->{'sequence'};

		if ( $gsize < 300000000 )
		{
			$cfg->{'Parameters'}->{'soft_mask'} = 1000;
		}
		else
		{
			$cfg->{'Parameters'}->{'soft_mask'} = 100;
		}

		print "soft_mask is in the 'auto' mode. soft_mask was set to: $cfg->{'Parameters'}->{'soft_mask'}\n" if $v;
	}
}
# ------------------------------------------------
sub DataReport
{
	print "data report\n" if $v;
	chdir $work_dir;
	
	RunCom( "$bin/probuild  --seq data/dna.fna  --allow_x  --stat info/dna.general  --details" );
	RunCom( "$bin/probuild  --seq data/dna.fna  --allow_x  --stat_fasta info/dna.multi_fasta" );
	RunCom( "$bin/probuild  --seq data/dna.fna  --allow_x  --substring_n_distr info/dna.gap_distr" );
	RunCom( "$bin/gc_distr.pl --in data/dna.fna  --out info/dna.gc.csv  --w 1000,8000" );
};
# ------------------------------------------------
sub CommitData
{
	print "commit input data\n" if $v;
	chdir $work_dir;

	CommitSequnece() if $cfg->{'Parameters'}->{'sequence'};
	CommitEvidence() if $cfg->{'Parameters'}->{'evidence'};

	CommitDataET()  if $cfg->{'Parameters'}->{'ET'};
	CommitDataEP()  if $cfg->{'Parameters'}->{'EP'};
	CommitDataETP() if $cfg->{'Parameters'}->{'ETP'};
}
# ------------------------------------------------
sub CommitSequnece
{
	# replace original FASTA defline
	# keep old/new relationship in trace file
	# check for valid alphabet
	# if softmask is defined, then hardmask lower case letters according the softmask value
	# uppercase
	
	my $name = $cfg->{'Parameters'}->{'sequence'};

	ResolveAutoSoftMask();
	my $str = '';
	$str = ('--mask_soft '. $cfg->{'Parameters'}->{'soft_mask'} .' ' ) if $cfg->{'Parameters'}->{'soft_mask'};

	RunCom("$bin/probuild --reformat_fasta --uppercase --allow_x --letters_per_line 60 --out data/dna.fna --label _dna --trace info/dna.trace --in $name  $str" );
}
# ------------------------------------------------
sub CommitEvidence
{
	# synchronize values in sequence name column with new FASTA defline

	my $name = $cfg->{'Parameters'}->{'evidence'};
	my $str = '';
	$str = '--quiet' if !$debug;
	RunCom( "$bin/reformat_gff.pl --out data/evidence.gff  --trace info/dna.trace  --in $name  $str" );
}
# ------------------------------------------------
sub CommitDataET
{
	my $name = $cfg->{'Parameters'}->{'ET'};
	my $str = '';
	$str = '--quiet' if !$debug;
	RunCom( "$bin/reformat_gff.pl --out data/et.gff  --trace info/dna.trace  --in $name  $str" );
}
# ------------------------------------------------
sub CommitDataEP
{
	my $name = $cfg->{'Parameters'}->{'EP'};

	if ( $name =~ /\S/ )
	{
		my $str = '';
		$str = '--quiet' if !$debug;
		RunCom( "$bin/reformat_gff.pl --out data/ep.gff  --trace info/dna.trace  --in $name  $str" );
	}
}
# ------------------------------------------------
sub CommitDataETP
{
	my $name = $cfg->{'Parameters'}->{'EP'};
	my $str = '';
	$str = '--quiet' if !$debug;
	RunCom( "$bin/reformat_gff.pl --out data/etp.gff  --trace info/dna.trace  --in $name  $str" );
}
# ------------------------------------------------
sub RunCom
{
	my ( $com, $mess ) = @_;
	my $res = system( $com );
	$mess = $com if !defined $mess;
	writeToLog($mess);
	if( $res ) { print "error on call: $mess\n"; writeToLog("error"); exit 1; }
}
# ------------------------------------------------
sub CreateDirectories
{
	print "create directories\n" if $v;
	chdir $work_dir;
	
	my @list =
	(
		'data',
		'info',
		'data/training',
		'run',
		'output'
	);

	make_path( @list,{ verbose => $debug } );
}
# ------------------------------------------------
sub ReadParameters
{
	SetDefaultValues();
	ReadCfgFile( $cfg->{'Config'}->{'def_cfg'}, $bin );
	$cfg->{'Config'}->{'gm_hmm'} = ResolvePath( $cfg->{'Config'}->{'gm_hmm'}, $bin );
	$cfg->{'Config'}->{'version'} .= "_lic" if (GetHMMVersion() =~ /lic/);
#       $cfg->{'Config'}->{'version'} .= GetHMMVersion();
#	$cfg->{'Config'}->{'version'} .= GetProtHintVersion();
	# Usage() prints some parameters from default configuration on screen
	# Read configuration file before outputting the Usage()
	Usage() if ( @ARGV < 1 );
	ParseCMD();
	CheckBeforeRun();
	SetLogger();
	SaveRunningConfiguration();

	print Dumper($cfg) if $debug;
}
# ------------------------------------------------
sub SaveRunningConfiguration
{
	my $name = $cfg->{'Config'}->{'run_cfg'};
	if ( !$name ) { print "error, 'run_cfg' configuration file name is missing: $0\n"; exit 1; }
	$name = abs_path( File::Spec->catfile( $work_dir, $name ) );
	$cfg->{'Config'}->{'run_cfg'} = $name;
	open( my $OUT, ">$name") or die "error on open file $0: $name\n$!\n";
	print $OUT Dump($cfg);
	close $OUT;
}
# ------------------------------------------------
sub SetLogger
{
	my $name = $cfg->{'Config'}->{'log_file'};
	if ( !$name ) { print "error, 'log_file' file name is missing: $0\n"; exit 1; }
	$name = abs_path( File::Spec->catfile( $work_dir, $name ) );
	$cfg->{'Config'}->{'log_file'} = $name;

	open($log, '>', $name) or die "error on creating file $name: $!\n";
	$logMutex = MCE::Mutex->new;
	writeToLog(Dumper($cfg)) if $debug;
}
# ------------------------------------------------
sub writeToLog
{
	my $logText = shift;
	$logMutex->lock;
	print $log "$bin/gmes_petap.pl : [" . localtime() . "] " . $logText . "\n";
	$logMutex->unlock;
}
# ------------------------------------------------
sub CheckOutOfRangeInt
{
	my $key = shift;
	my $min = shift;
	my $max = shift;

	if( ! exists $cfg->{'Parameters'}->{$key} )
		{ print "error, label is not found in parameter set: $key\n"; exit 1; }

	if ( $cfg->{'Parameters'}->{$key} !~ /^\s*([+-]?\d+\.?\d*)\s*$/ )
		{ print "error, numerical value is expected for $key: $cfg->{'Parameters'}->{$key}\n"; exit 1; }

	if ( defined $min )
	{
		if ( $cfg->{'Parameters'}->{$key} < $min )
			{ print "error, out of range value was specified for $key: $cfg->{'Parameters'}->{$key}\n"; exit 1; }
	}

	if ( defined $max )
	{
		if ( $cfg->{'Parameters'}->{$key} > $max )
			{ print "error, out of range value was specified for $key: $cfg->{'Parameters'}->{$key}\n"; exit 1; }
	}
}
# ------------------------------------------------
sub CheckBeforeRun
{
	print "check before run\n" if $v;
	
	# check sequence pre-processing parameters
	CheckOutOfRangeInt( 'max_contig', 100000 );
	CheckOutOfRangeInt( 'min_contig', 0 );
	CheckOutOfRangeInt( 'max_gap', 0 );
	CheckOutOfRangeInt( 'max_mask', 0 );

	if ( $cfg->{'Parameters'}->{'max_contig'} <  $cfg->{'Parameters'}->{'min_contig'} ) 
		{ print "error, minimum contig length is more than maximum:$cfg->{'Parameters'}->{'min_contig'}  $cfg->{'Parameters'}->{'max_contig'}\n"; exit 1; }

	CheckOutOfRangeInt( 'min_contig_in_predict', 0 );
	CheckOutOfRangeInt( 'min_gene_in_predict', 0 );

	CheckOutOfRangeInt( 'max_intron', 0 );
	CheckOutOfRangeInt( 'max_intergenic', 0 );

	if (  $cfg->{'Parameters'}->{'soft_mask'} ne "auto" )
	{
		CheckOutOfRangeInt( 'soft_mask', 0 );
	}

	if( $cfg->{'Parameters'}->{'cores'} < 1 or $cfg->{'Parameters'}->{'cores'} > 128 )
		{ print "error, out of range values specified for cores: $0\n"; exit 1; }

	# move to abs path
	$bin = ResolvePath( $bin );
	$cfg->{'Config'}->{'bin'} = $bin;
	$work_dir = ResolvePath( $work_dir );
	$cfg->{'Config'}->{'work_dir'} = $work_dir;
	
	if ( $bin eq $work_dir )
		{ print "error, code cannot be executed in installation folder\n"; exit 1; }
	
	$cfg->{'Config'}->{'heu_dir'} = ResolvePath( $cfg->{'Config'}->{'heu_dir'}, $bin );
	
	$cfg->{'Parameters'}->{'sequence'}     = ResolvePath( $cfg->{'Parameters'}->{'sequence'} );
	$cfg->{'Parameters'}->{'ET'}           = ResolvePath( $cfg->{'Parameters'}->{'ET'} );
	$cfg->{'Parameters'}->{'predict_with'} = ResolvePath( $cfg->{'Parameters'}->{'predict_with'} );
	$cfg->{'Parameters'}->{'evidence'}     = ResolvePath( $cfg->{'Parameters'}->{'evidence'} );

	if ( $cfg->{'Parameters'}->{'EP'} )
	{
		$cfg->{'Parameters'}->{'EP'} = ResolvePath( $cfg->{'Parameters'}->{'EP'} );
	}
	else
	{
		$cfg->{'Parameters'}->{'EP'} = " ";
		$cfg->{'Parameters'}->{'dbep'} = ResolvePath( $cfg->{'Parameters'}->{'dbep'} );
	}

	$cfg->{'Parameters'}->{'ini_mod'}      = ResolvePath( $cfg->{'Parameters'}->{'ini_mod'} );
	$cfg->{'Parameters'}->{'test_set'}     = ResolvePath( $cfg->{'Parameters'}->{'test_set'} );
	$cfg->{'Parameters'}->{'usr_cfg'}      = ResolvePath( $cfg->{'Parameters'}->{'usr_cfg'} );
	
	# check input sequence file 
	if( $cfg->{'Run'}->{'commit_input_data'} )
	{
	 	if( !$cfg->{'Parameters'}->{'sequence'} or !-e $cfg->{'Parameters'}->{'sequence'} )
			{ print "error, file with input sequence not found: $0\n"; exit 1; }
		
		if( !-f $cfg->{'Parameters'}->{'sequence'} )
			{ print "error, input not a file $0: $cfg->{'Parameters'}->{'sequence'}\n"; exit 1; }
	}
	
	# check training mode
	if( !$cfg->{'Parameters'}->{'ES'} and !$cfg->{'Parameters'}->{'ET'} and !$cfg->{'Parameters'}->{'EP'} and !$cfg->{'Parameters'}->{'predict_with'} and !$cfg->{'Parameters'}->{'ETP'} )
		{ print "error, only one of hmm, ES, ET, EP or ETP should be specified: $0\n"; exit 1; }

	if( $cfg->{'Parameters'}->{'training_only'} and $cfg->{'Parameters'}->{'prediction_only'} )
		{ print "error, only one of two command line parameters should be specified: training_only or prediction_only\n"; exit 1; }

	if ( $cfg->{'Parameters'}->{'ETP'} and !$cfg->{'Parameters'}->{'ini_mod'} )
	{
		print "error, initiation model file was not specified in ETP mode\n";
		print "# model file generated by ET is recommended as initial model\n";
		exit 1;
	}

	# misc
	if( ! $key_bin )
	{
		my $file_name = glob('~/.gm_key');
		$key_bin = 1 if ( ! -e $file_name );
		print "test key: $key_bin\n" if $debug;
	}
	
	CheckOutOfRangeInt( 'et_score', 0 );
	CheckEPscore( $cfg->{'Parameters'}->{'ep_score'} ); 

	CheckGCdonorOption();
}
# ------------------------------------------------
sub CheckEPscore
{
	my $value = shift;

	if ( $value =~ /\s*(\d+),([+-]?\d+\.?\d*)\s*/ )
	{
		;
	}
	else
	{
		CheckOutOfRangeInt( 'ep_score', 0 ) 
	}
}
# ------------------------------------------------
sub CheckGCdonorOption
{
	if ( $cfg->{'Parameters'}->{'gc_donor'} =~ /^\s*([+-]?\d+\.?\d*)\s*$/ )
	{
		if ($1 < 0 or $1 > 1)
			{ print "error, parameter value for --gc_donor is out of the allowed range 0..1: $1\n"; exit 1; }
		
		if ( $cfg->{'Parameters'}->{'gc_donor'} == 0 )
		{
			$cfg->{'donor_GC'}->{'status'} = 0;
			$cfg->{'donor_GT'}->{'status'} = 1;
		}
		elsif ( $cfg->{'Parameters'}->{'gc_donor'} == 1 )
		{
			$cfg->{'donor_GC'}->{'status'} = 1;
			$cfg->{'donor_GT'}->{'status'} = 0;
		}
	}
	elsif ( $cfg->{'Parameters'}->{'gc_donor'} eq "off" )
	{
		$cfg->{'donor_GC'}->{'status'} = 0;
		$cfg->{'donor_GT'}->{'status'} = 1;

		$cfg->{'Parameters'}->{'gc_donor'} = 0;
	}
	elsif ( $cfg->{'Parameters'}->{'gc_donor'} eq "auto" )
		{ print "error, option is not supported yet\n"; exit 1; }
	else
		{ print "error, unexpected value for parameter --gc_donor was found: $cfg->{'Parameters'}->{'gc_donor'}\n"; exit 1; }
}
# ------------------------------------------------
sub ParseCMD
{
	my $cmd = $0;
	foreach my $str (@ARGV) { $cmd .= ( ' '. $str ); } 
	
	my %h;
	my $opt_results = GetOptions
	(
		\%h,
		
		'sequence=s',
		
		'ES',
		'ET=s',
		'et_score=f',
		'EP:s{,1}',
		'ep_score=s',
		'dbep=s',
		'predict_with=s',
		'ETP=s',

		'fungus',
		'evidence=s',
		
		'cores=i',
		'pbs',
		
		'max_contig=i',
		'min_contig=i',
		'max_gap=i',
		'max_mask=i',
		'soft_mask=s',
		
		'max_intron=i',
		'max_intergenic=i',
		'min_contig_in_predict=i',
		'min_gene_in_predict=i',
		'gc_donor=s',

		'gc3=f',
		'training_only',
		'prediction_only',
		'usr_cfg=s',
		'ini_mod=s',
		'test_set=s',
		'mask_penalty=f',
		
		'key_bin' => \$key_bin,
		'verbose' => \$v,
		'debug'   => \$debug,
	);

	if( !$opt_results ) { print "error on command line: $0\n"; exit 1; }
	if( @ARGV > 0 ) { print "error, unexpected argument found on command line: $0 @ARGV\n"; exit 1; }
	
	$v = 1 if $debug;
	
	# user may specify additional configuration file on command line
	# parse user specified file first (if any) and then parse other command line parameters

	ReadCfgFile( $h{'usr_cfg'} ) if exists $h{'usr_cfg'};
	
	# fungi section: some default fungi parameters differ from default euk parameters
	# overwrite default euk parameters by fungi specific
	# do this before processing all other command lines

	if ( exists $h{'fungus'} )
	{
		$cfg->{'intron_DUR'}->{'max'}    = $cfg->{'Fungi'}->{'max_intron'};
		$cfg->{'prespacer_DUR'}->{'max'} = $cfg->{'Fungi'}->{'max_intron'};
		$cfg->{'Parameters'}->{'min_gene_in_predict'} = $cfg->{'Fungi'}->{'min_gene_in_predict'};
		$cfg->{'intergenic_DUR'}->{'max'} = $cfg->{'Fungi'}->{'max_intergenic'};
		$cfg->{'Parameters'}->{'min_contig'} = $cfg->{'Fungi'}->{'min_contig'};
	}
	
	# transfer command line parameters to cfg
	foreach my $key ( keys %h )
	{
		$cfg->{'Parameters'}->{$key} = $h{$key};
	}
	
	UpdateRunStatus();
	
	# save informaton for debug
	$cfg->{'Parameters'}->{'cmd'}     = $cmd;
	$cfg->{'Parameters'}->{'v'}       = $v;
	$cfg->{'Parameters'}->{'debug'}   = $debug;
	$cfg->{'Parameters'}->{'key_bin'} = $key_bin;

	# some cpmmand line parameters belong also to sections other than 'Parameters'
	# copy such values to corresponding sections

	# set score for initial ET/EP model
	if ( $cfg->{'Parameters'}->{'ET'} )
	{
		$cfg->{'ET_ini'}->{'score'} = $cfg->{'Parameters'}->{'et_score'};
	}
	
	if ( $cfg->{'Parameters'}->{'EP'} )
	{
		$cfg->{'EP_ini'}->{'score'} = $cfg->{'Parameters'}->{'ep_score'};
	}

	# if custom intron length was specified:
	if( $cfg->{'Parameters'}->{'max_intron'} > 0 )
	{
		$cfg->{'intron_DUR'}->{'max'}    = $cfg->{'Parameters'}->{'max_intron'};
		$cfg->{'prespacer_DUR'}->{'max'} = $cfg->{'intron_DUR'}->{'max'};
	}
	
	# if custom intergenic length was specified:
	if( $cfg->{'Parameters'}->{'max_intergenic'} > 0 )
	{
		$cfg->{'intergenic_DUR'}->{'max'} = $cfg->{'Parameters'}->{'max_intergenic'};
	}
}
# ------------------------------------------------
sub UpdateRunStatus
{
	if( $cfg->{'Parameters'}->{'predict_with'} )
	{
		$cfg->{'Run'}->{'commit_training_data'} = 0;
		$cfg->{'Run'}->{'training_data_report'} = 0;
		$cfg->{'Run'}->{'prepare_ini_mod'}      = 0;
		$cfg->{'Run'}->{'run_training'}         = 0;
		$cfg->{'Run'}->{'training_report'}      = 0;
	}
	elsif( $cfg->{'Parameters'}->{'training_only'} and !$cfg->{'Parameters'}->{'prediction_only'} )
	{
		$cfg->{'Run'}->{'run_prediction'}       = 0;
		$cfg->{'Run'}->{'prediction_report'}    = 0;
	}
	elsif( !$cfg->{'Parameters'}->{'training_only'} and $cfg->{'Parameters'}->{'prediction_only'} )
	{
		$cfg->{'Run'}->{'set_dirs'}             = 0;
		$cfg->{'Run'}->{'commit_input_data'}    = 0;
		$cfg->{'Run'}->{'input_data_report'}    = 0;
		$cfg->{'Run'}->{'commit_training_data'} = 0;
		$cfg->{'Run'}->{'training_data_report'} = 0;
		$cfg->{'Run'}->{'prepare_ini_mod'}      = 0;
		$cfg->{'Run'}->{'run_training'}         = 0;
		$cfg->{'Run'}->{'training_report'}      = 0;
	}
}
# ------------------------------------------------
sub ReadCfgFile
{
	my( $name, $path ) = @_;
	return '' if !$name;
	print "reading configuration file: $name\n" if $v;
	$name = ResolvePath( $name, $path );
	Hash::Merge::set_behavior( 'RIGHT_PRECEDENT' );
	
	my $cfg_from_file = YAML::LoadFile( $name );
	
	if( defined $cfg_from_file )
	{
		%$cfg = %{ merge( $cfg, $cfg_from_file) };
	}
	else { print "warning, configuration file is empty: $name\n"; }
	
	return $name;
}
# ------------------------------------------------
sub ResolvePath
{
	my( $name, $path ) = @_;
	return '' if !$name;
	$name = File::Spec->catfile( $path, $name ) if ( defined $path and $path );
	if( ! -e $name ) { print "error, file not found $0: $name\n"; exit 1; }
	return abs_path( $name );
}
# ------------------------------------------------
# default values may be overwritten by values from config file
# ------------------------------------------------
sub SetDefaultValues
{
	# basic configuration
	$cfg->{'Config'}->{'version'}  = "4.65";
	$cfg->{'Config'}->{'heu_dir'}  = "heu_dir";
	$cfg->{'Config'}->{'def_cfg'}  = "gmes.cfg";
	$cfg->{'Config'}->{'gm_hmm'}   = "gmhmme3";
	$cfg->{'Config'}->{'log_file'} = "gmes.log";
	$cfg->{'Config'}->{'run_cfg'}  = "run.cfg";
	$cfg->{'Config'}->{'bin'}      = $bin;
	$cfg->{'Config'}->{'work_dir'} = $work_dir;

	# switch algorithm steps "ON" and "OFF"
	 # data setup steps
	$cfg->{'Run'}->{'set_dirs'}             = 1;  # bool true
	$cfg->{'Run'}->{'commit_input_data'}    = 1;  # bool true
	$cfg->{'Run'}->{'input_data_report'}    = 1;  # bool true
	$cfg->{'Run'}->{'commit_training_data'} = 1;  # bool true
	 # training steps
	$cfg->{'Run'}->{'training_data_report'} = 1;  # bool true
	$cfg->{'Run'}->{'prepare_ini_mod'}      = 1;  # bool true
	$cfg->{'Run'}->{'run_training'}         = 1;  # bool true
	$cfg->{'Run'}->{'training_report'}      = 0;  # false as it is not implemented yet
	 # prediction steps
	$cfg->{'Run'}->{'run_prediction'}       = 1;  # bool true
	$cfg->{'Run'}->{'prediction_report'}    = 0;  # false as it is not implemented yet

	# command line parameters are in 'Parameters' section

	# optional algorithm run modes
	$cfg->{'Parameters'}->{'training_only'}   = 0;  # bool false
	$cfg->{'Parameters'}->{'prediction_only'} = 0;  # bool false
	$cfg->{'Parameters'}->{'ini_mod'} = '';
	$cfg->{'Parameters'}->{'usr_cfg'} = '';
	$cfg->{'Parameters'}->{'test_set'} = '';  # false as it is not implemented yet

	# genome sequence processing
	$cfg->{'Parameters'}->{'max_contig'} = 5000000;
	$cfg->{'Parameters'}->{'min_contig'} =   50000;
	$cfg->{'Parameters'}->{'max_gap'}    =    5000;
	$cfg->{'Parameters'}->{'max_mask'}   =    5000;
	$cfg->{'Parameters'}->{'soft_mask'}  =  "auto";
	$cfg->{'Parameters'}->{'min_contig_in_predict'} = 500;
	$cfg->{'Parameters'}->{'min_gene_in_predict'}   = 300;

	# algorithm mode selection and input data
	$cfg->{'Parameters'}->{'sequence'} = '';
	
	$cfg->{'Parameters'}->{'ES'}       = 0;
	$cfg->{'Parameters'}->{'ET'}       = '';
	$cfg->{'Parameters'}->{'et_score'} = 10;
	$cfg->{'Parameters'}->{'EP'}       = '';
	$cfg->{'Parameters'}->{'ep_score'} = "4,0.25";
	$cfg->{'Parameters'}->{'dbep'}     = '';
	$cfg->{'Parameters'}->{'ETP'}      = '';
	$cfg->{'Parameters'}->{'predict_with'} = '';
	
	$cfg->{'Parameters'}->{'fungus'}   = 0;
	$cfg->{'Parameters'}->{'evidence'} = '';

	# code execution mode
	$cfg->{'Parameters'}->{'pbs'}   = 0;  # bool false
	$cfg->{'Parameters'}->{'cores'} = 1;

	# species parameters
	$cfg->{'Parameters'}->{'max_intergenic'} = 0;  # if 0 then use default value set in intergenic duration
	$cfg->{'Parameters'}->{'max_intron'}     = 0;  # if 0 then use default value set in intron duration
	
	$cfg->{'Parameters'}->{'gc_donor'}       = 0.001;
	$cfg->{'Parameters'}->{'gc3'} = 100; # with 100% value GC3 is switched off
	$cfg->{'Parameters'}->{'mask_penalty'}   = 0;

	# misc
	$cfg->{'Parameters'}->{'v'}       = $v;
	$cfg->{'Parameters'}->{'debug'}   = $debug;
	$cfg->{'Parameters'}->{'key_bin'} = $key_bin;
}
# ------------------------------------------------
sub GetHMMVersion
{
	my $hmm_version = "";
	my $hmm_usage = "";

	if ( -e $cfg->{'Config'}->{'gm_hmm'} )
	{
		$hmm_usage = `$cfg->{'Config'}->{'gm_hmm'}; 2>&1`;
	}

	if ( $hmm_usage  =~ /version\s*(\S+)/ )
	{
		$hmm_version = "_". $1;
	}

	return $hmm_version;
}
# ------------------------------------------------
sub GetProtHintVersion
{
	if ( ! -e "$Bin/ProtHint/bin/prothint.py" )
	{
		return '';
	}
	
	my $com = ResolvePath( "$Bin/ProtHint/bin/prothint.py" );

	my $prot_hint_version = `$com --version`;

	$prot_hint_version =~ s/^\S+\s+/_/;
	
	return $prot_hint_version;
}
# ------------------------------------------------
# To run GeneMark-ETP with hints from both protein and transcriptome splice alignments
# --ETP
# ------------------------------------------------
sub Usage
{
	print qq(# -------------------
Usage:  $0  [options]  --sequence [filename]

GeneMark-ES Suite version $cfg->{'Config'}->{'version'}
Suite includes GeneMark.hmm, GeneMark-ES, GeneMark-ET and GeneMark-EP algorithms.

Input sequence/s should be in FASTA format.

Select one of the gene prediction algorithm

To run GeneMark-ES self-training algorithm
  --ES

To run GeneMark-ET with hints from transcriptome splice alignments
  --ET           [filename]; file with intron coordinates from RNA-Seq read splice alignment in GFF format
  --et_score     [number]; default $cfg->{'Parameters'}->{'et_score'}; minimum score of intron in initiation of the ET algorithm

To run GeneMark-EP with hints from protein splice alignments
  --EP           
  --dbep         [filename]; file with protein database in FASTA format
  --ep_score     [number,number]; default $cfg->{'Parameters'}->{'ep_score'}; minimum score of intron in initiation of the EP algorithm
or
  --EP           [filename]; file with intron coordinates from protein splice alignment in GFF format

To run GeneMark.hmm predictions using previously derived model
  --predict_with [filename]; file with species specific gene prediction parameters

To run ES, ET or EP with branch point model. This option is most useful for fungal genomes
  --fungus

To run hmm, ES, ET or EP in PLUS mode (prediction with hints)
  --evidence     [filename]; file with hints in GFF format

Masking option
  --soft_mask    [number] or [auto]; default $cfg->{'Parameters'}->{'soft_mask'}; to indicate that lowercase letters stand for repeats;
                 masks only lowercase repeats longer than specified length
                 In 'auto' mode length is adjusted based on the size of the input genome

Run options
  --cores        [number]; default $cfg->{'Parameters'}->{'cores'}; to run program with multiple threads
  --pbs          to run on cluster with PBS support
  --v            verbose

Optional sequence pre-processing parameters
  --max_contig   [number]; default $cfg->{'Parameters'}->{'max_contig'}; will split input genomic sequence into contigs shorter then max_contig
  --min_contig   [number]; default $cfg->{'Parameters'}->{'min_contig'} ($cfg->{'Fungi'}->{'min_contig'} fungi); will ignore contigs shorter then min_contig in training 
  --max_gap      [number]; default $cfg->{'Parameters'}->{'max_gap'}; will split sequence at gaps longer than max_gap
                 Letters 'n' and 'N' are interpreted as standing within gaps 
  --max_mask     [number]; default $cfg->{'Parameters'}->{'max_mask'}; will split sequence at repeats longer then max_mask
                 Letters 'x' and 'X' are interpreted as results of hard masking of repeats

Optinal algorithm parameters
  --max_intron            [number]; default $cfg->{'intron_DUR'}->{'max'} ($cfg->{'Fungi'}->{'max_intron'} fungi); maximum length of intron
  --max_intergenic        [number]; default $cfg->{'intergenic_DUR'}->{'max'} ($cfg->{'Fungi'}->{'max_intergenic'} fungi); maximum length of intergenic regions
  --min_contig_in_predict [number]; default $cfg->{'Parameters'}->{'min_contig_in_predict'}; minimum allowed length of contig in prediction step
  --min_gene_in_predict   [number]; default $cfg->{'Parameters'}->{'min_gene_in_predict'} ($cfg->{'Fungi'}->{'min_gene_in_predict'} fungi); minimum allowed gene length in prediction step
  --gc_donor              [value];  default $cfg->{'Parameters'}->{'gc_donor'}; transition probability to GC donor in the range 0..1; 'auto' mode detects probability from training; 'off' switches GC donor model OFF

Developer options
  --gc3          [number]; GC3 cutoff in training for grasses
  --training     to run only training step of algorithms; applicable to ES, ET or EP
  --prediction   to run only prediction step of algorithms using species parameters from previously executed training; applicable to ES, ET or EP
  --usr_cfg      [filename]; use custom configuration from this file
  --ini_mod      [filename]; use this file with parameters for algorithm initiation
  --test_set     [filename]; to evaluate prediction accuracy on the given test set
  --mask_penalty [number]; default $cfg->{'Parameters'}->{'mask_penalty'};
  --key_bin
  --debug
# -------------------
);
	exit 1;
}
# ================== END sub =====================

