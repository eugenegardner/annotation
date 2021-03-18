#!/usr/bin/env perl
# =================================
# Move files from eukaryotic GeneMark.hmm LST format to GTF format.
# Last update 8 2019
#
# Alex Lomsadze
# Georgia Institute of Technology
# alexl@gatech.edu
# =================================

use warnings;
use strict;
use Getopt::Long;
use File::Spec;

use FindBin qw( $RealBin );
use lib $FindBin::Bin ."/lib";
use ReadBackwards;

# ---------------------------------
my $v = 0;
my $debug = 0;

my $infile = '';
my $outfile = '';
my $app = 0;
my $min = 0;
my $format = 'gtf';
# ---------------------------------

Usage() if ( @ARGV < 1 );
ParseCMD();
CheckInput();

# temporary fix for cases of no output from HMM

if ( ! -s $infile )
{
	exit 0;
}

print "#\n" if $v;
print "# input file: $infile\n" if $v;

my $defline = GetDefline( $infile );
my ( $seq_id, $shift ) = ParseDefline( $defline );

my %gene_info = ();
PopulateGeneInfo( $infile, \%gene_info );
FilterGenesByGeneLength( \%gene_info, $min );

# in append mode, get last gene ID from existing output file
my $last_gene_id = GetLastId( $outfile );

my $source = "GeneMark.hmm";
my $type   = "";
my $start  = 0;
my $end    = 0;
my $score  = ".";
my $strand = ".";
my $phase  = ".";
my $attributes = ".";

# tmp to print start/stop codon positions
my $pos;

# id if gene in the output file
my $current_gene_count;

# attributes for different $type
my $att_start = '';
my $att_stop = '';
my $att_exon = '';
my $att_CDS = '';
my $att_intron = '';

# tag for name in $attributes field
my $tag_g;
my $tag_t;

# evidence labels
my $evi_start = 0;
my $evi_end = 0;
my $att_evi = '';

# other helpers
my $gmhmm_gene_id;
my $new_gene_id;
my $cds_id;
my $previous_exon_right = 0;
my $previous_exon_phase = 0;
my $current_gene_id = 0;

open( my $IN, "$infile") || die( "error on open input file $0: $infile, $!\n" );
my $OUT;
if ( $app )
{
        open( $OUT, ">>$outfile") || die( "error on open output file $0: $outfile, $!\n" );
}
else
{
        open( $OUT, ">$outfile") || die( "error on open output file $0: $outfile, $!\n" );
}
while(<$IN>)
{
  # GeneMark.hmm output format
  # gene_id   cds_id   strand   type   start   end   length   start_frame   end_frame   supported(10-11)
  #   1          2        3       4       5      6       7         8            9      10       11
  if( /^\s*(\d+)\s+(\d+)\s+([\+-])\s+(\w+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([-+])\s+([-+])\s*/ )
  {
  	$gmhmm_gene_id = $1;
	if( ! exists $gene_info{$gmhmm_gene_id}{'new_gene_id'} ) {die "error in the code, check gene ID: $gmhmm_gene_id\n";}

	$new_gene_id = $gene_info{$gmhmm_gene_id}{'new_gene_id'};
	next if( ! $new_gene_id );

	$cds_id = $2;
	$strand = $3;
	$start  = $5 + $shift - 1;
	$end    = $6 + $shift - 1;

	if ( $strand eq '+' )
	{
		$phase = (3 - $8 + 1)%3;
	}
	else
	{
		$phase = (3 - $9 + 1)%3;
	}
 
	$current_gene_count = $last_gene_id + $new_gene_id;

	$tag_g =  $current_gene_count ."_g";
	$tag_t =  $current_gene_count ."_t";

	if( $10 eq '+' ) {$evi_start = 1;} else {$evi_start = 0;}
	if( $11 eq '+' ) {$evi_end = 1;}   else {$evi_end = 0;}
	$att_evi = $evi_start ."_". $evi_end;

	$attributes = "gene_id \"". $tag_g ."\"; transcript_id \"". $tag_t ."\";";

	$att_start = $attributes ." count \"1_1\"\;";
	$att_stop  = $attributes ." count \"1_1\"\;";
	$att_exon  = $attributes;
	$att_intron = $attributes;
	$att_CDS   = $attributes;
    
 	if( $att_evi ne '0_0' )
	{
		$att_CDS  = $attributes . " evidence \"". $att_evi ."\";";
	}

	if ( $current_gene_id ne $current_gene_count )
	{
		$current_gene_id = $current_gene_count;
		$previous_exon_right = 0;
	}

	if ( $strand eq "+" )
	{
		if ( $previous_exon_right > 0 )
		{
			my $intron_phase = 3 - $phase;
			$intron_phase = 0 if ( $intron_phase == 3);

			$att_intron .= " count \"". ($cds_id -1) ."_". ($gene_info{$gmhmm_gene_id}{'cds_count'} -1) ."\";";

			print $OUT "$seq_id\t$source\tintron\t". ($previous_exon_right + 1) ."\t". ($start - 1) ."\t0\t$strand\t$intron_phase\t$att_intron\n";
		}
		$previous_exon_right = $end;
		$previous_exon_phase = $phase;

		if( $4 eq "Initial" )
		{
			print $OUT "$seq_id\t$source\texon\t$start\t$end\t0\t$strand\t.\t$att_CDS\n";
			$pos = $start + 2;
			print $OUT "$seq_id\t$source\tstart_codon\t$start\t$pos\t$score\t$strand\t0\t$att_start\n";
			$att_CDS  .= " cds_type \"Initial\"\; count \"". $cds_id ."_". $gene_info{$gmhmm_gene_id}{'cds_count'} ."\";";
			print $OUT "$seq_id\t$source\tCDS\t$start\t$end\t$score\t$strand\t$phase\t$att_CDS\n";

		}
		elsif( $4 eq "Internal")
		{
			print $OUT "$seq_id\t$source\texon\t$start\t$end\t0\t$strand\t.\t$att_CDS\n";
			$att_CDS  .= " cds_type \"Internal\"\; count \"". $cds_id ."_". $gene_info{$gmhmm_gene_id}{'cds_count'} ."\";";
			print $OUT "$seq_id\t$source\tCDS\t$start\t$end\t$score\t$strand\t$phase\t$att_CDS\n";
		}
		elsif( $4 eq "Terminal" )
		{
			print $OUT "$seq_id\t$source\texon\t$start\t$end\t0\t$strand\t.\t$att_CDS\n";
			$att_CDS  .= " cds_type \"Terminal\"\; count \"". $cds_id ."_". $gene_info{$gmhmm_gene_id}{'cds_count'} ."\";";
			if ( $format eq 'gtf2' )
			{
				if ( $end - $start + 1 > 3 )
				{
					print $OUT "$seq_id\t$source\tCDS\t$start\t". ($end -3) ."\t$score\t$strand\t$phase\t$att_CDS\n";
				}
			}
			else
			{
				print $OUT "$seq_id\t$source\tCDS\t$start\t$end\t$score\t$strand\t$phase\t$att_CDS\n";
			}

			$pos = $end - 2;
			print $OUT "$seq_id\t$source\tstop_codon\t$pos\t$end\t$score\t$strand\t0\t$att_stop\n";
		}
		elsif ( $4 eq "Single" )
		{
			print $OUT "$seq_id\t$source\texon\t$start\t$end\t0\t$strand\t.\t$att_CDS\n";
			$pos = $start + 2;
                        print $OUT "$seq_id\t$source\tstart_codon\t$start\t$pos\t$score\t$strand\t0\t$att_start\n";
			$att_CDS  .= " cds_type \"Single\"\; count \"". $cds_id ."_". $gene_info{$gmhmm_gene_id}{'cds_count'} ."\";";
			if ( $format eq 'gtf2' )
                        {
				if ( $end - $start + 1 > 3 )
                        	{
                                	print $OUT "$seq_id\t$source\tCDS\t$start\t". ($end -3) ."\t$score\t$strand\t$phase\t$att_CDS\n";
                        	}
			}
			else
			{
				print $OUT "$seq_id\t$source\tCDS\t$start\t$end\t$score\t$strand\t$phase\t$att_CDS\n";
			}

                        $pos = $end - 2;
                        print $OUT "$seq_id\t$source\tstop_codon\t$pos\t$end\t$score\t$strand\t0\t$att_stop\n";
		}
	}
	elsif ( $strand eq "-" )
	{
		if ( $previous_exon_right > 0 )
		{
			my $intron_phase = 3 - $previous_exon_phase;
			$intron_phase = 0 if ( $intron_phase == 3);

			$att_intron .= " count \"". $cds_id ."_". ($gene_info{$gmhmm_gene_id}{'cds_count'} -1) ."\";";

			print $OUT "$seq_id\t$source\tintron\t". ($previous_exon_right + 1) ."\t". ($start - 1) ."\t0\t$strand\t$intron_phase\t$att_intron\n";
		}
		$previous_exon_right = $end;
		$previous_exon_phase = $phase;

		if( $4 eq "Terminal" )
		{
			print $OUT "$seq_id\t$source\texon\t$start\t$end\t0\t$strand\t.\t$att_CDS\n";
			$pos = $start + 2;
			print $OUT "$seq_id\t$source\tstop_codon\t$start\t$pos\t$score\t$strand\t0\t$att_stop\n";
			$att_CDS .= " cds_type \"Terminal\"\; count \"". $cds_id ."_". $gene_info{$gmhmm_gene_id}{'cds_count'} ."\";";
			if ( $format eq 'gtf2' )
			{
				if( $end - $start + 1 > 3)
				{
					print $OUT "$seq_id\t$source\tCDS\t". ($start +3)."\t$end\t$score\t$strand\t$phase\t$att_CDS\n";
				}
			}
			else
			{
				print $OUT "$seq_id\t$source\tCDS\t$start\t$end\t$score\t$strand\t$phase\t$att_CDS\n";
			}
		}
		elsif( $4 eq "Internal" )
		{
			print $OUT "$seq_id\t$source\texon\t$start\t$end\t0\t$strand\t.\t$att_CDS\n";
			$att_CDS .= " cds_type \"Internal\"\; count \"". $cds_id ."_". $gene_info{$gmhmm_gene_id}{'cds_count'} ."\";";
			print $OUT "$seq_id\t$source\tCDS\t$start\t$end\t$score\t$strand\t$phase\t$att_CDS\n";
		}
		elsif( $4 eq "Initial" )
		{
			print $OUT "$seq_id\t$source\texon\t$start\t$end\t0\t$strand\t.\t$att_CDS\n";
			$att_CDS .= " cds_type \"Initial\"\; count \"". $cds_id ."_". $gene_info{$gmhmm_gene_id}{'cds_count'} ."\";";
			print $OUT "$seq_id\t$source\tCDS\t$start\t$end\t$score\t$strand\t$phase\t$att_CDS\n";
			$pos = $end - 2;
			print $OUT "$seq_id\t$source\tstart_codon\t$pos\t$end\t$score\t$strand\t0\t$att_start\n";
		}
		elsif( $4 eq "Single" )
		{
			print $OUT "$seq_id\t$source\texon\t$start\t$end\t0\t$strand\t.\t$att_CDS\n";
			$pos = $start + 2;
                        print $OUT "$seq_id\t$source\tstop_codon\t$start\t$pos\t$score\t$strand\t0\t$att_stop\n";
			$att_CDS .= " cds_type \"Single\"\; count \"". $cds_id ."_". $gene_info{$gmhmm_gene_id}{'cds_count'} ."\";";
			if ( $format eq 'gtf2' )
			{
	                        if( $end - $start + 1 > 3)
        	                {
                	                print $OUT "$seq_id\t$source\tCDS\t". ($start +3)."\t$end\t$score\t$strand\t$phase\t$att_CDS\n";
                        	}
			}
			else
			{
				print $OUT "$seq_id\t$source\tCDS\t$start\t$end\t$score\t$strand\t$phase\t$att_CDS\n";
			}
			$pos = $end - 2;
                        print $OUT "$seq_id\t$source\tstart_codon\t$pos\t$end\t$score\t$strand\t0\t$att_start\n";
		}
	}
  }
}

close $IN;
close $OUT;

print "# done\n" if $v;

exit 0;

# ----------------- sub --------------------------
sub GetLastId
{
	my $name = shift;
	my $i = 0;
	
	if( ! -e $name )
	{
		return $i;
	}
	
	my $fp =  File::ReadBackwards->new( $name ) or die "error on open $name $!" ;
	while( !$fp->eof )
	{
		if ( $fp->readline() =~ /gene_id \"(\d+)_g\";/ )
		{
			$i = $1;
			last;
		}
	}

	if ( $v )
	{
		print "# number of genes currently in the output file: $i\n";
	}
	
	return $i;
}
# ------------------------------------------------
sub ParseDefline
{
	my $line = shift;
	
	my $id = "seq";
	my $coord = 1;
	
	if(  $line =~ /^>?\S+\s+(.*)\s+(\d+)\s+\d+\s*$/ )
	{
		$id = $1;
		$id =~ s/^\s+//;
		$id =~ s/\s+$//;
		
		$coord = $2;
	}
	
	if ($v)
	{
		print "# global sequence ID: $id\n";
		print "# to get global coordinates shift coordinates by: $coord\n";
	}

	return ( $id, $coord );
}
# ------------------------------------------------
sub GetDefline
{
	my $name = shift;

	my $defline = '';

	open( my $IN, $name ) or die "error on open file $0: $name\n$!\n";
	while(<$IN>)
	{
		if ( /^FASTA defline:\s+(.*)$/ )
		{
			$defline = $1;
			last;
		}
	}
	close $IN;

	if ( ! $defline )
		{ die "error, defline is missing in the input file: $name\n"; } 

	print "# defline: $defline\n" if $v;

	return $defline;
}
# ------------------------------------------------
sub FilterGenesByGeneLength
{
	my $ref = shift;
	my $min = shift;

	my $gene_counter = 0;
	my $gene_count_excluded_by_length = 0;

	foreach my $key ( sort{$a<=>$b} keys %{$ref})
	{
		if ( $ref->{$key}{'gene_length'} < $min  &&  ! $ref->{$key}{'evidence'} )
		{
			$gene_count_excluded_by_length += 1;
			$ref->{$key}{'new_gene_id'} = 0;
		}
		else
		{
			$gene_counter += 1;
			$ref->{$key}{'new_gene_id'} = $gene_counter;
		}
	}

	if ($v)
	{
		print "# genes excluded by length: $gene_count_excluded_by_length\n";
	}
}
# ------------------------------------------------
sub PopulateGeneInfo
{
	my $name = shift;
	my $ref = shift;
	
	my $counter = 0;
	
	open( my $IN, $name ) or die "error on open file $0: $name\n$!\n";
	while(<$IN>)
	{	
		if( /^\s*(\d+)\s+(\d+)\s+([\+-])\s+(\w+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([-+])\s+([-+])\s*/ )
		{
			my $gene_id = $1;
			my $cds_id = $2;
			my $cds_length = $7;
			my $L_evi = $10;
			my $R_evi = $11;

			$ref->{$gene_id}{'cds_count'} += 1;
			$ref->{$gene_id}{'gene_length'} += $cds_length;

			$ref->{$gene_id}{'evidence'} += 0;

			if ( $L_evi eq "+" or $R_evi eq "+" )
			{
				$ref->{$gene_id}{'evidence'} += 1;
			}

			if ( $debug )
			{
				if ( ! exists $ref->{$gene_id}{'debug_cds_count'} )
				{
					$ref->{$gene_id}{'debug_cds_count'} = $cds_id;
				}
				elsif ( $ref->{$gene_id}{'debug_cds_count'} < $cds_id )
				{
					$ref->{$gene_id}{'debug_cds_count'} = $cds_id;
				}
			}
		}

		if ( /^>gene_(\d+)\|GeneMark.hmm\|(\d+)_nt/ )
		{
			my $gene_id = $1;
			my $gene_length = $2;

			if ( $debug )
			{
				$ref->{$gene_id}{'debug_gene_length'} = $gene_length;
			}
		}
	}
	close $IN;

	if ( $debug )
	{
		foreach my $key (keys %{$ref})
		{
			if ( $ref->{$key}{'debug_cds_count'} != $ref->{$key}{'cds_count'} )
				{ die "error, number of CDS in a gene differs: $ref->{$key}{'debug_cds_count'} $ref->{$key}{'cds_count'}\n"; }

			if ( $ref->{$key}{'debug_gene_length'} !=  $ref->{$key}{'gene_length'} )
				{ die "error, number of CDS in a gene differs: $ref->{$key}{'debug_gene_length'} $ref->{$key}{'gene_length'}\n"; }
		}
	}

	if ($v)
	{
		print "# genes in: ". (scalar (keys %{$ref}) ) ."\n";
	}
}
# ------------------------------------------------
sub CheckInput
{
	if( !$infile or ! -e $infile )  { print "error, input file is misssing $0\n"; exit 1; }
	if( !$outfile  ) { print "error, outfile file name is misssing $0\n"; exit 1; }
	if ( $infile eq $outfile ) { die "error, input file name is idetical with output file name: $infile $outfile\n"; }
}
# ------------------------------------------------
sub ParseCMD
{
	my $cmd = $0;
	foreach my $str (@ARGV) { $cmd .= ( ' '. $str ); }
	
	my %h;
	my $opt_result = GetOptions
	(
		\%h,
		'infile=s'  => \$infile,
		'outfile=s' => \$outfile,
		'app'       => \$app,
		'min=i'     => \$min,
		'verbose'   => \$v,
		'debug'     => \$debug,
		'format=s'  => \$format
	);
	
	if( !$opt_result ) { print "error on command line\n"; exit 1; }
	if( @ARGV > 0 ) { print "error, unexpected argument found on command line: @ARGV\n"; exit 1; }
	$v = 1 if $debug;
}
# -------------------------------------------------------------
sub Usage
{
	my $text =
"
Usage: $0 --in [input file] --out [output file] [optional]

   [input file]  - name of file with predictions by GeneMark.hmm eukaryotic
   [output file] - name of file to save predictions in GTF format

Optional

   --app      - append to output file
   --min      - [$min] filter out short gene predictions
   --format   - [$format] output format; 
                'gtf'  includes stop codon into CDS
                'gtf2' excludes stop codon from CDS
";
	print $text;
	
	exit 1;
}
# -------------------------------------------------------------

