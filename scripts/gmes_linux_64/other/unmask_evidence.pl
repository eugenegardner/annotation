#!/usr/bin/perl
# --------------------------------------------
# Alex Lomsadze
# Georgia Institute of Technology
# 2019
#
# Input is 
#   * softmasked genome in FASTA format
#   * evidence in GFF format
#
# This script unmasks interval borders
# --------------------------------------------

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;

my $Version = "1.02";
# --------------------------------------------
my $in = '';
my $out = '';
my $gff = '';

my $v = 0;
my $debug = 0;
# --------------------------------------------
Usage( $Version ) if ( @ARGV < 1 );
ParseCMD();
CheckBeforeRun();

my %pos = ();
LoadGFF( $gff, \%pos );

my %seq = ();
my %deflines = ();
my @order = ();
LoadFasta( $in, \%seq, \%deflines, \@order );

Unmask( \%seq, \%pos );
PrintOut( $out, \%seq, \%deflines, \@order );

exit 0;

# ================== sub =====================
sub Unmask
{
	my $seqh = shift;
	my $posh = shift;

	foreach my $id ( keys %{$seqh} )
	{
		my $max = length( $seqh->{$id} );

		next if ( !exists $posh->{$id} );

		my $size = scalar ( @{$posh->{$id}} );

		for( my $i = 0; $i < $size; $i += 2 )
		{
			my $L = $posh->{$id}[$i];
			my $R = $posh->{$id}[$i + 1];

			$L = 1 if ( $L < 1 );
			$R = $max if ( $R > $max );

			my $old = substr( $seqh->{$id}, $L - 1, $R - $L + 1 );
			my $new = uc $old;

			if ( $old ne $new )
			{
				if ( $debug )
				{
					print "#\n";
					print "$old\n";
					print "$new\n";
				}
			
				substr( $seqh->{$id}, $L - 1, $R - $L + 1, $new );
			}
		}
	}
}
# --------------------------------------------
sub PrintOut
{
	my $name = shift;
	my $seqh = shift;
	my $defh = shift;
	my $order_a = shift;

	open( my $OUT, ">", $name ) or die "error on open file $name: $!\n";

	foreach my $id (@{$order_a} )
	{
		print $OUT $defh->{$id};

		my $size = length( $seqh->{$id} );
		for( my $i = 0; $i < $size; $i += 50 )
		{
			print $OUT substr( $seqh->{$id}, $i, 50 );
			print $OUT "\n";
		}
	}

	close $OUT;
}
# --------------------------------------------
sub LoadGFF
{
	my $name = shift;
	my $ref = shift; 

	open( my $IN, $name ) or die "error on open file $name: $!\n";
	while( my $line = <$IN> )
	{
		next if( $line =~ /^\s*#/ );
		next if( $line =~ /^\s*$/ );

		if( $line =~ /^(\S+)\t\S+\t(\S+)\t(\d+)\t(\d+)\t\S+\t([-+.])\t(\S+)(\t(.*)|\s*)$/ )
		{
			my $id     = $1;
			my $type   = $2;
			my $start  = $3;
			my $end    = $4;
			my $strand = $5;
			my $ph     = $6;
			my $attr   = $7;

			if ( $type =~ /^[Ii]ntron$/ )
			{
				if ( $strand eq "+" )
				{
					push @{$ref->{$id}}, $start - 3;
					push @{$ref->{$id}}, $start + 5;

					push @{$ref->{$id}}, $end - 19;
					push @{$ref->{$id}}, $end + 1;
				}
				elsif ( $strand eq "-" )
				{
					push @{$ref->{$id}}, $start - 1;
					push @{$ref->{$id}}, $start + 19;

					push @{$ref->{$id}}, $end - 5;
					push @{$ref->{$id}}, $end + 3;
				}

				# PushData( $start, \@{$ref->{$id}} );
				# PushData( $end, \@{$ref->{$id}} );
			}
			elsif ( $type =~ /^[Ss]tart_codon$/ )
			{
				if ( $strand eq "+" )
				{
					push @{$ref->{$id}}, $start - 6;
					push @{$ref->{$id}}, $start + 5;
				}
				elsif ( $strand eq "-" )
				{
					push @{$ref->{$id}}, $end - 5;
					push @{$ref->{$id}}, $end + 6;
				}
				# PushData( $start, \@{$ref->{$id}} );
			}
			elsif ( $type =~ /^[Ss]top_codon$/ )
			{
				if ( $strand eq "+" )
				{
					push @{$ref->{$id}}, $start - 3;
					push @{$ref->{$id}}, $start + 8;
				}
				elsif ( $strand eq "-" )
				{
					push @{$ref->{$id}}, $end - 8;
					push @{$ref->{$id}}, $end + 3;
				}
				# PushData( $start, \@{$ref->{$id}} );
			}
			elsif ( $type =~ /^CDS$/ )
			{
				push @{$ref->{$id}}, $start - 19;
				push @{$ref->{$id}}, $end + 19;
			}
		}
		else { die "error, unexpected GFF line format found: $line\n"; }
	}
	close $IN;

	if ($v)
	{
		my $size = scalar ( keys %{$ref} );
		print "# number of records in GFF file: $size\n";

		foreach my $id ( keys %{$ref} )
		{
			$size = scalar @{$ref->{$id}};
			print "# number of positions in $id: $size\n";
		}
	}

	print Dumper($ref) if $debug;
}
# --------------------------------------------
sub PushData
{
	my $pos = shift;
	my $ref = shift;

	my $L = $pos - 20;
	$L = 1 if ($L < 1);

	push @{$ref}, $L;
}
# --------------------------------------------
sub LoadFasta
{
	my $name = shift;
	my $ref_s = shift;
	my $ref_d = shift;
	my $ref_order = shift;

	my $id = '';

	open( my $IN, $name ) or die "error on open file $0: $name, $!\n";
	while( my $line = <$IN> )
	{
		# skip empty
		if ( $line =~ /^\s*$/ ) { next; }

		# if defline was found
		if ( $line =~ /^>/ )
		{
			if ( $line =~ /^>\s*(\S+)\s*/ )
			{
				$id = $1;
				$ref_d->{$id} = $line;
				$ref_s->{$id} = "";
				push @{$ref_order}, $id;
			}
			else
				{ die "error, unexpected defline format in FASTA file: $line\n"; }
		}
		else
		{
			# parse sequence
			
			# remove white spaces, etc
			$line =~ tr/0123456789\n\r\t //d;
			
			$ref_s->{$id} .= $line;
		}
	}
	close $IN;
	
	if ( !$id ) { die "error, defline is not found in: $name\n"; }

	if ($v)
	{
		my $size = scalar ( keys %{$ref_d} );
		print "# number of records in FASTA file: $size\n";
	}

	if ($debug)
	{
		print Dumper($ref_d);
		print Dumper($ref_s);
		print Dumper($ref_order);
	}
}
# --------------------------------------------
sub CheckBeforeRun
{
	if ( !$out ) { print "error, output file name is missing: $0\n"; exit 1; }

	if ( !$in) { print "error, input file name is missing: $0\n"; exit 1; }
	if ( !-e $in )
		{ print "error, file not found $0: $in\n"; exit 1; }

	if ( !$gff) { print "error, input GFF file name is missing: $0\n"; exit 1; }
	if ( !-e $gff )
		{ print "error, file not found $0: $gff\n"; exit 1; }

	if (($out eq $in) or ($out eq $gff))
		{ die "error, input and output names are identivall: $out $in $gff\n"; }
}
# --------------------------------------------
sub ParseCMD
{
	my $opt_result = GetOptions
	(
	  'in=s'    => \$in,
	  'out=s'   => \$out,
	  'gff=s'   => \$gff,
	  'verbose' => \$v,
	  'debug'   => \$debug
	);

	if ( !$opt_result ) { die "error on command line: $0\n"; }
}
# --------------------------------------------
sub Usage
{
	my $version = shift;

	print qq(
Usage: $0 --out [name] --gff [name] --in [name]

  --in    softmasked sequence in FASTA format
  --gff   evidence hints in GFF format
  --out   softmaked sequence with unmasked regions

Optional
  --verbose
  --debug

Version $version
);
	exit 1;
}
# --------------------------------------------

