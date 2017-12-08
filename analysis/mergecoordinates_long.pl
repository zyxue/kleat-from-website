#!/bin/perl

=pod

=head1 NAME

mergecoordinates - merge coordinates belonging to the same label from a list

=head1 SYNOPSIS

  cat coordinates.txt | mergecoordinates [-depth] [-chr N] [-unique]

=head1 DESCRIPTION

Given a set of coordinates in the form

  chr start end [id]
  chr start end [id]
  chr start end [id]
  chr start end [id]

this script computes the union of all coordinates for a given chromosome and reports all contiguous regions covered by the coordinates. The ID associated with each coordinate is optional and can be, for example, a clone name. 

Each contiguous region, formed by the overlap of individual coordinates, is represented on a single line.

  chr regionstart regionend num_coord [id,id,...]
  chr regionstart regionend num_coord [id,id,...]

If you've supplied IDs with your coordinates, each line will be annotated with the IDs of the coordinates that overlap with a given region. Regardless, the number of such coordinates will be reported (num_coord).

This script is similar to cumulcoverage and also uses Set::IntSpan.

=head2 Depth

You can ask that mergecoordinates report the depth of coverage by using -depth. In this case, an additional field will be added to each covered region - the depth.

 cat bacend.txt | mergecoordinates -depth
 ...
 22 31138 47761719 47792856 5
 22 14826 47792857 47807682 4
 22 2846 47807683 47810528 3
 22 26683 47810529 47837211 4
 22 19031 47837212 47856242 5
 ...

Note that you'll get many more covers than if you don't use -depth, because depth varies more quickly than the binary yes/no coverage.

=head1 EXAMPLE

Suppose you have a list of BAC end coordinates for some clones on chr 22.

  22     38566064     38595824      CTD-3167C10     D3167C10 AQ183431,AQ149440
  22     38566070     38595859      CTD-3167A12     D3167A12 AQ188101,AQ149429
  22     38705436     38753675      CTD-3047E10     D3047E10 AQ133614,AQ133785
  22     40001912     40107721      CTD-2224M12     M2224M12 AQ033449,AQ146103

Many of these clones overlap, but they don't form a single contiguous partition on the chromosome. You want to find out what regions these clones provide coverage for.

  cat bacend.txt | mergecoordinates
  ...
  22 687936 16454715 17142650
  22 1327013 17251046 18578058
  22 331253 18691105 19022357
  22 663918 19057081 19720998
  22 80788 19746635 19827422
  ...

Using -depth splits the coverage coordinates further by depth of coverage. In this mode, the script takes significantly longer to run.

  cat bacend.txt | mergecoordinates -depth
  22 3755 17767926 17771680 5
  22 5622 17771681 17777302 4
  22 18884 17777303 17796186 5
  22 11489 17796187 17807675 6
  22 7258 17807676 17814933 7


=head1 HISTORY

=item * 10 Apr 2007

Added -unique to report only unique ids.

=item * 20 Jan 2005

Added -depth

=head1 BUGS

=head1 AUTHOR

Martin Krzywinski

$Id: script,v 1.7 2003/07/16 17:51:14 martink Exp $

=head1 CONTACT

  Martin Krzywinski
  Genome Sciences Centre
  Vancouver BC Canada
  www.bcgsc.ca
  martink@bcgsc.ca

=cut

################################################################
#
# Copyright 2002, 2003 Martin Krzywinski
#
# This file is part of the Genome Sciences Centre Perl code base.
#
# This script is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This script is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this script; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
################################################################

################################################################
#
#
################################################################
#                           Martin Krzywinski (martink@bcgsc.ca)
#                                                           2003
################################################################
# $Id: script,v 1.7 2003/07/16 17:51:14 martink Exp $
################################################################

use strict;
use Config::General;
use Data::Dumper;
use File::Basename;
use FindBin;
use Getopt::Long;
use IO::File;
use Math::VecStat qw(sum min max average);
use Pod::Usage;
use Set::IntSpan;
use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/lib";
use vars qw(%OPT %CONF);

use lib "/home/martink/cvs/setcover";
use SetCover;

################################################################
#
# *** YOUR MODULE IMPORTS HERE
#
################################################################

GetOptions(\%OPT,
	   "depth",
	   "chr=s",
	   "unique",
	   "mask=s",
	   ### ADD ANY ADDITIONAL COMMAND LINE FLAGS HERE
	   "configfile=s","delim=s","help","man","debug+","sep=s");

pod2usage() if $OPT{help};
pod2usage(-verbose=>2) if $OPT{man};
loadconfiguration($OPT{configfile});
populateconfiguration(); # copy command line options to config hash
validateconfiguration(); 
if($CONF{debug} > 1) {
  $Data::Dumper::Pad = "debug parameters";
  $Data::Dumper::Indent = 1;
  $Data::Dumper::Quotekeys = 0;
  $Data::Dumper::Terse = 1;
  print Dumper(\%CONF);
}

################################################################
#
# *** YOUR CODE HERE ***
#
################################################################

my %sets;
my $sep = exists($OPT{sep})? $OPT{sep} : ",";

# we're going to read in the coordinates and keep a hash of them, keyed by chr
while(<>) {
  chomp;
#  my @tok = split(/[$CONF{delim}]+/,$_);
  my ($chr,$start,$end,$id) = $_ =~ /^(\S+)\s(\S+)\s(\S+)\s(.*)/;
  next unless $chr;
  $id ||= "-";
  next if $CONF{chr} && $chr ne $CONF{chr};
  ($start,$end) = ($end,$start) if $end < $start;
  my $set = Set::IntSpan->new("$start-$end");
  push(@{$sets{$chr}},{set=>$set,id=>$id});
}

foreach my $chr (sort keys %sets) {
  my @sets = sort { ($a->{set}->min <=> $b->{set}->min) || ($a->{set}->max <=> $b->{set}->max) } @{$sets{$chr}};
  my $currset = undef;
  my @currid;
  my @currsets;
  for my $set (@sets) {
      if(! defined $currset || $set->{set}->intersect($currset)->cardinality) {
	  push(@currid,$set->{id}) if $set->{id};
	  push(@currsets,[$set->{set},$set->{id}]) if $CONF{depth};
	  $currset ||= Set::IntSpan->new();
	  $currset = $currset->union($set->{set});
      } else {
	  if($CONF{depth}) {
	      my @covers = SetCover::GetCovers(\@currsets);
	      foreach my $cover (@covers) {
		  my ($set,$id) = @$cover;
		  my $covertext;
		  my $ncovertext;
		  if($CONF{unique}) {
		    my %ids;
		    map { $ids{$_}++ } @$id;
		    $covertext = join($sep,sort keys %ids);
		    $ncovertext = int(keys %ids);
		  } else {
		    $covertext = join($sep,sort @$id);
		    $ncovertext = @$id;
		  }
		  printinfo($chr,
			    $set->min,$set->max,
			    $set->cardinality,
			    $ncovertext,$covertext);
	      }
	  } else {
	    my $covertext;
	    my $ncovertext;
	    if($CONF{unique}) {
	      my %ids;
	      map { $ids{$_}++ } @currid;
	      $covertext = join($sep,sort keys %ids);
	      $ncovertext = int(keys %ids);
	    } else {
	      $covertext = join($sep,sort @currid);
	      $ncovertext = @currid;
	    }
	    printinfo($chr,$currset->min,$currset->max,$currset->cardinality,$ncovertext,$covertext);
	  }
	  $currset = Set::IntSpan->new($set->{set}->run_list);
	  @currid = ($set->{id});
	  @currsets = ([$set->{set},$set->{id}]);
      }
  }
  if($CONF{depth} && @currsets) {
      my @covers = SetCover::GetCovers(\@currsets);
      foreach my $cover (@covers) {
	  my ($set,$id) = @$cover;
	  $id ||= [];
	  my $covertext;
	  my $ncovertext;
	  if($CONF{unique}) {
	    my %ids;
	    map { $ids{$_}++ } @$id;
	    $covertext = join($sep,sort keys %ids);
	    $ncovertext = int(keys %ids);
	  } else {
	    $covertext = join($sep,sort @$id);
	    $ncovertext = @$id;
	  }
	  printinfo($chr,
		    $set->min,$set->max,
		    $set->cardinality,
		    $ncovertext,$covertext);
      }
  } elsif ($currset) {
    my $covertext;
    my $ncovertext;
    if($CONF{unique}) {
      my %ids;
      map { $ids{$_}++ } @currid;
      $covertext = join($sep,sort keys %ids);
      $ncovertext = int(keys %ids);
    } else {
      $covertext = join($sep,sort @currid);
      $ncovertext = @currid;
    }
    printinfo($chr,$currset->min,$currset->max,$currset->cardinality,$ncovertext,$covertext);
  }
}

exit 0;

################################################################
#
# *** YOUR FUNCTIONS HERE ***
#
################################################################


################################################################
#
# *** VALIDATE YOUR CONFIGURATION PARAMETERS HERE
#
################################################################

sub validateconfiguration {
  # $CONF{test} ||= 5;
  $CONF{delim} ||= "\\s";
}

################################################################
#
# *** DO NOT EDIT BELOW THIS LINE ***
#
################################################################
################################################################
################################################################
################################################################


sub populateconfiguration {
  foreach my $key (keys %OPT) {
    $CONF{$key} = $OPT{$key};
  }
}

sub loadconfiguration {
  my $file = shift;
  my ($scriptname) = fileparse($0);
  if(-e $file && -r _) {
    # great the file exists
  } elsif (-e "/home/$ENV{LOGNAME}/.$scriptname.conf" && -r _) {
    $file = "/home/$ENV{LOGNAME}/.$scriptname.conf";
  } elsif (-e "$FindBin::RealBin/$scriptname.conf" && -r _) {
    $file = "$FindBin::RealBin/$scriptname.conf";
  } elsif (-e "$FindBin::RealBin/etc/$scriptname.conf" && -r _) {
    $file = "$FindBin::RealBin/etc/$scriptname.conf";
  } else {
    return undef;
  }
  $OPT{configfile} = $file;
  my $conf = new Config::General(-ConfigFile=>$file,
				 -AllowMultiOptions=>"yes",
				 -LowerCaseNames=>1,
				 -AutoTrue=>1);
  %CONF = $conf->getall;
}

sub printdebug {
  printinfo("debug",@_)  if $CONF{debug};
}

sub printinfo {
  my $message = shift;
  $message .= " " if $message;
  for (@_) {
    if (/^-?\d+\.[\de-]{5,}$/) {
      if(/e/) {
	$_ = sprintf("%0.4e",$_);
      } else {
	$_ = sprintf("%0.4f",$_);
      }
    }
  }
  printf("%s%s\n",$message,join(" ",@_));
}

