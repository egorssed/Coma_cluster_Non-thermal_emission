#! /usr/bin/perl
#require "utils.pl";
#require "getopts.pl";

#&Getopts('i:o:n:c:e');

#############################################
$alpha=2.1;
$norm=10;


$mean_rms=0.73;


$inp=$ARGV[0];
$out=$inp;
open INPP,">$inp.dat";
@spec=`cat $inp`;
chomp @spec;
$cnt1=0;
foreach $sp (@spec){
    ($null,$e1[$cnt1],$e2[$cnt1], $rate[$cnt1], $drate[$cnt1], $null)=split / +/,$sp;
    $rate[$cnt1]=$rate[$cnt1]/1000;
    $drate[$cnt1]=$drate[$cnt1]/1000;
    $crrate[$cnt1]=1;
    $dcrrate[$cnt1]=0.001;
    print INPP "$cnt1 $rate[$cnt1] $drate[$cnt1]\n";
    $cnt1++;
}
close(INPP);



$TFORM2=sprintf "%dE",$cnt1;
$rmff="$out.rmf";
open COLS,">cols1";
print COLS "ENERG_LO 1e keV
ENERG_HI 1e keV
N_GRP I
F_CHAN 1I
N_CHAN 1I
MATRIX $TFORM2
";
close COLS;
open COLS,">cols2";
print COLS "CHANNEL J
E_MIN E keV
E_MAX E keV
";
close COLS;

open FITS, "|fcreate cols1  datafile=\"-\" outfile=\"$rmff\" tbltype=\"binary\" extname=\"SPECRESP MATRIX\" clobber=yes";


$beta=1-$alpha;
for($i=0;$i<$cnt1;){
    $ener=($e1[$i]+$e2[$i])/2;
    $crfl[$i]=$norm*(exp($beta*log($e2[$i]))-exp($beta*log($e1[$i])))/$beta;
    $rmf[$i]=$crrate[$i]/$crfl[$i];
    print"$i $e1[$i] $e2[$i] $rmf[$i] $crfl[$i]\n";
    $num=$i+1;
    print FITS "$e1[$i] $e2[$i] 1 0 $cnt1";
    for($j=0;$j<$cnt1;){
	$koef=0;
	if($i==$j){$koef=$rmf[$i];}
	print FITS " $koef";
	$j++;
    }

    print FITS "\n";

    $i++;

}

close FITS;

open FITS, "|fcreate cols2  datafile=\"-\" outfile=\"temp.fits\" tbltype=\"binary\" extname=\"EBOUNDS\" clobber=yes";

for($i=0;$i<$cnt1;){
    $num=$i+1;
    print FITS "$i $e1[$i] $e2[$i]\n";

    $i++;

}

close FITS;

system "fappend temp.fits[1] $rmff";
system "rm -f temp.fits";

system "fparkey OGIP $rmff\[1] HDUCLASS add=yes";
system "fparkey RESPONSE $rmff\[1] HDUCLAS1 add=yes";
system "fparkey RSP_MATRIX $rmff\[1] HDUCLAS2 add=yes";
system "fparkey $cnt1 $rmff\[1] DETCHANS add=yes";
system "fparkey PHA $rmff\[1] CHANTYPE add=yes";
system "fparkey 1.0 $rmff\[1] EFFAREA add=yes";
system "fparkey FULL_RES $rmff\[1] BINNING add=yes";
system "fparkey 1992a $rmff\[1] RMFVERSN add=yes";
system "fparkey INTEGRAL $rmff\[1] TELESCOP add=yes";
system "fparkey ISGRI  $rmff\[1] INSTRUME add=yes";
system "fparkey DAL_TABLE $rmff\[1] BASETYPE add=yes";
system "fparkey 1 $rmff\[1] EXTVER add=yes";
system "fparkey 1 $rmff\[1] GRPID1 add=yes";
system "fparkey 8 $rmff\[1] TLMAX4 add=yes";
system "fparkey 0 $rmff\[1] TLMIN4 add=yes";

system "fparkey OGIP $rmff\[2] HDUCLASS add=yes";
system "fparkey RESPONSE $rmff\[2] HDUCLAS1 add=yes";
system "fparkey EBOUNDS $rmff\[2] HDUCLAS2 add=yes";
system "fparkey $cnt1 $rmff\[2] DETCHANS add=yes";
system "fparkey PHA $rmff\[2] CHANTYPE add=yes";
system "fparkey 1.0 $rmff\[2] EFFAREA add=yes";
system "fparkey DAL_TABLE $rmff\[2] BASETYPE add=yes";
system "fparkey 1992a $rmff\[2] RMFVERSN add=yes";
system "fparkey 1.3.0 $rmff\[2] HDUVERS add=yes";
system "fparkey 1.0.0 $rmff\[2] HDUVERS1 add=yes";
system "fparkey 1.3.0 $rmff\[2] HDUVERS2 add=yes";
system "fparkey 8 $rmff\[2] TLMAX1 add=yes";
system "fparkey 0 $rmff\[2] TLMIN1 add=yes";

system "fparkey 8 $rmff\[0] BITPIX add=yes";

#system "cphead crab_1.15.rmf $rmff";
#system "cphead crab_1.15.rmf\[1] $rmff\[1]";
#system "cphead crab_1.15.rmf\[2] $rmff\[2]";
#system "fparkey $cnt1 $rmff\[2] DETCHANS add=yes";


###### pha file #########


open(COLS,">cols");
print COLS "CHANNEL                    1I\n";
print COLS "RATE                       1D                  counts\n";
print COLS "STAT_ERR                   1D                  counts\n";

close(COLS);

open(HEADER,">header");
print HEADER <<EOHEAD;
EXTNAME = 'SPECTRUM'           / name of this binary table extension
HDUCLASS= 'OGIP    '           / format conforms to OGIP/GSFC standards
HDUCLAS1= 'SPECTRUM'           / Extension contains a Spectrum
HDUCLAS2= 'TOTAL   '           / Extension contains a Spectrum
HDUCLAS3= 'COUNTS   '           / Extension contains counts
HDUVERS1= '1.1.0   '           / Version number of the format
POISSERR=                    F / Are Poisson Distribution errors assumed.
SYS_ERR =                    0 / No systematic error was specified
GROUPING=                    0 / No grouping data has been specified
QUALITY =                    0 / No data quality information specified
TELESCOP= 'INTEGRAL'           / Telescope (mission) name
INSTRUME= 'ISGRI   '           / Instrument name
FILTER  = 'NONE    '           / Instrument filter in use
EXPOSURE= 1.00000000000000E+00 / Exposure time
AREASCAL=       1.00000000E+00 / Nominal effective area
BACKSCAL=       1.00000000E+00 / Background scale factor
CORRSCAL=       0.00000000E+00 / Correlation scale factor
BACKFILE= 'NONE    '           / Background FITS file for this object
CORRFILE= 'NONE    '           / Correlation FITS file for this object
RESPFILE= '$rmff'         / Redistribution matrix file (RMF)
ANCRFILE= 'NONE    '           / Ancillary response file (ARF)
XFLT0001= 'NONE    '           / XSPEC selection filter description
CHANTYPE= 'PHA     '           / Channels assigned by detector electronics
LONGSTRN= 'OGIP 1.0'           / The HEASARC Long String Convention may be used.
DETCHANS= $cnt1              / Number of channels in file
OBJECT  = '${object}'           / OBJECT from the FIRST input file
ROWID1  = 'XeCnt   '           / Column Name processed
ORIGIN  = 'NASA/GSFC'          / origin of fits file
CREATOR = 'do_pha_igr version 2.0' / Program name that produced this file
DATE    = '2008-12-05T13:58:43' / file creation date (YYYY-MM-DDThh:mm:ss UTC)
RA_OBJ  =       0.0            / RA of First input object
DEC_OBJ =       0.0            / DEC of First input object
EQUINOX =              2000.00 / Equinox of the FIRST object
RADECSYS= 'FK5     '           / Co-ordinate frame used for equinox
FREQUEN=       $fcentr          / Central frequency of the range, Hz
ER_FREQ=       $ef             / the width of the frequency range, Hz
PHAVERSN= '1992a   '           / OGIP memo number for file format
TIMESYS = 'TT      '           / The time system is MJD
GAINAPP =                    T / Gain all ready subracted
COMMENT   This file contents the spectrum of the source in the definite frequency range
EOHEAD
close(HEADER);

system "fcreate headfile=header cdfile=cols datafile=$inp.dat outfile=$out.pha clobber=yes";
