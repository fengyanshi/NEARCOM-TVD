# --- exit, if already done
if ( -f "swanmain.f" or -f "swanmain.F" or -f "swanmain.for" )
{
  exit;
}

# --- parsing arguments
$mpi = "FALSE";
$f95 = "FALSE";
$omp = "FALSE";
$dos = "FALSE";
$unx = "FALSE";
$cry = "FALSE";
$sgi = "FALSE";
$cvi = "FALSE";
while ( $ARGV[0]=~/-.*/ )
   {
   if ($ARGV[0]=~/-mpi/) {$mpi="TRUE";shift;}
   if ($ARGV[0]=~/-f95/) {$f95="TRUE";shift;}
   if ($ARGV[0]=~/-omp/) {$omp="TRUE";shift;}
   if ($ARGV[0]=~/-dos/) {$dos="TRUE";shift;}
   if ($ARGV[0]=~/-unix/) {$unx="TRUE";shift;}
   if ($ARGV[0]=~/-cray/) {$cry="TRUE";shift;}
   if ($ARGV[0]=~/-sgi/) {$sgi="TRUE";shift;}
   if ($ARGV[0]=~/-cvis/) {$cvi="TRUE";shift;}
   }

# --- make a list of all files
@files = ();
foreach (@ARGV) {
   @files = (@files , glob );
}

# --- change each file if necessary
foreach $file (@files)
{
  open file or die "can't open $file\n";
  if ($unx=~/TRUE/ && $omp=~/TRUE/)
  {
    ($tempf)=split(/.ftn/, $file);
    $outfile = join(".",$tempf,"F");
  }
  elsif ($unx=~/TRUE/)
  {
    ($tempf)=split(/.ftn/, $file);
    $outfile = join(".",$tempf,"f");
  }
  else
  {
    ($tempf)=split(/.ftn/, $file);
    $outfile = join(".",$tempf,"for");
  }
  open(OUTFILE,">".$outfile);
  while ($line=<file>)
  {
    $newline=$line;
    if ($mpi=~/TRUE/) {$newline=~s/^!MPI//;}
    if ($f95=~/TRUE/) {$newline=~s/^!F95//;}
    if ($omp=~/TRUE/) {$newline=~s/^!OMP//;}
    if ($dos=~/TRUE/) {$newline=~s/^!DOS//;}
    if ($unx=~/TRUE/) {$newline=~s/^!UNIX//;}
    if ($cry=~/TRUE/) {$newline=~s/^!\/Cray//;}
    if ($sgi=~/TRUE/) {$newline=~s/^!\/SGI//;}
    if ($cvi=~/TRUE/) {$newline=~s/^!CVIS//;}    
	print OUTFILE $newline;
  }
  close file;
  close(OUTFILE);
}
