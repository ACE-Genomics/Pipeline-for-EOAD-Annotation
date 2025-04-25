#!/usr/bin/perl
use strict;
use warnings;
use SLURMACE qw(send2slurm);
############################################
###### Variables de ejecucion ##############
############################################
# Tiempo maximo de ejecucion de cada proceso
my $time = '72:0:0';
# Numero de CPUs de cada proceso
my $cpus_per_proc =4;
# Memoria a usar por cada CPU 
# Si no estas seguro de lo que haces no lo toques
my $mem_per_cpu = '4G';
# Particion del cluster a usar
# Si no estas seguro,dejarlo en el defaulta
my $partition = 'fast';
# Directorio para almacenar los scripts y logs
my $wdir = 'slurm';
############################################
##### No editar a partir de aqui ###########
############################################
my $ifile = $ARGV[0];
die "Should supply input file\n" unless $ifile;
mkdir $wdir;
my $count = 0;
open IPDF, "<$ifile" or die "Could not open input file\n$!\n";

my %ptask = ( 'job_name' => $ifile,
	'cpus' => $cpus_per_proc,
	'mem_per_cpu' => $mem_per_cpu,
	'time' => $time, 
	'mailtype' => 'FAIL,TIME_LIMIT,STAGE_OUT',
	'partition' => $partition,	
);

while (<IPDF>) {
	$count++;
	my $ofile = sprintf ("%s_%04d", 'sorder', $count);
	$ptask{'filename'} = $wdir.'/'.$ofile.'.sh';
	$ptask{'output'} = $wdir.'/'.$ifile.'-%j'; 
	$ptask{'command'} = $_;
	send2slurm(\%ptask);
}
close IPDF;
my %warn = ('job_name' => $ifile,         
	'filename' => $wdir.'/tasks_end.sh',         
	'mailtype' => 'END',         
	'output' => $wdir.'/tasks_end',
	'dependency' => 'singleton' 
); 
send2slurm(\%warn); 

