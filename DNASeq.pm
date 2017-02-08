#!/usr/bin/perl
use lib '/crihome/zlei2/gatk_pipeline';
use FileSystem;
use Algorithm;

################################################################
# DNASeq Package											   #
# Author:   Zhengdeng Lei			                           #
# E-mail:   leizhengdeng@hotmail.com                           #
################################################################ 

# Example in caller DNASeq.pl
# Assume the paths for fastx_tk, maq, bwa, samtools are in ENV{PATH}; check $PATH setting $HOME/.bash_profile
# Assume perl and python are installed
#my $NGS_obj = new DNASeq ('sample_name'	=>	'Testis_T28',
#								'input_fq1'		=>	'Run1_testicular-28T_lane2_read1_sequence.txt', 
#								'input_fq2'		=>	'Run1_testicular-28T_lane2_read2_sequence.txt', 
#								'pe'			=>	1,
#								'reference_fa'	=>	'/data/public/reference_genomes/genome_hg18.fa',
#								'dbsnp_rod'		=>	'/data/nextgen1/pipeline/dbsnp_130_hg18.rod',
#								'picard_dir'	=>	'/home/gmslz/tools/picard-tools-1.44',
#								'gatk_dir'		=>	'/data/public/tools/gatk-v4333/Sting',
#								'tmp_dir'		=>	'/home/tmp',
#								'java_mem'		=>	'-Xmx4g'); 
#
#$NGS_obj -> check_fastq_format();
#$NGS_obj -> illumina2sanger();
#$NGS_obj -> qcstats();
#$NGS_obj -> fastq_trimmer(1, 74);
#$NGS_obj -> aln_bwa();					#alignment step, which includes bwa, sort, rm_pcrdup, add_rg
#$NGS_obj -> snp_analysis_gatk();		#gatk step, which includes realignment, recalibration, genotyping (indel, snp)

package DNASeq;
#use strict;
use warnings;
use Time::Local;
use Time::localtime;
use Time::Local 'timelocal_nocheck';
use Cwd();

# private members: 
# sample_name <sample name>
# input_fq1 <input fastq1>
# input_fq2 <input fastq2>
# sanger_fq1 <input fastq1 in sanger format>
# sanger_fq2 <input fastq1 in sanger format>
# reference_fa <reference genome> 
# fastq_format <fastq format>
	# "Sanger"						-- quality score from 0 to 93 using ASCII 33 to 126
	# "Solexa"		Illumina 1.0	-- quality score from -5 to 62 using ASCII 59 to 126
	# "Illumina"	Illumina 1.3+	-- quality score from 0 to 62 using ASCII 64 to 126
# pe <paired end>
	# 1	paired end
	# 0 single end
#my @attributes = qw/reference_fa sample_name input_fq1 input_fq2 fastq_format sanger_fq1 sanger_fq2 pe/;

######################################################################################################
# The subroutines "new" and "initialize" for the constructor
sub new   
{
	my $class_name = shift;		# Gets class name from parmlist
	my $this = {};				# Creates reference to object
	bless($this, $class_name);	# Associates reference with class
	$this->initialize(@_);		# Call local initialize subroutine
								# passing rest of parmlist
	return $this;				# Explicit return of value of $this
}

sub initialize
{
	my $this = shift;    # Receive an object reference as 1st param
	my %parms = @_;         # Get passed parameters from the call to "new"

	$this->{parameters_file}	= $parms{parameters_file}	|| 'parameters.csv';
    my @parameters = FileSystem::ReadFile($this->{parameters_file});
	chomp(@parameters);
	#print "start reading variables\n";

	for (my $i=1; $i<=$#parameters; $i++)
	{
		#print "$i: \n";
		my ($module, $attribute, $value, $default_value, $comment) =  split(/,/, $parameters[$i]);
		if (defined($attribute))
		{
			#print "attr = $attribute\n";
			$attribute = Algorithm::trim_space($attribute); #remove space from both ends of the string
			if (defined($value))
			{
				$value = Algorithm::trim_space($value);
				$value =~ s/""/'/g;
				$value =~ s/"//g;
				#$value =~ s/'//g;

			}
			if (defined($default_value)) 
			{
				$default_value = Algorithm::trim_space($default_value);
				$default_value =~ s/""/'/g;
				$default_value =~ s/"//g;
				#$value =~ s/'//g;

			}
			if ($attribute ne "") 
			{
				#$this->{$attribute} = $value || $default_value;
				if ($value ne "")       
				{       
						$this->{$attribute} = $value; 
				} 
				else 
				{ 
						 $this->{$attribute} = $default_value;
				}
			}
		}

	}

	if ($this->{java_mem} =~ /[g|G]$/)
	{
		$this->{java_mem} = '-Xmx'.$this->{java_mem};
	}
	else
	{
		$this->{java_mem} = '-Xmx'.$this->{java_mem}.'g';
	}

	if ($this->{dry_run})
	{
		$this->{cluster_mode} = 0;
	}
	#print "variables read\n";

	$this->{input_data_dir}	= FileSystem::NoSlashPath($this->{input_data_dir});
	$this->{output_project_dir}	= FileSystem::NoSlashPath($this->{output_project_dir});
	$this->{gatk_dir}	= FileSystem::NoSlashPath($this->{gatk_dir});
	$this->{picard_dir}	= FileSystem::NoSlashPath($this->{picard_dir});
	$this->{tmp_dir}	= FileSystem::NoSlashPath($this->{tmp_dir});

	$this->{picard_MarkDuplicates}			= "java $this->{java_mem} -jar $this->{picard_dir}/MarkDuplicates.jar";
	$this->{picard_AddOrReplaceReadGroups}	= "java $this->{java_mem} -jar $this->{picard_dir}/AddOrReplaceReadGroups.jar";
	$this->{picard_SortSam}					= "java $this->{java_mem} -jar $this->{picard_dir}/SortSam.jar";
	$this->{picard_BuildBamIndex}			= "java $this->{java_mem} -jar $this->{picard_dir}/BuildBamIndex.jar";
	$this->{picard_FixMateInformation}		= "java $this->{java_mem} -Djava.io.tmpdir=$this->{tmp_dir}/ -jar $this->{picard_dir}/FixMateInformation.jar";
	$this->{trimmomatic}					= "java $this->{java_mem} -jar $this->{trimmomatic_path}";
	$this->{gatk_GenomeAnalysisTK}			= "/cri/home2/zlei2/tools/java/jre1.7.0_65/bin/java $this->{java_mem} -Djava.io.tmpdir=$this->{tmp_dir}/ -jar $this->{gatk_dir}/GenomeAnalysisTK.jar";
	$this->{gatk_makeIndelMask_py}			= "$this->{gatk_dir}/python/makeIndelMask.py";


	#$bam_file = $this->{output_project_dir}."/".$this->{bam_sub_dir}."/".$this->{samples}[$i].".bam";
	$this->{final_recalibrated_bam_file} = "./bam/$this->{sample_name}.final.recalibrated.bam"; #This is the final recalibrated bam file after after gatk

	#print $ENV{'PATH'}, "\n\n\n";
	$this->{log_dir}			= $this->{output_project_dir}."/$this->{log_sub_dir}";
	$this->{bam_dir}			= $this->{output_project_dir}."/$this->{bam_sub_dir}";
	$this->{realigned_bam_dir}	= $this->{output_project_dir}."/$this->{realigned_bam_sub_dir}";

	my @folders = ($this->{output_project_dir});
	push(@folders, $this->{output_project_dir}."/$this->{log_sub_dir}");
	push(@folders, $this->{output_project_dir}."/$this->{before_qc_report_sub_dir}");
	push(@folders, $this->{output_project_dir}."/$this->{after_qc_report_sub_dir}");
	push(@folders, $this->{output_project_dir}."/$this->{qc_fastq_data_sub_dir}");
	push(@folders, $this->{output_project_dir}."/$this->{bam_sub_dir}");
	push(@folders, $this->{output_project_dir}."/$this->{realigned_bam_sub_dir}");
	push(@folders, $this->{output_project_dir}."/$this->{snp_sub_dir}");

	#my @subfolders = qw/qcstat bam snp_analysis_gatk log/;
	foreach my $dir (@folders)
	{
		if (!(-e $dir))
		{
			mkdir($dir);
		}
	}

	$this->{samples} = ();
	my @samples_names = split(/;/, $this->{sample_name});
	$this->{num_of_samples} = $#samples_names + 1;
	for (my $i=0; $i<=$#samples_names; $i++)
	{
		$this->{samples}[$i] = Algorithm::trim_space($samples_names[$i]);
	}

	$this->{fq1_files} = ();
	$this->{fq2_files} = ();
	$this->{qc_fq1_files} = ();
	$this->{qc_fq2_files} = ();

	my @fq1_files = split(/;/, $this->{input_fq1});
	for (my $i=0; $i<=$#fq1_files; $i++)
	{
		$fq1_files[$i] = Algorithm::trim_space($fq1_files[$i]);
		$this->{fq1_files}[$i] = $this->{input_data_dir}."/".$fq1_files[$i];
		$this->{qc_fq1_files}[$i] = $this->{output_project_dir}."/$this->{qc_fastq_data_sub_dir}"."/".$fq1_files[$i];
		$this->{qc_fq1_files}[$i] =~ s/fastq$/qc\.fastq/;

	}
	if ($this->{pe})
	{
		my @fq2_files = split(/;/, $this->{input_fq2});
		for (my $i=0; $i<=$#fq2_files; $i++)
		{
			$fq2_files[$i] = Algorithm::trim_space($fq2_files[$i]);
			$this->{fq2_files}[$i] = $this->{input_data_dir}."/".$fq2_files[$i];
			$this->{qc_fq2_files}[$i] = $this->{output_project_dir}."/$this->{qc_fastq_data_sub_dir}"."/".$fq2_files[$i];
			$this->{qc_fq2_files}[$i] =~ s/fastq$/qc\.fastq/;
		}
	}

	my $current_time = localtime();
	my ($year, $mon, $day, $hour, $min, $sec) = ($current_time->year + 1900, $current_time->mon+1, $current_time->mday, $current_time->hour, $current_time->min, $current_time->sec);
	$this->{pipeline_start_time} = sprintf("%04d%02d%02d_%02d%02d%02d", $year, $mon, $day, $hour, $min, $sec);
	if (!$this->{pipeline_time_stamp})
	{
		$this->{pipeline_time_stamp} = $this->{pipeline_start_time};
	}
	$this->{log_file}	= FileSystem::NoSlashPath($this->{output_project_dir}."/$this->{log_sub_dir}/DNASeq_$this->{pipeline_start_time}.log");
	$this->log_it("####################################PIPELINE STARTS####################################");
	my $log_parameter_csv = $this->{output_project_dir}."/$this->{log_sub_dir}/".FileSystem::GetFileNameFromPath($this->{parameters_file});
	$log_parameter_csv =~ s/\.csv/_$this->{pipeline_start_time}\.csv/;
	
	$this->run_cmd("cp $this->{parameters_file} $log_parameter_csv");
	$this->{report_html} = $this->{output_project_dir}."/$this->{log_sub_dir}/report_$this->{pipeline_time_stamp}.html";

#	print "initialization done\n";

}

#This finalize function is obselete
sub finalize
{
	my $this = shift;  # Figure out who I am
	my $tail = '</div><div class="footer">Produced by Center for Research Informatics, UIC (zlei2@uic.edu)</div></body></html>';
	FileSystem::AppendFile($this->{report_html}, $tail);
	$this->log_it("####################################PIPELINE ENDS####################################");
}

# e.g. $DNASeq_obj -> set ('java_mem'	=>	'8g'); 
sub set  #  modify the object
{
	my $this = shift;  # Figure out who I am
	my %parms = @_;       # Make hash out of rest of parm list
	for(keys %parms)  
	{
		$this->{$_} = $parms{$_}; # Replace with new value 
	}
}

# $DNASeq_obj -> get ('java_mem');
sub get  #  modify the object
{
	my $this = shift;  # Figure out who I am
	my $parm = shift;       # get parameter
	return	$this->{$parm}; 
}



sub is_fastq_sanger
{
	my $this = shift;
	my $sanger_format = 1;
	my $counter = 2000; # count 2000 lines of quality scores
	my $min_qualscore = 1e6;
	my $max_qualscore = -1e6;
	my $tmp;
	my $log_msg;

	if ($this->{dry_run})
	{
		return;
	}


	$this->log_it("####################################SANGER FORMAT CHECKING STARTS####################################");

	my @fq_files = ();
	for (my $i=0; $i < $this->{num_of_samples}; $i++)
	{
		push(@fq_files, $this->{fq1_files}[$i]);
	}
	if($this->{pe}) 
	{
		for (my $i=0; $i < $this->{num_of_samples}; $i++)
		{
			push(@fq_files, $this->{fq2_files}[$i]);
		}
	}

	for (my $i=0; $i<=$#fq_files; $i++)
	{
		$min_qualscore = 1e6;
		$max_qualscore = -1e6;

		unless (-e $fq_files[$i])
		{
			$this->log_it("ERROR (File NOT exist): $fq_files[$i]");
			exit;
		} 
		open(FASTQ, $fq_files[$i]) or die "can't open $fq_files[$i]: $@\n";
		while ($counter > 0) 
		{
			$tmp = <FASTQ>; # @read name
			$tmp = <FASTQ>; # base calls
			$tmp = <FASTQ>; # +read name
			my $scores = <FASTQ>; # quality scores
			if (!$scores) { last; }
			#print $scores;
			chomp($scores);
			for my $chr (map(ord, split(//, $scores))) 
			{
				if ($chr < $min_qualscore) 
				{
					$min_qualscore = $chr;
				}
				if ($chr > $max_qualscore) 
				{
					$max_qualscore = $chr;
				}
			}
			$counter--;
		}
		# Phred+33 means quality score + 33 = ASCII
		if ($min_qualscore >= 33 && $max_qualscore <= 126)
		{
			$log_msg = $fq_files[$i]." is in Sanger format, PASSED!";
			#print $log_msg;
			#FileSystem::WriteFile($this->{log_file}, $log_msg);
			$this->log_it($log_msg);
		}
		else
		{
			$sanger_format = 0;
			$log_msg = $fq_files[$i]." is NOT in Sanger format, FAILED!";
			#print $log_msg;
			#FileSystem::WriteFile($this->{log_file}, $log_msg);
			$this->log_it($log_msg);
		}
	}
	if (!$sanger_format)
	{
		print "One of fastq files is not Sanger format, program terminated\n";
		print "Please check log file", $this->{log_file}, "\n";
		exit;
	}
	$this->log_it("####################################SANGER FORMAT CHECKING ENDS####################################");

}


########################################################################################################
# 2. fastqc
sub fastqc_before_QC
{
	my $this = shift;
	$this->log_it("####################################FASTQC BEFORE QC STARTS####################################");

	my @cmds = ();


	if($this->{pe})
	{
		@cmds = ();
		for (my $i=0; $i < $this->{num_of_samples}; $i++)
		{
			push(@cmds, "fastqc $this->{fq1_files}[$i] -t $this->{cpu_threads} -o ".$this->{output_project_dir}."/$this->{before_qc_report_sub_dir}");

			push(@cmds, "fastqc $this->{fq2_files}[$i] -t $this->{cpu_threads} -o ".$this->{output_project_dir}."/$this->{before_qc_report_sub_dir}");
			$this->run_multiple_cmds(@cmds);

			#$this->run_cmd("fastqc $this->{fq2_files}[$i] -o ".$this->{output_project_dir}."/$this->{before_qc_report_sub_dir}"); #option -Q 33 is for Sanger format
		}
	}
	else
	{
		for (my $i=0; $i < $this->{num_of_samples}; $i++)
		{
			$this->run_cmd("fastqc $this->{fq1_files}[$i] -o ".$this->{output_project_dir}."/$this->{before_qc_report_sub_dir}"); #option -Q 33 is for Sanger format
		}
	}

	#$this->log_it("####################################Testing paralellel####################################");
	
	$this->log_it("####################################FASTQC BEFORE QC ENDS####################################");


	if ($this->{run_pipeline} == 1)
	{
		print "QC only\nExit now\n";
		exit 0;
	}
}

sub fastqc_after_QC
{
	my $this = shift;
	$this->log_it("####################################FASTQC AFTER QC STARTS####################################");
	my @cmds = ();

#	for (my $i=0; $i < $this->{num_of_samples}; $i++)
#	{
#		$this->run_cmd(@cmds, "fastqc $this->{qc_fq1_files}[$i] -t $this->{cpu_threads} -o ".$this->{output_project_dir}."/$this->{after_qc_report_sub_dir}");
#		#$this->run_cmd("fastqc $this->{qc_fq1_files}[$i] -o ".$this->{output_project_dir}."/$this->{after_qc_report_sub_dir}"); #option -Q 33 is for Sanger format
#	}
#	if($this->{pe})
#	{
#		for (my $i=0; $i < $this->{num_of_samples}; $i++)
#		{
#
#			push(@cmds, "fastqc $this->{qc_fq1_files}[$i] -t $this->{cpu_threads} -o ".$this->{output_project_dir}."/$this->{after_qc_report_sub_dir}");
#			push(@cmds, "fastqc $this->{qc_fq2_files}[$i] -t $this->{cpu_threads} -o ".$this->{output_project_dir}."/$this->{after_qc_report_sub_dir}");
#			$this->run_multiple_cmds(@cmds);
#
#			#$this->run_cmd("fastqc $this->{qc_fq2_files}[$i] -o ".$this->{output_project_dir}."/$this->{after_qc_report_sub_dir}"); #option -Q 33 is for Sanger format
#		}
#	}



	if($this->{pe})
	{
		@cmds = ();
		for (my $i=0; $i < $this->{num_of_samples}; $i++)
		{
			push(@cmds, "fastqc $this->{qc_fq1_files}[$i] -t $this->{cpu_threads} -o ".$this->{output_project_dir}."/$this->{after_qc_report_sub_dir}");
			push(@cmds, "fastqc $this->{qc_fq2_files}[$i] -t $this->{cpu_threads} -o ".$this->{output_project_dir}."/$this->{after_qc_report_sub_dir}");
			$this->run_multiple_cmds(@cmds);

		}
	}
	else
	{
		for (my $i=0; $i < $this->{num_of_samples}; $i++)
		{
			$this->run_cmd("fastqc $this->{qc_fq1_files}[$i] -t $this->{cpu_threads} -o ".$this->{output_project_dir}."/$this->{after_qc_report_sub_dir}"); #option -Q 33 is for Sanger format
		}
	}



	$this->log_it("####################################FASTQC AFTER QC ENDS####################################");
}





#	java -jar trimmomatic-0.30.jar PE --phred33 input_forward.fq.gz input_reverse.fq.gz 
#output_forward_paired.fq.gz output_forward_unpaired.fq.gz 
#output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz 
#ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

sub trimmomatic
{

	my $this = shift;
	$this->log_it("####################################TRIMMOMATIC STARTS####################################");
	my @cmds = ();

	for (my $i=0; $i < $this->{num_of_samples}; $i++)
	{
		if($this->{pe})
		{

			$this->run_cmd(@cmds, "$this->{trimmomatic} PE -threads $this->{cpu_threads} -phred33 ".
				"$this->{fq1_files}[$i] $this->{fq2_files}[$i] ".
				"$this->{qc_fq1_files}[$i] $this->{qc_fq1_files}[$i].unpaired ".
				"$this->{qc_fq2_files}[$i] $this->{qc_fq2_files}[$i].unpaired ".
				"LEADING:$this->{trimmomatic_leading} TRAILING:$this->{trimmomatic_trailing} ".
				"SLIDINGWINDOW:$this->{trimmomatic_slidingwindow} MINLEN:$this->{trimmomatic_minlen}");  
		}
		else
		{
			$this->run_cmd(@cmds, "$this->{trimmomatic} SE -threads $this->{cpu_threads} -phred33 ".
				"$this->{fq1_files}[$i] $this->{qc_fq1_files}[$i] ".
				"LEADING:$this->{trimmomatic_leading} TRAILING:$this->{trimmomatic_trailing} ".
				"SLIDINGWINDOW:$this->{trimmomatic_slidingwindow} MINLEN:$this->{trimmomatic_minlen}");  
		}
	}

	#$this->run_multiple_cmds(@cmds);
	$this->log_it("####################################TRIMMOMATIC ENDS####################################");
}

#fastq_quality_filter -Q 33 -q 17 -p 80 -i $f -o $fname_out

sub generate_qc_fastq_data
{
	my $this = shift;
	$this->log_it("####################################GENERATE QC DATA STARTS####################################");

	if ($this->{trim_or_filter} == 1)
	{
		$this->qc_trim();
	}elsif ($this->{trim_or_filter} == 2)
	{
		$this->qc_filter();
	}elsif ($this->{trim_or_filter} == 0)
	{
		$this->trimmomatic();
	}
	$this->log_it("####################################GENERATE QC DATA ENDS####################################");


}

sub qc_no_action
{
}


sub qc_trim
{
	my $this = shift;
	$this->log_it("####################################QC TRIM STARTS####################################");
	my @cmds = ();

	if($this->{pe})
	{
		for (my $i=0; $i < $this->{num_of_samples}; $i++)
		{
			push(@cmds, "fastq_quality_filter -Q 33 -q $this->{filter_quality_score} -p $this->{filter_percent} -i $this->{fq1_files}[$i] -o $this->{qc_fq1_files}[$i].tmp.fastq");
			push(@cmds, "fastq_quality_filter -Q 33 -q $this->{filter_quality_score} -p $this->{filter_percent} -i $this->{fq2_files}[$i] -o $this->{qc_fq2_files}[$i].tmp.fastq");
			$this->run_multiple_cmds(@cmds);
			#$this->run_cmd("fastq_quality_filter -Q 33 -q $this->{filter_quality_score} -p $this->{filter_percent} -i $this->{fq1_files}[$i] -o $this->{qc_fq1_files}[$i].tmp.fastq"); #option -Q 33 is for Sanger format
			#$this->run_cmd("fastq_quality_filter -Q 33 -q $this->{filter_quality_score} -p $this->{filter_percent} -i $this->{fq2_files}[$i] -o $this->{qc_fq2_files}[$i].tmp.fastq"); #option -Q 33 is for Sanger format
			$this->run_cmd("$this->{fix_paired_reads} $this->{qc_fq1_files}[$i].tmp.fastq $this->{qc_fq2_files}[$i].tmp.fastq ".
				"$this->{qc_fq1_files}[$i] $this->{qc_fq2_files}[$i] $this->{qc_fq1_files}[$i].unpaired.fastq $this->{qc_fq2_files}[$i].unpaired.fastq");

		}


	} else #single-end
	{
		for (my $i=0; $i < $this->{num_of_samples}; $i++)
		{
			$this->run_cmd("fastq_quality_filter -Q 33 -q $this->{filter_quality_score} -p $this->{filter_percent} -i $this->{fq1_files}[$i] -o $this->{qc_fq1_files}[$i]"); #option -Q 33 is for Sanger format
		}
	}
	$this->log_it("####################################QC TRIM ENDS####################################");

}


sub qc_filter
{

	my $this = shift;
	$this->log_it("####################################QC FILTER STARTS####################################");
	my @cmds = ();

	if($this->{pe})
	{
		for (my $i=0; $i < $this->{num_of_samples}; $i++)
		{
			push(@cmds, "fastq_quality_filter -Q 33 -q $this->{filter_quality_score} -p $this->{filter_percent} -i $this->{fq1_files}[$i] -o $this->{qc_fq1_files}[$i].tmp.fastq");
			push(@cmds, "fastq_quality_filter -Q 33 -q $this->{filter_quality_score} -p $this->{filter_percent} -i $this->{fq2_files}[$i] -o $this->{qc_fq2_files}[$i].tmp.fastq");
			$this->run_multiple_cmds(@cmds);
			#$this->run_cmd("fastq_quality_filter -Q 33 -q $this->{filter_quality_score} -p $this->{filter_percent} -i $this->{fq1_files}[$i] -o $this->{qc_fq1_files}[$i].tmp.fastq"); #option -Q 33 is for Sanger format
			#$this->run_cmd("fastq_quality_filter -Q 33 -q $this->{filter_quality_score} -p $this->{filter_percent} -i $this->{fq2_files}[$i] -o $this->{qc_fq2_files}[$i].tmp.fastq"); #option -Q 33 is for Sanger format
			$this->run_cmd("$this->{fix_paired_reads} $this->{qc_fq1_files}[$i].tmp.fastq $this->{qc_fq2_files}[$i].tmp.fastq ".
				"$this->{qc_fq1_files}[$i] $this->{qc_fq2_files}[$i] $this->{qc_fq1_files}[$i].unpaired.fastq $this->{qc_fq2_files}[$i].unpaired.fastq");

		}


	} else #single-end
	{
		for (my $i=0; $i < $this->{num_of_samples}; $i++)
		{
			$this->run_cmd("fastq_quality_filter -Q 33 -q $this->{filter_quality_score} -p $this->{filter_percent} -i $this->{fq1_files}[$i] -o $this->{qc_fq1_files}[$i]"); #option -Q 33 is for Sanger format
		}
	}
	$this->log_it("####################################QC FILTER ENDS####################################");

}


sub fastq_qc_summary
{
	my $this = shift;
	my $fastq_filename;
	my $fastqc_data_file;


	if ($this->{dry_run})
	{
		return;
	}



	$this->log_it("####################################FASTQ QC SUMMARY STARTS####################################");
	FileSystem::WriteFile($this->{report_html}, get_header());

	#print "$this->{qc_summary_header}\n";
	#my @header = FileSystem::ReadFile($this->{qc_summary_header});

	my $data =	"<tr><td>FASTQ_FILE</td>\n".
				"<td>TOTAL_READS1</td><td>QC_CHECK1</td><td>TOTAL_READS2</td><td>QC_CHECK2</td>\n".
				"<td>PERCENT_READS</font></td></tr>\n\n";

	my $warning = '<font color=red>WARNING<br>(lower quartile < 20)</font>';
	my $passed = '<font color=green>PASSED</font>';
	my $qc_header;


	if ($this->{trim_or_filter} == 1)
	{
		$qc_header = '<div class="main"><div class="module"><h2 id="M0">QC using simple trim'."</h2></div>\n";
	}elsif ($this->{trim_or_filter} == 2)
	{
		$qc_header = '<div class="main"><div class="module"><h2 id="M0">QC filtering: '.
		$this->{filter_percent}.'% of bases have quality score '.$this->{filter_quality_score}."</h2></div>\n";
	}elsif ($this->{trim_or_filter} == 0)
	{
		$qc_header = '<div class="main"><div class="module"><h2 id="M0">QC using trimmomatic'."</h2></div>\n";
	}


	my $qc_table = q(<table>
<tr>
<th>Fastq file</th>
<th>Before QC<br>Total reads</th>
<th>Before QC<br>Base quality</th>
<th>After QC<br>Total reads</th>
<th>After QC<br>Base quality</th>
<th>% of reads remaining</th>
</tr>
);
	FileSystem::AppendFile($this->{report_html}, $qc_header);
	FileSystem::AppendFile($this->{report_html}, $qc_table);

	for (my $i=0; $i < $this->{num_of_samples}; $i++)
	{
		my $total_reads1;
		my $total_reads2;
		my $pct_reads;
		my $QC_CHECK;
		my $current_data = $data;

		$fastq_filename = FileSystem::GetFileNameFromPath($this->{fq1_files}[$i]);
		$current_data =~ s/FASTQ_FILE/$fastq_filename/;

		$fastq_filename =~ s/\.fastq/_fastqc/;
		$fastqc_data_file = $this->{output_project_dir}."/$this->{before_qc_report_sub_dir}/".$fastq_filename."/fastqc_data.txt";

		($total_reads1, $QC_CHECK) = $this->parse_fastqc_data($fastqc_data_file);
		$current_data =~ s/TOTAL_READS1/$total_reads1/;
		if ($QC_CHECK)
		{
			#print "passed\n";
			$current_data =~ s/QC_CHECK1/$passed/;
		}
		else
		{
			#print "warning\n";
			$current_data =~ s/QC_CHECK1/$warning/;
		}



		$fastq_filename = FileSystem::GetFileNameFromPath($this->{qc_fq1_files}[$i]);
		$fastq_filename =~ s/\.fastq/_fastqc/;
		$fastqc_data_file = $this->{output_project_dir}."/$this->{after_qc_report_sub_dir}/".$fastq_filename."/fastqc_data.txt";

		($total_reads2, $QC_CHECK) = $this->parse_fastqc_data($fastqc_data_file);
		$current_data =~ s/TOTAL_READS2/$total_reads2/;
		if ($QC_CHECK)
		{
			$current_data =~ s/QC_CHECK2/$passed/;
		}
		else
		{
			$current_data =~ s/QC_CHECK2/$warning/;
		}

		$pct_reads = sprintf("%.1f", $total_reads2/$total_reads1*100).'%';
		$current_data =~ s/PERCENT_READS/$pct_reads/;

		FileSystem::AppendFile($this->{report_html}, $current_data);
	}


	if($this->{pe})
	{

		for (my $i=0; $i < $this->{num_of_samples}; $i++)
		{
			my $total_reads1;
			my $total_reads2;
			my $pct_reads;
			my $QC_CHECK;
			my $current_data = $data;

			$fastq_filename = FileSystem::GetFileNameFromPath($this->{fq2_files}[$i]);
			$current_data =~ s/FASTQ_FILE/$fastq_filename/;

			$fastq_filename =~ s/\.fastq/_fastqc/;
			$fastqc_data_file = $this->{output_project_dir}."/$this->{before_qc_report_sub_dir}/".$fastq_filename."/fastqc_data.txt";

			($total_reads1, $QC_CHECK) = $this->parse_fastqc_data($fastqc_data_file);
			$current_data =~ s/TOTAL_READS1/$total_reads1/;
			if ($QC_CHECK)
			{
				#print "passed\n";
				$current_data =~ s/QC_CHECK1/$passed/;
			}
			else
			{
				#print "warning\n";
				$current_data =~ s/QC_CHECK1/$warning/;
			}



			$fastq_filename = FileSystem::GetFileNameFromPath($this->{qc_fq2_files}[$i]);
			$fastq_filename =~ s/\.fastq/_fastqc/;
			$fastqc_data_file = $this->{output_project_dir}."/$this->{after_qc_report_sub_dir}/".$fastq_filename."/fastqc_data.txt";

			($total_reads2, $QC_CHECK) = $this->parse_fastqc_data($fastqc_data_file);
			$current_data =~ s/TOTAL_READS2/$total_reads2/;
			if ($QC_CHECK)
			{
				$current_data =~ s/QC_CHECK2/$passed/;
			}
			else
			{
				$current_data =~ s/QC_CHECK2/$warning/;
			}

			$pct_reads = sprintf("%.1f", $total_reads2/$total_reads1*100).'%';
			$current_data =~ s/PERCENT_READS/$pct_reads/;

			FileSystem::AppendFile($this->{report_html}, $current_data);
		}
	}
	FileSystem::AppendFile($this->{report_html}, '</table>');

	my $tail = '</div><div class="footer">Produced by Center for Research Informatics, UIC (zlei2@uic.edu)</div></body></html>';
	FileSystem::AppendFile($this->{report_html}, $tail);
	$this->log_it("####################################FASTQ QC SUMMARY ENDS####################################");
}

sub parse_fastqc_data
{
	my $this = shift;
	my $fastqc_data_file = shift;

	my @output = ();
	my @fastqc_data = FileSystem::ReadFile($fastqc_data_file);
	chomp(@fastqc_data);
	my ($TotalSequences_string, $total_reads) = split(/\t/, $fastqc_data[6]);

	my $QC_CHECK = 1;
	my $j = 0;
	my $k = 0;
	for (; $j<=$#fastqc_data; $j++)
	{
		if ($fastqc_data[$j] =~ /Base	Mean	Median	Lower Quartile	Upper Quartile	10th Percentile	90th Percentile/)
		{
			last;
		}
	}
	for ($k=$j+1; $k<=$#fastqc_data; $k++)
	{
		if ($fastqc_data[$k] =~ /END_MODULE/)
		{
			last;
		}
		my ($base_pos, $mean, $median, $low_quantile) = split(/\t/, $fastqc_data[$k]);
		if ($low_quantile < 20)
		{
			$QC_CHECK = 0;
			last;
		}
	}
	return ($total_reads, $QC_CHECK);
}

########################################################################################################
# 2. QC step: qcstats
# use fastx tool kit to process
# http://hannonlab.cshl.edu/fastx_toolkit/commandline.html
# check linux version by "uname -a"
# To set PATH=$PATH:$HOME/bin:$HOME/tools/fastx_tk
# vi .bash_profile
# source .bash_profile
#sub qc
#{
#	my $this = shift;
#	my $qc_stat_file;
#	my $qc_txt_file;
#	my $qc_png_file;
#
#	for (my $i=0; $i < $this->{num_of_samples}; $i++)
#	{
#		$qc_stat_file = $this->{fq1_files}[$i];
#		$qc_stat_file =~ s/$this->{output_project_dir}/$this->{output_project_dir}\/$this->{qc_sub_dir}/;
#		$qc_txt_file = $qc_stat_file."_QC.txt";
#		$qc_qs_png_file = $qc_stat_file."_QC_quality.score.png";
#		$qc_nu_png_file = $qc_stat_file."_QC_nucleotide.dist.png";
#		$this->run_cmd("fastx_quality_stats -Q 33 -i $this->{fq1_files}[$i] -o $qc_txt_file"); #option -Q 33 is for Sanger format
#		$this->run_cmd("fastq_quality_boxplot_graph.sh -i $qc_txt_file -o $qc_qs_png_file -t \"$this->{fq1_files}[$i]\"");
#		$this->run_cmd("fastx_nucleotide_distribution_graph.sh -i $qc_txt_file -o $qc_nu_png_file -t \"$this->{fq1_files}[$i] nucleotide distribution\"");
#	}
#
#	if($this->{pe})
#	{
#		for (my $i=0; $i < $this->{num_of_samples}; $i++)
#		{
#			$qc_stat_file = $this->{fq2_files}[$i];
#			$qc_stat_file =~ s/$this->{output_project_dir}/$this->{output_project_dir}\/$this->{qc_sub_dir}/;
#			$qc_txt_file = $qc_stat_file."_QC.txt";
#			$qc_qs_png_file = $qc_stat_file."_QC_quality.score.png";
#			$qc_nu_png_file = $qc_stat_file."_QC_nucleotide.dist.png";
#			$this->run_cmd("fastx_quality_stats -Q 33 -i $this->{fq2_files}[$i] -o $qc_txt_file"); #option -Q 33 is for Sanger format
#			$this->run_cmd("fastq_quality_boxplot_graph.sh -i $qc_txt_file -o $qc_qs_png_file -t \"$this->{fq2_files}[$i]\"");
#			$this->run_cmd("fastx_nucleotide_distribution_graph.sh -i $qc_txt_file -o $qc_nu_png_file -t \"$this->{fq2_files}[$i] nucleotide distribution\"");
#		}
#	}
#}


#sub fastq_trimmer
#{
#	my $this = shift;
#	my ($first_base_pos, $last_base_pos) = @_;
#	print "reads at positions($first_base_pos, $last_base_pos) selected\n";
#	my $trimmed_fq1 = $this->{sanger_fq1};
#	$trimmed_fq1 =~ s/fq$/trimmed.fq/;
#	$this->run_cmd("fastx_trimmer -Q 33 -f $first_base_pos -l $last_base_pos -i $this->{sanger_fq1} -o $trimmed_fq1");
#	$this->{sanger_fq1} = $trimmed_fq1;
#
#	if($this->{pe}) 
#	{
#		my $trimmed_fq2 = $this->{sanger_fq2};
#		$trimmed_fq2 =~ s/fq$/trimmed.fq/;
#		$this->run_cmd("fastx_trimmer -Q 33 -f $first_base_pos -l $last_base_pos -i $this->{sanger_fq2} -o $trimmed_fq2");
#		$this->{sanger_fq2} = $trimmed_fq2;
#	}
#}


#http://qualimap.bioinfo.cipf.es
# /crihome/zlei2/tools/qualimap_v0.7.1/qualimap --java-mem-size=8G bamqc -bam ILM_002_Feinstein_MSF1_1_GTAGAG_L001.dedup.bam -outdir ILM_002_Feinstein_MSF1_1_GTAGAG_L001.dedup
# /crihome/zlei2/tools/qualimap_v0.7.1/qualimap --java-mem-size=8G bamqc -bam ILM_002_Feinstein_MSF1_1_GTAGAG_L001.bam -outdir ILM_002_Feinstein_MSF1_1_GTAGAG_L001


# /crihome/zlei2/tools/qualimap_v0.7.1/qualimap --java-mem-size=8G bamqc -bam ILM_002_Feinstein_MSF_4_GTCCGC_L002.dedup.bam -outdir ILM_002_Feinstein_MSF_4_GTCCGC_L002.dedup
# /crihome/zlei2/tools/qualimap_v0.7.1/qualimap --java-mem-size=8G bamqc -bam ILM_002_Feinstein_MSF_4_GTCCGC_L002.bam -outdir ILM_002_Feinstein_MSF_4_GTCCGC_L002


sub qualimap
{
	my $this = shift;
	my $bam_file = shift;

	my $out_dir = $bam_file;
	$out_dir =~ s/\.bam//;

	my $QM_cmd = "$this->{qualimap_path} --java-mem-size=8G bamqc -bam $bam_file -outdir $out_dir";
	if ($this->{target_bed} ne "")
	{
		$QM_cmd .= " -gff $this->{target_bed}";
	}
	$this->run_cmd($QM_cmd);
}


sub bam_qc_summary
{
	my $this = shift;
	my $bam_file;
	my $dedup_bam_file;
	my $realigned_bam_file;
	my $realigned_sorted_fixed_bam_file;
	my $realigned_recal_bam_file;
	my $reduced_bam_file;
	$this->log_it("####################################BAM FILE QC STARTS####################################");

	for (my $i=0; $i < $this->{num_of_samples}; $i++)
	{
		$bam_file = $this->{bam_dir}."/".$this->{samples}[$i].".bam";
		$aligned_bam_file = $this->{bam_dir}."/".$this->{samples}[$i].".aligned.bam";
		$realigned_bam_file = $this->{realigned_bam_dir}."/".$this->{samples}[$i].".realigned.bam";
		$realigned_sorted_fixed_bam_file = $this->{realigned_bam_dir}."/".$this->{samples}[$i].".realigned.sorted.fixed.bam";
		$realigned_recal_bam_file = $this->{realigned_bam_dir}."/".$this->{samples}[$i].".realigned.recal.bam";
		$realigned_recal_bam_file = $this->{realigned_bam_dir}."/".$this->{samples}[$i].".realigned.recal.bam";
		$reduced_bam_file = $this->{realigned_bam_dir}."/".$this->{samples}[$i].".reduced.bam";

		#my @bam_files = ($bam_file, $dedup_bam_file, $realigned_bam_file, $realigned_sorted_fixed_bam_file, $realigned_recal_bam_file, $reduced_bam_file);
		my @bam_files = ($bam_file, $aligned_bam_file);

		foreach my $bam (@bam_files)
		{
			if (-e $bam)
			{
				$this->qualimap($bam);
			}
		}
	}
	$this->log_it("####################################BAM FILE QC ENDS####################################");

}


#This function 'bam_summary' is obselete
sub bam_summary
{
	my $this = shift;

	my $bam_file;
	my $dedup_bam_file;
	my $realigned_bam_file;
	my $realigned_sorted_fixed_bam_file;
	my $realigned_recal_bam_file;
	my $reduced_bam_file;

	my $data =	"<tr><td>_BAM_FILE</td><td>_TOTAL_READS</td>\n".
			"<td>_MAPPED_READS</td><td>_UNMAPPED_READS</td><td>_PERCENT_MAPPED_READS</td><td>_UNIQUE_MAPPED_READS</td>\n".
			"<td>_PERCENT__UNIQUE_MAPPED_READS</font></td></tr>\n\n";


	my $bam_header = '<div class="module"><h2 id="M1">Bam file(s) statistics</h2></div>'."\n";
	my $bam_table = q(<table>
<tr>
<th>Bam file</th>
<th>Total reads</th>
<th>Mapped reads</th>
<th>Unmapped reads</th>
<th>% mapped reads</th>
<th>Unique mapped reads</th>
<th>% unique mapped reads</th>

</tr>
);

	FileSystem::AppendFile($this->{report_html}, $bam_header);
	FileSystem::AppendFile($this->{report_html}, $bam_table);


	for (my $i=0; $i < $this->{num_of_samples}; $i++)
	{
		$bam_file = $this->{bam_dir}."/".$this->{samples}[$i].".bam";
		$dedup_bam_file = $this->{bam_dir}."/".$this->{samples}[$i].".dedup.bam";
		$realigned_bam_file = $this->{realigned_bam_dir}."/".$this->{samples}[$i].".realigned.bam";
		$realigned_sorted_fixed_bam_file = $this->{realigned_bam_dir}."/".$this->{samples}[$i].".realigned.sorted.fixed.bam";
		$realigned_recal_bam_file = $this->{realigned_bam_dir}."/".$this->{samples}[$i].".realigned.recal.bam";
		$realigned_recal_bam_file = $this->{realigned_bam_dir}."/".$this->{samples}[$i].".realigned.recal.bam";
		$reduced_bam_file = $this->{realigned_bam_dir}."/".$this->{samples}[$i].".reduced.bam";

		#my @bam_files = ($bam_file, $dedup_bam_file, $realigned_bam_file, $realigned_sorted_fixed_bam_file, $realigned_recal_bam_file, $reduced_bam_file);
		my @bam_files = ($bam_file, $dedup_bam_file);

		foreach my $bam (@bam_files)
		{
			if (-e $bam)
			{
				my $current_data = $data;
				my ($total_reads, $mapped_reads, $unmapped_reads, $unique_mapped_reads) = $this->bam_stat($bam);
				my $pct_mapped_reads = sprintf("%.1f", $mapped_reads/$total_reads*100).'%';
				my $pct_unique_mapped_reads = sprintf("%.1f", $unique_mapped_reads/$total_reads*100).'%';

				$bam = FileSystem::GetFileNameFromPath($bam);

				$current_data =~ s/_BAM_FILE/$bam/;
				$current_data =~ s/_TOTAL_READS/$total_reads/;
				$current_data =~ s/_MAPPED_READS/$mapped_reads/;
				$current_data =~ s/_UNMAPPED_READS/$unmapped_reads/;
				$current_data =~ s/_PERCENT_MAPPED_READS/$pct_mapped_reads/;
				$current_data =~ s/_UNIQUE_MAPPED_READS/$unique_mapped_reads/;
				$current_data =~ s/_PERCENT__UNIQUE_MAPPED_READS/$pct_unique_mapped_reads/;

				FileSystem::AppendFile($this->{report_html}, $current_data);
			}
		}
	}
}

#This function 'bam_stat' is obselete
sub bam_stat
{
	my $this = shift;
	my $bam_file = shift;
	$this->log_it("####################################$bam_file STAT STARTS####################################");

	my $total_reads = FileSystem::GetCmdRes("samtools view -c $bam_file");
	$this->log_it("Total reads: $total_reads");

	my $mapped_reads = FileSystem::GetCmdRes("samtools view  -F0x4 -c $bam_file");
	$this->log_it("Mapped reads: $mapped_reads");

	#my $unmapped_reads = FileSystem::GetCmdRes("samtools view  -f0x4 -c $bam_file");
	my $unmapped_reads = $total_reads - $mapped_reads;
	$this->log_it("Unmapped reads: $unmapped_reads");

	my $unique_mapped_reads = FileSystem::GetCmdRes("samtools view $bam_file | grep XT:A:U | wc -l");
	$this->log_it("Unique mapped reads: $unique_mapped_reads");

	#my $multiple_mapped_reads = FileSystem::GetCmdRes("samtools view $bam_file | grep XT:A:R | wc -l");
	#$this->log_it("Multiple mapped reads: $multiple_mapped_reads");

	$this->log_it("####################################$bam_file STAT ENDS####################################");

	return ($total_reads, $mapped_reads, $unmapped_reads, $unique_mapped_reads);

}



sub bwa
{
	my $this = shift;
	my $fastq_file;
	my $sai_file1;
	my $sai_file2;
	my $sam_file;
	my $bam_file;
	my $aligned_bam_file;

	my @cmds = ();

	$this->log_it("####################################ALIGNEMENT STARTS####################################");

	for (my $i=0; $i < $this->{num_of_samples}; $i++)
	{
		@cmds = ();
		$fastq_file1 = $this->{qc_fq1_files}[$i];
		$sai_file1 = $fastq_file1;
		$sai_file1 =~ s/$this->{output_project_dir}\/$this->{qc_fastq_data_sub_dir}/$this->{bam_dir}/;
		$sai_file1 .= ".sai";



		$sam_file = $this->{bam_dir}."/".$this->{samples}[$i].".sam";
		$bam_file = $this->{bam_dir}."/".$this->{samples}[$i].".bam";
		$aligned_bam_file = $this->{bam_dir}."/".$this->{samples}[$i].".aligned.bam";

		if (!$this->{remove_duplicates})
		{
			$bam_file = $aligned_bam_file;
		}


		if($this->{pe})  ########## paired end ##################
		{
			$fastq_file2 = $this->{qc_fq2_files}[$i];
			$sai_file2 = $fastq_file2;
			$sai_file2 =~ s/$this->{output_project_dir}\/$this->{qc_fastq_data_sub_dir}/$this->{bam_dir}/;
			$sai_file2 .= ".sai";

			push(@cmds, "bwa aln -t $this->{cpu_threads} $this->{reference_fa} $fastq_file1 > $sai_file1");
			push(@cmds, "bwa aln -t $this->{cpu_threads} $this->{reference_fa} $fastq_file2 > $sai_file2");
			#$this->run_cmd("bwa aln -t $this->{cpu_threads} $this->{reference_fa} $fastq_file2 > $sai_file2");
			$this->run_multiple_cmds(@cmds);


			#$this->run_cmd("bwa sampe $this->{reference_fa} $sai_file1 $sai_file2 $fastq_file1 $fastq_file2 > $sam_file"); #pipe | to bam
			$this->run_cmd("bwa sampe -r \"\@RG\tID:$this->{samples}[$i]\tLB:$this->{samples}[$i]\tSM:$this->{samples}[$i]\tPL:$this->{platform}\" ".
				"$this->{reference_fa} $sai_file1 $sai_file2 $fastq_file1 $fastq_file2 > $sam_file"); #pipe | to bam

			#add read group
			#bwa sampe -P -r "@RG\tID:sample\tLB:chr_WGS\tSM:sample\tPL:ILLUMINA" /bwa/indexes/hg19 A_R1.sai A_R2.sai A_R1.fastq A_R2.fastq > A.sam
		}
		else ########## single end ##################
		{

			if ($this->{long_reads})
			{
				$this->run_cmd("bwa mem -R \"\@RG\tID:$this->{samples}[$i]\tLB:$this->{samples}[$i]\tSM:$this->{samples}[$i]\tPL:$this->{platform}\" ".
				"-t $this->{cpu_threads} $this->{reference_fa} $fastq_file1 > $sam_file");
			}
			else
			{

				$this->run_cmd("bwa aln -t $this->{cpu_threads} $this->{reference_fa} $fastq_file1 > $sai_file1");
				$this->run_cmd("bwa samse -r \"\@RG\tID:$this->{samples}[$i]\tLB:$this->{samples}[$i]\tSM:$this->{samples}[$i]\tPL:$this->{platform}\" ".
					"$this->{reference_fa} $sai_file1 $fastq_file1 > $sam_file");
			}

		}


		$this->run_cmd("$this->{picard_SortSam} SO=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true ".
			"I=$sam_file ".
			"O=$bam_file");
		$this->run_cmd("samtools flagstat $bam_file > $bam_file.stat");


		#my $sort_rg_bam_file = $this->{bam_dir}."/".$this->{samples}[$i].".sort.rg.bam";
		#$this->run_cmd("$this->{picard_AddOrReplaceReadGroups} RGID=$this->{samples}[$i] RGLB=$this->{samples}[$i] RGPL=$this->{platform} RGPU=run RGSM=$this->{samples}[$i] CREATE_INDEX=TRUE VALIDATION_STRINGENCY=SILENT ".
		#"I=$sort_bam_file ".
		#"O=$sort_rg_bam_file");

		#my $sort_rg_bai_file = $sort_rg_bam_file;
		#$sort_rg_bai_file =~ s/m$/i/;

		#my $sort_rg_dedup_bam_file = $this->{output_project_dir}."/".$this->{bam_sub_dir}."/".$this->{samples}[$i]."sort.rg.dedup.bam";
		#if (lc($this->{platform}) ne "iontorrent")
		#{
		#}

		if ($this->{remove_duplicates})
		{
			$this->run_cmd("$this->{picard_MarkDuplicates} VALIDATION_STRINGENCY=SILENT M=$bam_file.PCR_duplicates.txt REMOVE_DUPLICATES=true AS=true ".
			"I=$bam_file ".
			"O=$aligned_bam_file");
			$this->run_cmd("samtools flagstat $aligned_bam_file > $aligned_bam_file.stat");
		}


		###################################3
		# Is the following step necessary??
		if($this->{pe})
		{
			$this->run_cmd("$this->{picard_FixMateInformation} SO=coordinate VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=2000000 CREATE_INDEX=true ".
			"I=$aligned_bam_file"); #The output file to write to. If no output file is supplied, the input file is overwritten.
		}

		my $output_bai_file = $aligned_bam_file;
		$output_bai_file =~ s/m$/i/;
		$this->run_cmd("$this->{picard_BuildBamIndex} I=$aligned_bam_file O=$output_bai_file VALIDATION_STRINGENCY=SILENT"); #use picard

		#unlink($sai_file1, $sai_file2, $sam_file); #$sort_bam_file, $sort_rg_bam_file, $sort_rg_bai_file
		$this->log_it("DELETE the following files if they exist: $sai_file1; $sai_file2; $sam_file");#; $sort_bam_file; $sort_rg_bam_file; $sort_rg_bai_file
	}
	$this->log_it("OUTPUT of BAM file(s) (sorted, read group added, PCR duplicates removed (except iontorrent)) are in directory ".$this->{bam_dir});
	$this->log_it("####################################ALIGNEMENT ENDS####################################");

#	#//sed 's/[ \t]*$//g' $sample_name.sam >$sample_name.sed.sam
	#scp -r ./picard-tools-1.44/ NUSSTF\\gmslz@172.25.138.143:~/tools/
	#add $HOME/tools/picard-tools-1.44 to PATH in file .bash_profile 
	#source .bash_profile
	##### remove PCR duplicate using PICARD
	##### remove PRC duplicates by picard, use AS=true so that picard knows the bam is sorted.
}



#Tool	Full name				Type of traversal	NT	NCT	SG
#RTC	RealignerTargetCreator	RodWalker			+	-	-
#IR		IndelRealigner			ReadWalker			-	-	+
#BR		BaseRecalibrator		LocusWalker			-	+	+
#PR		PrintReads				ReadWalker			-	+	-
#RR		ReduceReads				ReadWalker			-	-	+
#UG		UnifiedGenotyper		LocusWalker			+	+	+
sub tmp
{
	my $this = shift;
	my $bam_file;
	my $realigned_bam_file;
	my $realigned_bam_stat_file;

	my $realigned_sorted_bam_file;
	my $realigned_sorted_fixed_bam_file;

	$this->log_it("####################################GATK REALIGNMENT STARTS####################################");


	for (my $i=0; $i < $this->{num_of_samples}; $i++)
	{
		$bam_file = $this->{bam_dir}."/".$this->{samples}[$i].".dedup.bam";
		#$bam_file = $this->{bam_dir}."/".$this->{samples}[$i].".bam";
		$realigned_bam_file = $this->{realigned_bam_dir}."/".$this->{samples}[$i].".realigned.bam";
		$realigned_sorted_bam_file = $this->{realigned_bam_dir}."/".$this->{samples}[$i].".realigned.sorted.bam";
		$realigned_sorted_fixed_bam_file = $this->{realigned_bam_dir}."/".$this->{samples}[$i].".realigned.sorted.fixed.bam";


		#step 3 --- picard FixMateInformation; Fixing the mate pairs of realigned reads #paired end data only?
		$this->run_cmd("$this->{picard_FixMateInformation} SO=coordinate VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=2000000 CREATE_INDEX=true ".
		"I=$realigned_sorted_bam_file ".
		"O=$realigned_sorted_fixed_bam_file");
		#my $realigned_sorted_fixed_bai_file = $realigned_sorted_fixed_bam_file;
		#$bai_file =~ s/m$/i/;
		#$this->run_cmd("$this->{picard_BuildBamIndex} I=$realigned_sorted_fixed_bam_file O=$realigned_sorted_fixed_bai_file VALIDATION_STRINGENCY=SILENT"); #use picard

		$this->run_cmd("samtools flagstat $realigned_sorted_fixed_bam_file > $realigned_sorted_fixed_bam_file.stat");

	}
	$this->log_it("####################################GATK REALIGNMENT ENDS####################################");

}


sub gatk_realign
{
	my $this = shift;
	my $bam_file;
	my $realigned_bam_file;
	my $realigned_bam_stat_file;

	my $realigned_sorted_bam_file;
	my $realigned_sorted_fixed_bam_file;

	$this->log_it("####################################GATK REALIGNMENT STARTS####################################");


	for (my $i=0; $i < $this->{num_of_samples}; $i++)
	{
		$bam_file = $this->{bam_dir}."/".$this->{samples}[$i].".aligned.bam";
		#$bam_file = $this->{bam_dir}."/".$this->{samples}[$i].".bam";
		$realigned_bam_file = $this->{realigned_bam_dir}."/".$this->{samples}[$i].".realigned.bam";
		$realigned_sorted_bam_file = $this->{realigned_bam_dir}."/".$this->{samples}[$i].".realigned.sorted.bam";
		#$realigned_sorted_fixed_bam_file = $this->{realigned_bam_dir}."/".$this->{samples}[$i].".realigned.sorted.fixed.bam";

		if ($this->{user_interval_file} ne "") 
		{
			$this->run_cmd("$this->{gatk_GenomeAnalysisTK} -T IndelRealigner -l INFO ". #not support -nt nor -nct
			"-R $this->{reference_fa} ".
			"-I $bam_file -targetIntervals $this->{user_interval_file} ".
			"-o $realigned_bam_file -stats $this->{realigned_bam_dir}/$this->{samples}[$i].realigned.stats.txt > $this->{log_dir}/$this->{samples}[$i].realigned.log");
			##$this->run_cmd("$this->{picard_BuildBamIndex} I=$realigned_bam_file O=$realigned_bam_file.bai"); #No need, bai file will be generated in the step 1.2

		}
		else
		{
			#step 1 --- Local realignment around indels
			#step 1.1 --- gatk -T RealignerTargetCreator; Determining (small) suspicious intervals around indels
			my $RTC_cmd = "$this->{gatk_GenomeAnalysisTK} -T RealignerTargetCreator -l INFO -nt $this->{data_threads} ".
					"-R $this->{reference_fa} -known $this->{Mills_KG_indel_vcf} ";
			if ($this->{target_bed} ne "")
			{
				$RTC_cmd = $RTC_cmd."-L $this->{target_bed} ";
			}
			$RTC_cmd = $RTC_cmd."-I $bam_file ".
					"-o $this->{realigned_bam_dir}/$this->{samples}[$i].intervals > $this->{log_dir}/$this->{samples}[$i].intervals.log";

			#step 1.2 --- gatk -T IndelRealigner; Running the realigner over those intervals
			my $IR_cmd = "$this->{gatk_GenomeAnalysisTK} -T IndelRealigner -l INFO ". #not support -nt nor -nct
				"-R $this->{reference_fa} -known $this->{Mills_KG_indel_vcf} ".
				"-I $bam_file -targetIntervals $this->{realigned_bam_dir}/$this->{samples}[$i].intervals ".
				"-o $realigned_bam_file -stats $this->{realigned_bam_dir}/$this->{samples}[$i].realigned.stats.txt > $this->{log_dir}/$this->{samples}[$i].realigned.log";
			$this->run_cmd($RTC_cmd);
			$this->run_cmd($IR_cmd);
			##$this->run_cmd("$this->{picard_BuildBamIndex} I=$realigned_bam_file O=$realigned_bam_file.bai"); #No need, bai file will be generated in the step 1.2
		}
		#step 2 --- picard SortSam
		$this->run_cmd("$this->{picard_SortSam} SO=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true ".
			"I=$realigned_bam_file ".
			"O=$realigned_sorted_bam_file");
		#my $realigned_sorted_bai_file = $realigned_sorted_bam_file;
		#$bai_file =~ s/m$/i/;
		#$this->run_cmd("$this->{picard_BuildBamIndex} I=$realigned_sorted_bam_file O=$realigned_sorted_bai_file VALIDATION_STRINGENCY=SILENT"); #use picard

		#step 3 --- picard FixMateInformation; Fixing the mate pairs of realigned reads #paired end data only?
		#$this->run_cmd("$this->{picard_FixMateInformation} SO=coordinate VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=2000000 CREATE_INDEX=true ".
		#"I=$realigned_sorted_bam_file ".
		#"O=$realigned_sorted_fixed_bam_file");
		#my $realigned_sorted_fixed_bai_file = $realigned_sorted_fixed_bam_file;
		#$bai_file =~ s/m$/i/;
		#$this->run_cmd("$this->{picard_BuildBamIndex} I=$realigned_sorted_fixed_bam_file O=$realigned_sorted_fixed_bai_file VALIDATION_STRINGENCY=SILENT"); #use picard

		#$this->run_cmd("samtools flagstat $realigned_sorted_bam_file > $realigned_sorted_bam_file.stat");

	}
	$this->log_it("####################################GATK REALIGNMENT ENDS####################################");

}

#base quality score recalibration
sub gatk_bqsr
{
	#DBSNP: ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.5/hg19/dbsnp_137.hg19.vcf.gz

	my $this = shift;
	my $realigned_bam_file;
	my $grp_file;
	my $realigned_recal_bam_file;

	#if (lc($this->{platform}) eq "iontorrent")
	#{
	#	for (my $i=0; $i < $this->{num_of_samples}; $i++)
	#	{
	#		$realigned_bam_file = $this->{realigned_bam_dir}."/".$this->{samples}[$i].".realigned.sorted.fixed.bam";
	#		$realigned_recal_bam_file = $this->{realigned_bam_dir}."/".$this->{samples}[$i].".realigned.recal.bam";
	#
	#		$this->log_it("####################################BASE QUALITY SCORE RECALIBRATION FOR ION TORRENT DATA --RENAME FILE ####################################");
	#		if (-e $realigned_bam_file) 
	#		{
	#			$this->run_cmd("mv $realigned_bam_file $realigned_recal_bam_file");
#
#			}
#
#			my $bai_file = $realigned_recal_bam_file;
#			$bai_file =~ s/m$/i/;
#			$this->run_cmd("$this->{picard_BuildBamIndex} I=$realigned_recal_bam_file O=$bai_file VALIDATION_STRINGENCY=SILENT"); #use picard
#		}
#		return;
#	}

 
	$this->log_it("####################################BASE QUALITY SCORE RECALIBRATION STARTS####################################");
	for (my $i=0; $i < $this->{num_of_samples}; $i++)
	{
		$realigned_bam_file = $this->{realigned_bam_dir}."/".$this->{samples}[$i].".realigned.sorted.bam";
		$grp_file = $this->{realigned_bam_dir}."/".$this->{samples}[$i].".CovarTable_Recal.grp";
		$realigned_recal_bam_file = $this->{realigned_bam_dir}."/".$this->{samples}[$i].".realigned.recal.bam";


		my $BR_cmd = "$this->{gatk_GenomeAnalysisTK} -T BaseRecalibrator -l INFO -nct $this->{cpu_threads} ". #
		"-R $this->{reference_fa} -knownSites $this->{dbsnp_vcf} ".

		if ($this->{target_bed} ne "")
		{
			$BR_cmd = $BR_cmd."-L $this->{target_bed} ";
		}
		$BR_cmd = $BR_cmd."-I $realigned_bam_file ".
			"-o $grp_file ";
		#http://gatkforums.broadinstitute.org/discussion/1241/dataprocessingpipeline-and-disable_indel_quals
		#this issue only affected GATK Lite, 
		#BaseRecalibrator does indeed recalibrate base qualities for indels in all versions - See more at: http://gatkforums.broadinstitute.org/discussion/1241/dataprocessingpipeline-and-disable_indel_quals#sthash.zMJKTzrv.dpuf
		if ($this->{gatk_full_version} == 0)
		{
			$BR_cmd .= "--disable_indel_quals";
		}
		$this->run_cmd($BR_cmd);

		$this->run_cmd("$this->{gatk_GenomeAnalysisTK} -T PrintReads -l INFO -nct $this->{cpu_threads} ". #not support 
		"-R $this->{reference_fa} ".
 		"-I $realigned_bam_file ".
		"-BQSR $grp_file ".
		"-o $realigned_recal_bam_file");

		$this->run_cmd("samtools flagstat $realigned_recal_bam_file > $realigned_recal_bam_file.stat");
	}
	$this->log_it("####################################BASE QUALITY SCORE RECALIBRATION ENDS####################################");
}
#java -jar GenomeAnalysisTK.jar \ -T PrintReads \ -R reference.fasta \ -I input.bam \ -BQSR recalibration_report.grp \ -o output.bam 
#- See more at: http://gatkforums.broadinstitute.org/discussion/44/base-quality-score-recalibration-bqsr#sthash.kUSxcTcq.dpuf


#GATK Lite does not support all of the features of the full version: the ReduceReads walker is available only in the full version of the GATK
sub gatk_reduce_bam
{
	my $this = shift;
	my $realigned_recal_bam_file;
	my $reduced_bam_file;

	if ($this->{gatk_full_version} == 0)
	{
		$this->log_it("####################################REDUCING BAM IS NOT SUPPORTED BY GATK LITE VERSION####################################");
		return;
	}

	$this->log_it("####################################REDUCING READS STARTS####################################");
	for (my $i=0; $i < $this->{num_of_samples}; $i++)
	{
		$realigned_recal_bam_file = $this->{realigned_bam_dir}."/".$this->{samples}[$i].".realigned.recal.bam";
		$reduced_bam_file = $this->{realigned_bam_dir}."/".$this->{samples}[$i].".reduced.bam";
		$this->run_cmd("$this->{gatk_GenomeAnalysisTK} -T ReduceReads -l INFO ". #not support -nt nor -nct
		"-R $this->{reference_fa} ".
 		"-I $realigned_recal_bam_file ".
		"-o $reduced_bam_file");

		my $reduced_bam_stat_file = $this->{realigned_bam_dir}."/".$this->{samples}[$i].".reduced.bam.stat";
		$this->run_cmd("samtools flagstat $reduced_bam_file > $reduced_bam_stat_file");
	}


	$this->log_it("####################################REDUCING READS ENDS####################################");

}

#http://gatkforums.broadinstitute.org/discussion/1530/how-to-use-unifiedgenotyper-annotation-option
#Re: your annotation question, the UG uses a certain number of annotations by default. The core set are called "standard annotations" and are used by default by all tools (unless otherwise specified by the 'exclude' argument).
#Currently, the standard annotations are the following:
#BaseQualityRankSumTest
#ChromosomeCounts
#Coverage   #DepthOfCoverage
#DepthPerAlleleBySample
#FisherStrand
#HaplotypeScore
#InbreedingCoeff
#MappingQualityRankSumTest
#MappingQualityZero
#QualByDepth
#ReadPosRankSumTest
#RMSMappingQuality
#SpanningDeletions
#TandemRepeatAnnotator

#http://www.icbi.at/svnsimplex/CommonBioCommands/tags/simplex-1.0/CommonBioCommands/clusterservice/cheops/scripts/gatk/snp_calling.sh
#-A                         annotation: comma separated list of specific annotations to apply to variant calls.
#  -T UnifiedGenotyper -A AlleleBalance -l INFO -nct 
#!!!!!!!!!!!!!!!!!The DepthOfCoverage annotation has been renamed to Coverage, because of a naming conflict with the DepthOfCoverage walker. - See more at: http://gatkforums.broadinstitute.org/discussion/2318/undocumented-change-in-2-4-a-depthofcoverage#sthash.lN0omXQB.dpuf

#$this->{genometyper_annotation} = -A HaplotypeScore,DepthOfCoverage,FisherStrand,ReadPosRankSumTest,MappingQualityRankSumTest,QualByDepth
# -T UnifiedGenotyper -glm BOTH \ -A Coverage \ -A AlleleBalance \ -A HaplotypeScore \ -A HomopolymerRun \ -A MappingQualityZero \ -A QualByDepth \ -A RMSMappingQuality \ -A SpanningDeletions \ -A MappingQualityRankSumTest \ -A ReadPosRankSumTest \ -A FisherStrand \ -A InbreedingCoeff - See more at: http://gatkforums.broadinstitute.org/discussion/2325/unifiedgenotyper-core-dump#sthash.SNKyZN1R.dpuf
sub gatk_raw_snp_indel
{
	my $this = shift;
	my $realigned_recal_bam_file;
	my $reduced_bam_file;
	my $raw_snp_indel_vcf = $this->{output_project_dir}."/$this->{snp_sub_dir}/raw_snp_indel.vcf";
	print "******$raw_snp_indel_vcf****\n";

	$this->log_it("####################################RAW SNP CALLING STARTS####################################");

	my $sample_list = "";

	for (my $i=0; $i < $this->{num_of_samples}; $i++)
	{
		if ($this->{gatk_full_version} == 0)
		{
			$sample_list .= "-I ".$this->{realigned_bam_dir}."/".$this->{samples}[$i].".realigned.recal.bam ";
		}
		else
		{
			$sample_list .= "-I ".$this->{realigned_bam_dir}."/".$this->{samples}[$i].".reduced.bam ";
		}
	}
	$this->log_it("####################################call gatk####################################");

	if ($this->{snp_caller} == 1)
	{
		$this->run_cmd("$this->{gatk_GenomeAnalysisTK} -T UnifiedGenotyper -l INFO -nt $this->{data_threads} -nct $this->{cpu_threads} ". #not suport 
			"-R $this->{reference_fa} ".
			"--dbsnp $this->{dbsnp_vcf} ".
			"-stand_call_conf $this->{stand_call_conf} -stand_emit_conf $this->{stand_emit_conf} ".
			"-dcov $this->{dcov} -glm $this->{glm} ".
			"-A HaplotypeScore -A Coverage -A FisherStrand -A ReadPosRankSumTest -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality ".
			$sample_list.
			"-o $raw_snp_indel_vcf");

	}elsif ($this->{snp_caller} == 2)
	{

		$this->run_cmd("$this->{gatk_GenomeAnalysisTK} -T HaplotypeCaller -l INFO ". #not support -nt nor -nct
			"-R $this->{reference_fa} ".
			"--dbsnp $this->{dbsnp_vcf} ".
			"-stand_call_conf $this->{stand_call_conf} -stand_emit_conf $this->{stand_emit_conf} ".
			$sample_list.
			"-o $raw_snp_indel_vcf");
	}

	$this->log_it("####################################RAW SNP CALLING ENDS####################################");
}

#variant quality score recalibration
# All filtering criteria are learned from the data itself.
#http://www.broadinstitute.org/gatk/guide/article?id=1259
#SNP specific recommendations
#For SNPs we use both HapMap v3.3 and the Omni chip array from the 1000 Genomes Project as training data. 
#In addition we take the highest confidence SNPs from the project's callset. 
#These datasets are available in the GATK resource bundle. Arguments for VariantRecalibrator command:
#   -percentBad 0.01 -minNumBad 1000 \
#   -resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.b37.sites.vcf \
#   -resource:omni,known=false,training=true,truth=true,prior=12.0 1000G_omni2.5.b37.sites.vcf \
#   -resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.vcf \
#   -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp.b37.vcf \
#   -an QD -an MQRankSum -an ReadPosRankSum -an FS -an DP \
#   -mode SNP \

#Note that, for the above to work, the input vcf needs to be annotated with the corresponding values (QD, FS, DP, etc.). 
#If any of these values are somehow missing, then VariantAnnotator needs to be run first so that VariantRecalibration can run properly.

#In our testing we've found that in order to achieve the best exome results one needs to use an exome SNP and/or indel callset with at least 30 samples. 
#For users with experiments containing fewer exome samples there are several options to explore:

#(1) Add additional samples for variant calling, either by sequencing additional samples or using publicly available exome bams from the 1000 Genomes Project 
#(this option is used by the Broad exome production pipeline) OR
#(2) Use the VQSR with the smaller variant callset but experiment with the precise argument settings 
#(try adding --maxGaussians 4 --percentBad 0.05 to your command line, for example)
#- See more at: http://gatkforums.broadinstitute.org/discussion/1259/what-vqsr-training-sets-arguments-should-i-use-for-my-specific-project#sthash.biLuT8Ud.dpuf

sub hard_or_soft_filtering
{
	my $this = shift;

	if ($this->{hard_filtering})
	{
		$this->hard_filtering();
	}
	else
	{
		$this->gatk_vqsr();

		$this->select_variant_by_VQSLOD();
	}
}

#http://www.broadinstitute.org/gatk/guide/article?id=1259
sub gatk_vqsr
{
	my $this = shift;
	my $raw_snp_indel_vcf = $this->{output_project_dir}."/$this->{snp_sub_dir}/raw_snp_indel.vcf";
	my $snp_recal_file		= $this->{output_project_dir}."/$this->{snp_sub_dir}/snp.output.recal";
	my $snp_tranches_file	= $this->{output_project_dir}."/$this->{snp_sub_dir}/snp.output.tranches";
	my $snp_rscript_file	= $this->{output_project_dir}."/$this->{snp_sub_dir}/snp.output.plots.R";

	my $indel_recal_file		= $this->{output_project_dir}."/$this->{snp_sub_dir}/indel.output.recal";
	my $indel_tranches_file	= $this->{output_project_dir}."/$this->{snp_sub_dir}/indel.output.tranches";
	my $indel_rscript_file	= $this->{output_project_dir}."/$this->{snp_sub_dir}/indel.output.plots.R";

	#snp.model <- BuildErrorModelWithVQSR(raw.vcf, SNP)
    #indel.model <- BuildErrorModelWithVQSR(raw.vcf, INDEL)
    #recalibratedSNPs.rawIndels.vcf <- ApplyRecalibration(raw.vcf, snp.model, SNP)
    #analysisReady.vcf <- ApplyRecalibration(recalibratedSNPs.rawIndels.vcf, indel.model, INDEL)

	my $recalibratedSNPs_rawIndels_vcf = $this->{output_project_dir}."/$this->{snp_sub_dir}/recalibratedSNPs_rawIndels.vcf";
	my $recalibrated_filtered_vcf = $this->{output_project_dir}."/$this->{snp_sub_dir}/RECALIBRATED.FILTERED.vcf";


	$this->log_it("####################################VARIANT QUALITY RECALIBRATION STARTS####################################");


	#InbreedingCoeff is a population level statistic that requires at least 10 samples in order to be calculated.

	#Depth of coverage (the DP annotation invoked by Coverage) should not be used when working with hybrid capture datasets 
	#since there is extreme variation in the depth to which targets are captured! 
	#In whole genome experiments this variation is indicative of error but that is not the case in capture experiments. 

	#The UnifiedGenotyper produces a statistic called the HaplotypeScore which should be used for SNPs.
	#This statistic isn't necessary for the HaplotypeCaller - See more at: http://gatkforums.broadinstitute.org/discussion/1259/what-vqsr-training-sets-arguments-should-i-use-for-my-specific-project#sthash.biLuT8Ud.dpuf

	if ($this->{snp_caller} == 1) #UnifiedGenotyper
	{
		if ($this->{snp_model_annotation} !~ /HaplotypeScore/)
		{
			$this->{snp_model_annotation} .= " -an HaplotypeScore";
		}
		#if ($this->{indel_model_annotation} !~ /HaplotypeScore/) #in 2.4 HaplotypeScore is no longer applied to indel 
		#{
		#	$this->{indel_model_annotation} .= " -an HaplotypeScore";
		#}
	}
	my $VR_snp_model_cmd ="$this->{gatk_GenomeAnalysisTK} -T VariantRecalibrator -l INFO -nt $this->{data_threads} ". #support nt
				"-R $this->{reference_fa} ".
				"$this->{snp_model_para} ".
				"-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $this->{hapmap_vcf} ".
				"-resource:omni,known=false,training=true,truth=false,prior=12.0 $this->{KG_omni_vcf} ".
				"-resource:1000G,known=false,training=true,truth=false,prior=10.0 $this->{KG_snp_vcf} ".
				"-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $this->{dbsnp_vcf} ".
				"$this->{snp_model_annotation} ". # -an InbreedingCoeff for sample size >= 10
				"-mode SNP ".
				"-input $raw_snp_indel_vcf ".
				"-recalFile $snp_recal_file -tranchesFile $snp_tranches_file -rscriptFile $snp_rscript_file ";
	$this->run_cmd($VR_snp_model_cmd);

	my $VR_indel_model_cmd ="$this->{gatk_GenomeAnalysisTK} -T VariantRecalibrator -l INFO -nt $this->{data_threads} ".
				"-R $this->{reference_fa} ".
				"$this->{indel_model_para} ".
				"-resource:mills,known=false,training=true,truth=true,prior=12.0 $this->{Mills_KG_indel_vcf} ".
				"-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $this->{dbsnp_vcf} ".
				"$this->{indel_model_annotation} ". # -an InbreedingCoeff for sample size >= 10
				"-mode INDEL ".
				"-input $raw_snp_indel_vcf ".
				"-recalFile $indel_recal_file -tranchesFile $indel_tranches_file -rscriptFile $indel_rscript_file";
	$this->run_cmd($VR_indel_model_cmd);




	$this->run_cmd("$this->{gatk_GenomeAnalysisTK} -T ApplyRecalibration -l INFO -nt $this->{data_threads} ". #support nt
				"-R $this->{reference_fa} ".
				"--ts_filter_level $this->{snp_filter_level} ".
				"-mode SNP ".
				"-recalFile $snp_recal_file -tranchesFile $snp_tranches_file ".
				"-input $raw_snp_indel_vcf ".
				"-o $recalibratedSNPs_rawIndels_vcf");

	$this->run_cmd("$this->{gatk_GenomeAnalysisTK} -T ApplyRecalibration -l INFO -nt $this->{data_threads} ". 
				"-R $this->{reference_fa} ".
				"--ts_filter_level $this->{indel_filter_level} ".
				"-mode INDEL ".
				"-recalFile $indel_recal_file -tranchesFile $indel_tranches_file ".
				"-input $recalibratedSNPs_rawIndels_vcf ".
				"-o $recalibrated_filtered_vcf");

	$this->log_it("####################################VARIANT QUALITY RECALIBRATION ENDS####################################");
}

#rpoplin Posts: 96GSA Official Member mod
#May 14
#Hey Manny,
#I'm glad that you are enjoying the VQSR. For most of our projects we are generally interested in operating in the regime where sensitivity is at or above 99% so we are happy to tolerate variants in which VQSLOD is near zero in order to get there. We use the VQSLOD simply as a way to rank-order the variants and then decide where to apply that cutoff. If you really want to be sure that your variants have a very low false positive rate then VQSLOD > 3 is a great place to start.
#I would love to see somebody including the VQSLOD in their statistical methods. I think there is a lot of untapped potential there.
#I hope that is helpful. Let me know if you have any other questions.
#Cheers,
#- See more at: http://gatkforums.broadinstitute.org/discussion/2639/vqsr-vqslod-and-indels#sthash.i71hDtql.dpuf


sub select_variant_by_VQSLOD
{

	my $this = shift;
	my $recalibrated_filtered_vcf = $this->{output_project_dir}."/$this->{snp_sub_dir}/RECALIBRATED.FILTERED.vcf";
	my $VQSLOD_vcf = $this->{output_project_dir}."/$this->{snp_sub_dir}/$this->{sample_name}.VQSLOD".$this->{VQSLOD_Threshold}.'.vcf';

	if ($this->{dry_run})
	{
		$this->log_it("Select variants from ".$recalibrated_filtered_vcf." to ".$VQSLOD_vcf);
		return;
	}


	$this->log_it("####################################VARIANT SELECTION BY VQSLOD STARTS####################################");


#	my @header = ();
#	my $i = 0;
	open(SRC_VCF, $recalibrated_filtered_vcf) || die("$!: $recalibrated_filtered_vcf: Make sure the file path is correct!!\n");
	while (my $line = <SRC_VCF>)
	{
		if ($line =~ /^\#/)
		{
				FileSystem::WriteFile($VQSLOD_vcf, $line);
		}
		else
		{	
			if ($line =~ /;VQSLOD=(.*);/)
			{
				my $VQSLOD=$1;
				#print "$VQSLOD\n";
				if ($VQSLOD >= $this->{VQSLOD_Threshold})
				{
					FileSystem::WriteFile($VQSLOD_vcf, $line);
				}
			}
		}
	}
	close(SRC_VCF);
	FileSystem::Close();
	$this->log_it("####################################VARIANT SELECTION BY VQSLOD ENDS####################################");
}


sub split_into_homo_hetero
{
	my $this = shift;

	my $VQSLOD_vcf = $this->{output_project_dir}."/$this->{snp_sub_dir}/$this->{sample_name}.VQSLOD".$this->{VQSLOD_Threshold}.'.vcf';

	if ($this->{hard_filtering})
	{
		$VQSLOD_vcf = $this->{output_project_dir}."/$this->{snp_sub_dir}/HARD.FILETERING.SNP.INDEL.vcf";

	}

	my $homozygous_vcf = $VQSLOD_vcf;
	my $heterozygous_vcf = $VQSLOD_vcf;
	my @header = ();
	my $line;
	$homozygous_vcf =~ s/vcf$/homozygous.vcf/;
	$heterozygous_vcf =~ s/vcf$/heterozygous.vcf/;

	if ($this->{dry_run})
	{
		$this->log_it("Split ".$$VQSLOD_vcf." into ".$homozygous_vcf."and ".$heterozygous_vcf);
		return;
	}


	$this->log_it("####################################HOMO HETERO VCF STARTS####################################");



	open(SRC_VCF, $VQSLOD_vcf) || die("$!: $VQSLOD_vcf: Make sure the file path is correct!!\n");
	while ($line = <SRC_VCF>)
	{
			if ($line =~ /^\#/)
			{
					push(@header, $line);
			}
			else
			{
					last;
			}
	}

	FileSystem::WriteFile($homozygous_vcf, @header);
	FileSystem::WriteFile($heterozygous_vcf, @header);

	my @fields = split(/\t/, $line);
	if ($fields[9] =~ /^1[\/\|]{1}1/)
	{
		FileSystem::WriteFile($homozygous_vcf, $line);
	}
	else
	{
		FileSystem::WriteFile($heterozygous_vcf, $line);
	}
	while ($line = <SRC_VCF>)
	{
		@fields = split(/\t/, $line);
		if ($fields[9] =~ /^1[\/\|]{1}1/)
		{
			FileSystem::WriteFile($homozygous_vcf, $line);
		}
		else
		{
			FileSystem::WriteFile($heterozygous_vcf, $line);
		}
	}
	FileSystem::Close();

	$this->log_it("####################################HOMO HETERO VCF ENDS####################################");

}




#sudo easy_install poster
#./seattle_seq_annotation.py snp MSF1_MSF4_intersect.vcf zlei2@uic.edu
#./seattle_seq_annotation.py indel MSF1_MSF4_intersect.vcf zlei2@uic.edu

###############################3
#nohup /crihome/zlei2/tools/ZL_seattle_seq_annotation.py both ILM_002_Feinstein_MSF_4_GTCCGC_L002.VQSLOD4.homozygous.vcf zlei2@uic.edu >homo.both.log &
#nohup /crihome/zlei2/tools/ZL_seattle_seq_annotation.py both interset_of_MSF1_EXOME_and_MSF4_EXOME.heterozygous.vcf zlei2@uic.edu >hetero.both.log &
#nohup /crihome/zlei2/tools/ZL_seattle_seq_annotation.py both interset_of_MSF1_EXOME_and_MSF4_EXOME.homozygous.vcf zlei2@uic.edu >homo.both.log &


sub seattle_seq_annotation
{
	my $this = shift;
	my $VQSLOD_vcf = $this->{output_project_dir}."/$this->{snp_sub_dir}/$this->{sample_name}.VQSLOD".$this->{VQSLOD_Threshold}.'.vcf';

	if ($this->{hard_filtering})
	{
		$VQSLOD_vcf = $this->{output_project_dir}."/$this->{snp_sub_dir}/HARD.FILETERING.SNP.INDEL.vcf";

	}

	$this->log_it("####################################SEATTLE SEQ STARTS####################################");

	#my $recalibrated_filtered_vcf = $this->{output_project_dir}."/$this->{snp_sub_dir}/recalibrated.filtered.$this->{pipeline_time_stamp}.vcf";
	#my $email = 'zlei2@uic.edu';

	#$this->run_cmd("$this->{seattle_path} snp $recalibrated_filtered_vcf $this->{email}");
	#$this->run_cmd("$this->{seattle_path} indel $recalibrated_filtered_vcf $this->{email}");

	$this->run_cmd("$this->{seattle_path} both $VQSLOD_vcf $this->{email}");
	$this->log_it("####################################SEATTLE SEQ ENDS####################################");
}




#MSF1
#java -Xmx4g -Djava.io.tmpdir=/tmp/ -jar /data1/rhel60/gatk_git20130205/dist/GenomeAnalysisTK.jar -T SelectVariants 
#-R /cri/pinal/reference/hg19_ucsc_broad/gatk_resources/ucsc_bwa/ucsc.hg19.fasta
#-V /crihome/zlei2/Feinstein_MSF1/snp/recalibrated.filtered.20130616_205555.vcf 
#-select "VQSLOD >= 4"
#-o MSF1_VQSLOD4.vcf



#java -Xmx4g -Djava.io.tmpdir=/tmp/ -jar /data1/rhel60/gatk_git20130205/dist/GenomeAnalysisTK.jar -T SelectVariants -R /cri/pinal/reference/hg19_ucsc_broad/gatk_resources/ucsc_bwa/ucsc.hg19.fasta -V /crihome/zlei2/Feinstein_MSF1/snp/recalibrated.filtered.20130616_205555.vcf -select "VQSLOD>4" -o MSF1_VQSLOD4.vcf




#MSF4

#java -Xmx4g -Djava.io.tmpdir=/tmp/ -jar /data1/rhel60/gatk_git20130205/dist/GenomeAnalysisTK.jar -T SelectVariants 
#-R /cri/pinal/reference/igenome_ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta_minusChrM/genome_minusChrM.fa 
#--variant /crihome/zlei2/Feinstein_MSF4/snp/recalibrated.filtered.20130616_205555.vcf 
#--select_expressions "VQSLOD >= 4"
#-o MSF4_VQSLOD4.vcf







#Hi Jeremy,
#
#In our study we used SelectVariant to extract variant position with VQSLOD > 4.0.
#
#$pathToGATK \
#  -T SelectVariants \
#  -R reference \
#  --variant inVCF \
#  --select_expressions "VQSLOD >= $lodCutoff" \
#  --out outName.vcf


#"java -Xmx4g -Djava.io.tmpdir=/tmp/ -jar /data1/rhel60/gatk_git20130205/dist/GenomeAnalysisTK.jar -T SelectVariants -R /cri/pinal/reference/igenome_ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta_minusChrM/genome_minusChrM.fa -V:variant union_MSF1_MSF4.vcf -select 'set == ""Intersection"";' -o intersect_MSF1_MSF4.vcf"



#http://gatkforums.broadinstitute.org/discussion/2806/howto-apply-hard-filters-to-a-call-set
sub hard_filtering
{

#http://gatkforums.broadinstitute.org/discussion/2806/howto-apply-hard-filters-to-a-call-set
	my $this = shift;
	my $raw_snp_indel_vcf = $this->{output_project_dir}."/$this->{snp_sub_dir}/raw_snp_indel.vcf";
	my $raw_snp_vcf	  = $this->{output_project_dir}."/$this->{snp_sub_dir}/raw_snp.vcf";
	my $raw_indel_vcf = $this->{output_project_dir}."/$this->{snp_sub_dir}/raw_indel.vcf";

	my $filtered_snp_vcf = $this->{output_project_dir}."/$this->{snp_sub_dir}/HARD.FILTERING.SNP.vcf";
	my $filtered_indel_vcf = $this->{output_project_dir}."/$this->{snp_sub_dir}/HARD.FILTERING.INDEL.vcf";
	my $filtered_snp_indel_vcf = $this->{output_project_dir}."/$this->{snp_sub_dir}/HARD.FILETERING.SNP.INDEL.vcf";

	$this->log_it("####################################HARD FILETERING STARTS####################################");

	$this->run_cmd("$this->{gatk_GenomeAnalysisTK} -T SelectVariants -l INFO ". #not suport -nt $this->{data_threads} 
		"-R $this->{reference_fa} -selectType SNP ".
		"-V $raw_snp_indel_vcf ".
		"-o $raw_snp_vcf");

	$this->run_cmd("$this->{gatk_GenomeAnalysisTK} -T SelectVariants -l INFO ". #not suport -nt $this->{data_threads} 
		"-R $this->{reference_fa} -selectType INDEL ".
		"-V $raw_snp_indel_vcf ".
		"-o $raw_indel_vcf");


#http://gatkforums.broadinstitute.org/discussion/2806/howto-apply-hard-filters-to-a-call-set

	$this->run_cmd("$this->{gatk_GenomeAnalysisTK} -T VariantFiltration -l INFO ". #not suport -nt $this->{data_threads} 
		"-R $this->{reference_fa} ".
		"-V $raw_snp_vcf ".
		"$this->{snp_filter} ".
		"-o $filtered_snp_vcf");
	
	$this->run_cmd("$this->{gatk_GenomeAnalysisTK} -T VariantFiltration -l INFO ". #not suport -nt $this->{data_threads} 
		"-R $this->{reference_fa} ".
		"-V $raw_indel_vcf ".
		"$this->{indel_filter} ".
		"-o $filtered_indel_vcf");

	$this->run_cmd("$this->{gatk_GenomeAnalysisTK} -T CombineVariants -l INFO ". #not suport -nt $this->{data_threads} 
		"-R $this->{reference_fa} ".
		"--variant $filtered_snp_vcf ".
		"--variant $filtered_indel_vcf ".
		"-o $filtered_snp_indel_vcf ".
		"-genotypeMergeOptions UNIQUIFY");


	$this->log_it("####################################HARD FILETERING ENDS####################################");



}
		

#################################################################################################################
##
#java -jar GenomeAnalysisTK.jar \ 
#    -T SelectVariants \ 
#    -R reference.fa \ 
#    -V raw_HC_variants.vcf \ 
#    -selectType INDEL \ 
#    -o raw_indels.vcf 
#6. Apply the filter to the Indel call set
#Run the following GATK command:
#java -jar GenomeAnalysisTK.jar \ 
#    -T VariantFiltration \ 
#    -R reference.fa \ 
#    -V raw_indels.vcf \ 
#    --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \ 
#   --filterName "my_indel_filter" \ 
#    -o filtered_indels.vcf 











#java -jar ./dist/GenomeAnalysisTK.jar -T VariantFiltration ref.fasta -o out.vcf -B:variant,VCF input.vcf 
#--filterExpression "QUAL<30.0" --filterName "LowQual"
#--filterExpression "SB>=-1.0" --filterName "StrandBias"
#--filterExpression "QD<1.0" --filterName "QualByDepth"
#--filterExpression "(MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1))" --filterName "HARD_TO_VALIDATE"
#--filterExpression "HRun>=15" --filterName "HomopolymerRun"!

#Example filter
#QD < 2.0 
#ReadPosRankSum < -8.0
#MQ < 40.0
#FS > 60.0 
#MQRandkSum < -12.5
#HaplotypeScore > 13.0
sub variant_filter #for small dataset, please use this hard filtering instead instead of vqsr
{


	$this->run_cmd("$this->{gatk_GenomeAnalysisTK} -T VariantFiltration -l INFO ". #not suport -nt $this->{data_threads} 
		'--clusterSize 3 --clusterWindowSize 10 '.
		'--filterExpression "DP <= 20" --filterName "Depth" '.
		'--filterExpression "SB > -0.10" --filterName "StrandBias" '.
		'--filterExpression "HRun > 8" --filterName "HRun8" '.
		'--filterExpression "QD < 5.0" --filterName "QD5" '.
		'--filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" --filterName "hard2validate" '.
		"-D $this->{dbsnp_rod} -R $this->{reference_fa} ".
		"-B:variant,VCF ./snp_analysis_gatk/$this->{sample_name}.raw.snps.vcf ".
		"-B:mask,Bed ./snp_analysis_gatk/$this->{sample_name}.indels.mask.bed --maskName InDel ".
		"-o ./snp_analysis_gatk/$this->{sample_name}.snp.filtered.vcf > ./log/$this->{sample_name}.snp.filter.log");

}



## filter poor quality & suspicious SNP calls 
#vcftools --vcf foo.gatk.VariantFiltration.snps.vcf \
#--remove-filtered  DP8 --remove-filtered StrandBias --remove-filtered  LowQual \
#--remove-filtered  hard2validate --remove-filtered  SnpCluster \
#--keep-INFO AC --keep-INFO AF --keep-INFO AN --keep-INFO DB \
#--keep-INFO DP --keep-INFO DS --keep-INFO Dels --keep-INFO HRun --keep-INFO HaplotypeScore --keep-INFO MQ --keep-INFO MQ0 --keep-INFO QD --keep-INFO SB --out foo.gatk.good.snps ;








sub coverage
{
		$this->run_cmd("$this->{gatk_GenomeAnalysisTK} -T DepthOfCoverage -l INFO ". #not support -nt $this->{data_threads} 
			"-D $this->{dbsnp_rod} -R $this->{reference_fa} -L $this->{target_bed} ".
			"-I $this->{aligned_bam_file} ".
			"-o ./snp/$this->{sample_name}.depthofcoverage > ./log/$this->{sample_name}.depthofcoverage.log");

}


		#	-[arg]:REFSEQ /crihome/zlei2/knownsites/refSeq.sorted.txt
# java -Xmx2g -Djava.io.tmpdir=/tmp/ -jar /data1/rhel60/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -T DepthOfCoverage -l INFO -nt 4 
#	 -R /cri/pinal/reference/hg19_ucsc_broad/gatk_resources/ucsc_bwa/ucsc.hg19.fasta 
#	 --omitDepthOutputAtEachBase --omitLocusTable --omitIntervalStatistics --includeDeletions 
#	 --summaryCoverageThreshold 10 --summaryCoverageThreshold 30 --summaryCoverageThreshold 50 
#	 --interval_merging OVERLAPPING_ONLY -geneList /crihome/zlei2/knownsites/refSeq.sorted.txt
#	 -I /crihome/zlei2/green_KADKOL_HALO_10/bam/green_KADKOL_HALO_10.dedup.bam 
#	 -o /crihome/zlei2/green_KADKOL_HALO_10/bam/DepthOfCoverage.out.txt 

# java -Xmx2g -Djava.io.tmpdir=/tmp/ -jar /data1/rhel60/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -T DepthOfCoverage -l INFO -nt 4 
#	 -R /cri/pinal/reference/hg19_ucsc_broad/gatk_resources/ucsc_bwa/ucsc.hg19.fasta 
#	 --omitDepthOutputAtEachBase --omitLocusTable --omitIntervalStatistics --includeDeletions 
#	 --summaryCoverageThreshold 10 --summaryCoverageThreshold 30 --summaryCoverageThreshold 50 
#	 --interval_merging OVERLAPPING_ONLY -geneList /crihome/zlei2/knownsites/refSeq.sorted.txt
#	 -I /crihome/zlei2/Feinstein_MSF1/bam/ILM_002_Feinstein_MSF1_1_GTAGAG_L001.dedup.bam
#	 -o /crihome/zlei2/Feinstein_MSF1/bam/DepthOfCoverage.out.txt 

#java -Xmx2g -Djava.io.tmpdir=/tmp/ -jar /data1/rhel60/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -T DepthOfCoverage -l INFO -nt 4 -R /cri/pinal/reference/hg19_ucsc_broad/gatk_resources/ucsc_bwa/ucsc.hg19.fasta --omitDepthOutputAtEachBase --omitLocusTable --omitIntervalStatistics --includeDeletions --summaryCoverageThreshold 10 --summaryCoverageThreshold 30 --summaryCoverageThreshold 50 --interval_merging OVERLAPPING_ONLY -geneList /crihome/zlei2/knownsites/refSeq.sorted.txt -I /crihome/zlei2/Feinstein_MSF1/bam/ILM_002_Feinstein_MSF1_1_GTAGAG_L001.dedup.bam -o /crihome/zlei2/Feinstein_MSF1/bam/DepthOfCoverage.out.txt 


# java -Xmx4g -Djava.io.tmpdir=/tmp/ -jar /data1/rhel60/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -T DepthOfCoverage -l INFO -nt 4 
#	 -R /cri/pinal/reference/hg19_ucsc_broad/gatk_resources/ucsc_bwa/ucsc.hg19.fasta 
#	 --omitDepthOutputAtEachBase --omitLocusTable --omitIntervalStatistics --includeDeletions 
#	 --summaryCoverageThreshold 10 --summaryCoverageThreshold 30 --summaryCoverageThreshold 50 
#	 --interval_merging OVERLAPPING_ONLY -geneList /crihome/zlei2/knownsites/refSeq.sorted.txt
#	 -I /crihome/zlei2/Feinstein_MSF4/bam/ILM_002_Feinstein_MSF_4_GTCCGC_L002.dedup.bam
#	 -o /crihome/zlei2/Feinstein_MSF4/bam/DepthOfCoverage.out.txt 

#nohup java -Xmx4g -Djava.io.tmpdir=/tmp/ -jar /data1/rhel60/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -T DepthOfCoverage -l INFO -nt 4 -R /cri/pinal/reference/hg19_ucsc_broad/gatk_resources/ucsc_bwa/ucsc.hg19.fasta --omitDepthOutputAtEachBase --omitLocusTable --omitIntervalStatistics --includeDeletions --summaryCoverageThreshold 10 --summaryCoverageThreshold 30 --summaryCoverageThreshold 50 --interval_merging OVERLAPPING_ONLY -geneList /crihome/zlei2/knownsites/refSeq.sorted.txt -I /crihome/zlei2/Feinstein_MSF4/bam/ILM_002_Feinstein_MSF_4_GTCCGC_L002.dedup.bam -o /crihome/zlei2/Feinstein_MSF4/bam/DepthOfCoverage.out.txt &




sub snp_analysis_gatk
{
#	GATK Documents:
#	http://www.broadinstitute.org/gatk/gatkdocs/
#	e.g.:
#	http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_filters_DuplicateReadFilter.html
#	http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_filters_VariantFiltration.html

	my $this = shift;  
	#chdir("../snp_analysis_gtk");
	my $max_depth = 10000;

#	goto EVAL1;
	$this->run_cmd("$this->{gatk_GenomeAnalysisTK} -T DepthOfCoverage -l INFO ". #not support -nt $this->{data_threads} 
			"-D $this->{dbsnp_rod} -R $this->{reference_fa} -L $this->{target_bed} ".
			"-I $this->{aligned_bam_file} ".
			"-o ./snp_analysis_gatk/$this->{sample_name}.depthofcoverage > ./log/$this->{sample_name}.depthofcoverage.log");


	#step 1 --- Local realignment around indels
	#step 1.1 --- gatk -T RealignerTargetCreator; Determining (small) suspicious intervals around indels
	$this->run_cmd("$this->{gatk_GenomeAnalysisTK} -T RealignerTargetCreator -l INFO ". #not support -nt $this->{data_threads} 
			"-D $this->{dbsnp_rod} -R $this->{reference_fa} ".
			"-I $this->{aligned_bam_file} ".
			"-o ./snp_analysis_gatk/$this->{sample_name}.intervals > ./log/$this->{sample_name}.intervals.log");
	#step 1.2 --- gatk -T IndelRealigner; Running the realigner over those intervals
	$this->run_cmd("$this->{gatk_GenomeAnalysisTK} -T IndelRealigner -l INFO ". #not support -nt $this->{data_threads} 
		"-D $this->{dbsnp_rod} -R $this->{reference_fa} ".
		"-I $this->{aligned_bam_file} -targetIntervals ./snp_analysis_gatk/$this->{sample_name}.intervals ".
		"-o ./bam/$this->{sample_name}.realigned.bam -stats ./bam/$this->{sample_name}.realign.stats.txt > ./log/$this->{sample_name}.realign.log");
	$this->run_cmd("samtools index ./bam/$this->{sample_name}.realigned.bam");
	#step 1.3 --- picard SortSam
	$this->run_cmd("$this->{picard_SortSam} SO=coordinate VALIDATION_STRINGENCY=SILENT ".
		"I=./bam/$this->{sample_name}.realigned.bam ".
		"O=./bam/$this->{sample_name}.realigned.sorted.bam");
	$this->run_cmd("samtools index ./bam/$this->{sample_name}.realigned.sorted.bam");
	#step 1.4 --- picard FixMateInformation; Fixing the mate pairs of realigned reads #paired end data only?
	$this->run_cmd("$this->{picard_FixMateInformation} SO=coordinate VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=2000000 ".
		"I=./bam/$this->{sample_name}.realigned.sorted.bam ".
		"O=./bam/$this->{sample_name}.realigned.sorted.fixed.bam");
	$this->run_cmd("samtools index ./bam/$this->{sample_name}.realigned.sorted.fixed.bam"); #use picard
	#step 1.5 --- skip picard MarkDuplicates


	#step 2 --- Recalibrate quality scores --> Final Recalibrated Bam #skip for ion torrent data
	#step 2.1 --- gatk -T CountCovariates;
	$this->run_cmd("$this->{gatk_GenomeAnalysisTK} -T CountCovariates -l INFO -nt $this->{data_threads} ".
		"--downsample_to_coverage $max_depth --default_platform Illumina ".
		"-cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate ".
		"-D $this->{dbsnp_rod} -R $this->{reference_fa} ".
 		"-I ./bam/$this->{sample_name}.realigned.sorted.fixed.bam ".
		"-recalFile ./snp_analysis_gatk/$this->{sample_name}.recal.csv > ./log/$this->{sample_name}.count.covariates.log");
	#step 2.2 --- gatk -T TableRecalibration  --> Final Recalibrated Bam
	$this->run_cmd("$this->{gatk_GenomeAnalysisTK} -T TableRecalibration -l INFO --default_platform Illumina ". #not support -nt $this->{data_threads} 
		"-R $this->{reference_fa} ".
 		"-I ./bam/$this->{sample_name}.realigned.sorted.fixed.bam -recalFile ./snp_analysis_gatk/$this->{sample_name}.recal.csv ".
		"--out $this->{final_recalibrated_bam_file} > ./log/$this->{sample_name}.recalibration.log");
	$this->run_cmd("samtools index $this->{final_recalibrated_bam_file}");
	# skip 2.3 -T CountCovariates
	# skip 2.4 AnalyzeCovariate

#############################################################################
#	You may merge multiple final_recalibrated_bam files before Genotyping	#
#############################################################################


	#step 3 --- Genotyping
	#step 3.1 --- gatk -T IndelGenotyperV2;
	#print STDERR "\nRunning the indel genotyper:\n";
    $this->run_cmd("$this->{gatk_GenomeAnalysisTK} -T IndelGenotyperV2 -l INFO --window_size 450 ". #not suport -nt $this->{data_threads}
		"-D $this->{dbsnp_rod} -R $this->{reference_fa} ".
		"-I $this->{final_recalibrated_bam_file} ".
		"-o ./snp_analysis_gatk/$this->{sample_name}.indels.vcf -bed ./snp_analysis_gatk/$this->{sample_name}.indels.bed > ./log/$this->{sample_name}.indels.log");

	# skip filter indels: perl /tools/filterSingleSampleCalls.pl
	#step 3.2 python/makeIndelMask.py
	# create intervals to mask snvs close to indels
	# The output of the IndelGenotyper is used to mask out SNPs near indels. 
	# The number 10 in this command stands for the number of bases that will be included on either side of the indel
	$this->run_cmd("python $this->{gatk_dir}/python/makeIndelMask.py ./snp_analysis_gatk/$this->{sample_name}.indels.bed 10 ./snp_analysis_gatk/$this->{sample_name}.indels.mask.bed");
	#step 3.3 --- gatk -T Unifiedgenotyper;
	$this->run_cmd("$this->{gatk_GenomeAnalysisTK} -T UnifiedGenotyper -l INFO ". #not suport -nt $this->{data_threads} 
		"-stand_call_conf 30 -stand_emit_conf 10 -pl SOLEXA -mmq 30 -mm40 3 ".
		"-D $this->{dbsnp_rod} -R $this->{reference_fa} ".
		"-I $this->{final_recalibrated_bam_file} ".
		"-o ./snp_analysis_gatk/$this->{sample_name}.raw.snps.vcf > ./log/$this->{sample_name}.snp.log");

	#step 3.4 --- gatk -T VariantFiltration; VariantFiltration is used to annotate suspicious calls from VCF files based on their failing given filters. 
	$this->run_cmd("$this->{gatk_GenomeAnalysisTK} -T VariantFiltration -l INFO ". #not suport -nt $this->{data_threads} 
		'--clusterSize 3 --clusterWindowSize 10 '.
		'--filterExpression "DP <= 8" --filterName "DP8" '.
		'--filterExpression "SB > -0.10" --filterName "StrandBias" '.
		'--filterExpression "HRun > 8" --filterName "HRun8" '.
		'--filterExpression "QD < 5.0" --filterName "QD5" '.
		'--filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" --filterName "hard2validate" '.
		"-D $this->{dbsnp_rod} -R $this->{reference_fa} ".
		"-B:variant,VCF ./snp_analysis_gatk/$this->{sample_name}.raw.snps.vcf ".
		"-B:mask,Bed ./snp_analysis_gatk/$this->{sample_name}.indels.mask.bed --maskName InDel ".
		"-o ./snp_analysis_gatk/$this->{sample_name}.snp.filtered.vcf > ./log/$this->{sample_name}.snp.filter.log");

## grep -i 'indel'  Testis_T28.snp.filtered.vcf|more
## grep -i 'pass'  Testis_T28.snp.filtered.vcf|more

#EVAL1:
	$this->run_cmd("$this->{gatk_GenomeAnalysisTK} -T VariantEval -l INFO ". #not suport -nt $this->{data_threads} 
		"-D $this->{dbsnp_rod} -R $this->{reference_fa} -L $this->{target_bed} ".
		#-G Indicates that your VCF file has genotypes; -A Print extended evaluation information; -V Print a list of interesting sites (FPs, FNs, etc).
		"-B:eval,VCF ./snp_analysis_gatk/$this->{sample_name}.snp.filtered.vcf ".
		"-o ./snp_analysis_gatk/$this->{sample_name}.snp.filtered.eval > ./log/$this->{sample_name}.snp.filter.eval.log");



#annotate variants Annovar 
#$cmd = "annotate_variation.pl -filter -batchsize 50m -dbtype snp$verdbsnp -buildver $buildver -outfile $outfile $queryfile $dbloc";


# Consensus Quality: Phred-scaled likelihood that the genotype is wrong
# SNP Quality: Phred-scaled likelihood that the genotype is identical to the reference,
#	which is also called `SNP quality'. Suppose the reference base is A and in alignment
#	we see 17 G and 3 A. We will get a low consensus quality because it is difficult to
#	distinguish an A/G heterozygote from a G/G homozygote. We will get a high SNP
#	quality, though, because the evidence of a SNP is very strong.
# RMS: root mean square (RMS) mapping quality, a measure of the variance of quality scores
# Coverage: reads covering the position
}

sub run_multiple_cmds
{
	my $this = shift;
	my @cmds = @_;
#	print "$this->{cluster_mode}        ***\n";
#	exit;

	if ($this->{dry_run})
	{
		foreach my $cmd (@cmds)
		{
			$this->log_it($cmd);
		}
		return;
	}


	if ($this->{cluster_mode})
	{
		$this->parallel_run_cmds(@cmds);
	}
	else
	{
		foreach my $cmd (@cmds)
		{
			$this->run_cmd($cmd);
		}
	}
}


sub parallel_run_cmds
{

	my $this = shift;
	my @cmds = @_;
	my $qsub_script_template = 	q(#!/bin/sh
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -m bea
#PBS -M EMAIL

cd PWD
CMD_TO_SUMMITTED

);
	my $current_script;
	my $script_file;
	my $pwd = FileSystem::GetCmdRes("pwd");
	$qsub_script_template =~ s/PWD/$pwd/;
	$qsub_script_template =~ s/EMAIL/$this->{email}/;

	my %running_job_ids = ();
	my $num_of_remaining_jobs = 0;

	my $start_time = localtime();
	my ($year, $mon, $day, $hour, $min, $sec) = ($start_time->year + 1900, $start_time->mon+1, $start_time->mday, $start_time->hour, $start_time->min, $start_time->sec);
	my $start_time_sec = timelocal_nocheck($sec, $min, $hour, $day, $mon-1, $year);
	#printf("%04d/%02d/%02d %02d:%02d:%02d", $year, $mon, $day, $hour, $min, $sec);
	#print("\tSTART @cmds\n");

	for (my $i=0; $i<=$#cmds; $i++)
	{
		my $cmd = $cmds[$i];
		$this->log_it("PARALLEL START $cmd");
		$current_script = $qsub_script_template;
		$current_script =~ s/CMD_TO_SUMMITTED/$cmd/;
		$script_file = "$this->{output_project_dir}/qsub_script.$i.sh";
		FileSystem::WriteFile($script_file, $current_script);
		FileSystem::Close();
		my $jobid = FileSystem::GetCmdRes("/usr/local/bin/qsub $script_file"); #better use full path of qsub

#		print "qsub $script_file\n";
#		print "job id = $jobid\n";
#		exit;

		if ($jobid =~ /^(\d+)\./)
		{
			$running_job_ids{$1} = 1;
			$this->log_it("qsub job id = $jobid");

		} else
		{
			$this->log_it("qsub failed for $script_file");
		}
		sleep(1);
	}

	sleep(20);

	while (1)
	{
		$num_of_remaining_jobs = 0;
		open(PIPE,'/usr/local/bin/qstat -i |') or die "qstat failed, make sure use the full path of qstat\n";   
		my @bjobs_list1=<PIPE>;           
		close PIPE;  

		open(PIPE,'/usr/local/bin/qstat -r |') or die "qstat failed, make sure use the full path of qstat\n";   
		my @bjobs_list2=<PIPE>;           
		close PIPE;  

		my @bjobs_list = (@bjobs_list1, @bjobs_list2);

		for (my $i=1; $i<=$#bjobs_list; $i++)
		{
			if ($bjobs_list[$i] =~ /^(\d+)\./)
			{
				my $jobid = $1;

				#$this->log_it("qstat $jobid");

				if (exists($running_job_ids{$jobid}))
				{
					$num_of_remaining_jobs++;
				}
			}
		}
		if ($num_of_remaining_jobs == 0)
		{
			last;
		}
		sleep(20);
	}



#	while (1)
#	{
#		$num_of_remaining_jobs = 0;
#		open(PIPE,'/usr/local/bin/qstat|') or die "qstat failed, make sure use the full path of qstat\n";   
#		my @bjobs_list=<PIPE>;           
#		close PIPE;  
#		for (my $i=1; $i<=$#bjobs_list; $i++)
#		{
#			if ($bjobs_list[$i] =~ /^(\d+)\./)
#			{
#				my $jobid = $1;
#
#				#$this->log_it("qstat $jobid");
#
#				if (exists($running_job_ids{$jobid}))
#				{
#					$num_of_remaining_jobs++;
#				}
#			}
#		}
#		if ($num_of_remaining_jobs == 0)
#		{
#			last;
#		}
#		sleep(10);
#	}



	my $end_time = localtime;
	($year, $mon, $day, $hour, $min, $sec) = ($end_time->year + 1900, $end_time->mon+1, $end_time->mday, $end_time->hour, $end_time->min, $end_time->sec);
	my $end_time_sec = timelocal_nocheck($sec, $min, $hour, $day, $mon-1, $year);
	my $difference = $end_time_sec - $start_time_sec;
	my $seconds    =  $difference % 60;
	$difference = ($difference - $seconds) / 60;
	my $minutes    =  $difference % 60;
	$difference = ($difference - $minutes) / 60;
	my $hours      =  $difference % 24;
	#printf("%04d/%02d/%02d %02d:%02d:%02d", $year, $mon, $day, $hour, $min, $sec);
	$this->log_it ("FINISH after running $hours hours $minutes minutes $seconds seconds");

}

sub run_cmd
{
	my $this = shift;
	my $cmd = shift;
	my $start_time = localtime();
	my ($year, $mon, $day, $hour, $min, $sec) = ($start_time->year + 1900, $start_time->mon+1, $start_time->mday, $start_time->hour, $start_time->min, $start_time->sec);
	my $start_time_sec = timelocal_nocheck($sec, $min, $hour, $day, $mon-1, $year);
	printf("%04d/%02d/%02d %02d:%02d:%02d", $year, $mon, $day, $hour, $min, $sec);

	if ($this->{dry_run})
	{
		$this->log_it($cmd);
		return;
	}




	print("\tSTART $cmd\n");
	$this->log_it("START $cmd");

	my $cmd_res = system($cmd);
	if ($cmd_res) 
	{

		my $error_time = localtime;
		($year, $mon, $day, $hour, $min, $sec) = ($error_time->year + 1900, $error_time->mon+1, $error_time->mday, $error_time->hour, $error_time->min, $error_time->sec);
		printf("%04d/%02d/%02d %02d:%02d:%02d", $year, $mon, $day, $hour, $min, $sec);
		print "\tERROR: $?\nEXIT!!\n";
		$this->log_it("ERROR: $?\nEXIT!!");
		exit 1;
	}
	else
	{
		my $end_time = localtime;
		($year, $mon, $day, $hour, $min, $sec) = ($end_time->year + 1900, $end_time->mon+1, $end_time->mday, $end_time->hour, $end_time->min, $end_time->sec);
		my $end_time_sec = timelocal_nocheck($sec, $min, $hour, $day, $mon-1, $year);
		my $difference = $end_time_sec - $start_time_sec;
		my $seconds    =  $difference % 60;
		$difference = ($difference - $seconds) / 60;
		my $minutes    =  $difference % 60;
		$difference = ($difference - $minutes) / 60;
		my $hours      =  $difference % 24;
		printf("%04d/%02d/%02d %02d:%02d:%02d", $year, $mon, $day, $hour, $min, $sec);
		print "\tSUCCESS after running $hours hours $minutes minutes $seconds seconds\n";
		$this->log_it ("SUCCESS after running $hours hours $minutes minutes $seconds seconds");
	}
}

sub log_it
{
	my $this = shift;
	my $log_msg = shift;

	my $current_time = localtime();
	my ($year, $mon, $day, $hour, $min, $sec) = ($current_time->year + 1900, $current_time->mon+1, $current_time->mday, $current_time->hour, $current_time->min, $current_time->sec);
	my $formatted_time = sprintf("%04d/%02d/%02d %02d:%02d:%02d", $year, $mon, $day, $hour, $min, $sec);
	FileSystem::AppendFile($this->{log_file}, "$formatted_time\t$log_msg\n");
	FileSystem::Close();
}

#!/usr/bin/perl
sub get_header
{
	my $r =q(<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Strict//EN">
<html>
<head><title>DNA-Seq pipeline report</title>

<style type="text/css">

 @media screen {
  div.summary {
    width: 18em;
    position:fixed;
    top: 3em;
    margin:1em 0 0 1em;
  }
  
  div.main {
    display:block;
    position:absolute;
    overflow:auto;
    height:auto;
    width:auto;
    top:4.5em;
    bottom:2.3em;
    left:0em;
    right:0;
    border-left: 1px solid #CCC;
    padding:0 0 0 1em;
    background-color: white;
    z-index:1;
  }
  
  div.header {
    background-color: #EEE;
    border:0;
    margin:0;
    padding: 0.5em;
    font-size: 200%;
    font-weight: bold;
    position:fixed;
    width:100%;
    top:0;
    left:0;
    z-index:2;

  }

  div.footer {
    background-color: #EEE;
    border:0;
    margin:0;
	padding:0.5em;
    height: 1.3em;
	overflow:hidden;
    font-size: 100%;
    font-weight: bold;
    position:fixed;
    bottom:0;
    width:100%;
    z-index:2;
  }
  
  img.indented {
    margin-left: 3em;
  }
 }
 
 @media print {
	img {
		max-width:100% !important;
		page-break-inside: avoid;
	}
	h2, h3 {
		page-break-after: avoid;
	}
	div.header {
      background-color: #FFF;
    }
	
 }
 
 body {    
  font-family: sans-serif;   
  color: #000;   
  background-color: #FFF;
  border: 0;
  margin: 0;
  padding: 0;
  }
  
  div.header {
  border:0;
  margin:0;
  padding: 0.5em;
  font-size: 200%;
  font-weight: bold;
  width:100%;
  }    
  
  #header_title {
  display:inline-block;
  float:left;
  clear:left;
  }
  #header_filename {
  display:inline-block;
  float:right;
  clear:right;
  font-size: 50%;
  margin-right:2em;
  text-align: right;
  }

  div.header h3 {
  font-size: 50%;
  margin-bottom: 0;
  }
  
  div.summary ul {
  padding-left:0;
  list-style-type:none;
  }
  
  div.summary ul li img {
  margin-bottom:-0.5em;
  margin-top:0.5em;
  }
	  
  div.main {
  background-color: white;
  }
      
  div.module {
  padding-bottom:1.5em;
  padding-top:1.5em;
  }
	  
  div.footer {
  background-color: #EEE;
  border:0;
  margin:0;
  padding: 0.5em;
  font-size: 100%;
  font-weight: bold;
  width:100%;
  }


  a {
  color: #000080;
  }

  a:hover {
  color: #800000;
  }
      
  h2 {
  color: #800000;
  padding-bottom: 0;
  margin-bottom: 0;
  clear:left;
  }

  table { 
  margin-left: 3em;
  text-align: center;
  }
  
  th { 
  text-align: center;
  background-color: #000080;
  color: #FFF;
  padding: 0.4em;
  }      
  
  td { 
  font-family: monospace; 
  text-align: left;
  background-color: #EEE;
  color: #000;
  padding: 0.4em;
  }

  img {
  padding-top: 0;
  margin-top: 0;
  border-top: 0;
  }

  
  p {
  padding-top: 0;
  margin-top: 0;
  }
  
</style>

</head>



<body>
<div class="header">
<div id="header_title">DNA-Seq Pipeline Report</div>
</div>
);
	return $r;
}

1;


#
# Reference:
# http://www.bbmriwiki.nl/wiki/PipelineCommands
# http://www.bbmriwiki.nl/wiki/SnpCallingPipeline/VariantCalling
# http://www.broadinstitute.org/gatk/gatkdocs/





########################################################################################################
# 0. check fastq format
#sub check_fastq_format
#{
#	my $this = shift;  # Figure out who I am
#	#open(FASTQ, $this->{input_fq1}) or die "can't open $this->{input_fq1}: $@\n";
#	
#
#	print "we are now checking fast_format\n\n";
#	print $this->{sample_name}, "\n";
#	print $this->{java_mem};
#	my $fastq_format;
#   my $counter = 2000; # count 2000 lines of quality scores
#	my $min_qualscore = 1e6;
#	my $max_qualscore = -1e6;
#	my $tmp;
#	while ($counter > 0) 
#	{
#		$tmp = <FASTQ>; # @read name
#		$tmp = <FASTQ>; # base calls
#		$tmp = <FASTQ>; # +read name
#		my $scores = <FASTQ>; # quality scores
#		if (!$scores) { last }
#		#print $scores;
#		chomp($scores);
#		for my $chr (map(ord, split(//, $scores))) 
#		{
#			if ($chr < $min_qualscore) 
#			{
#				$min_qualscore = $chr;
#			}
#			if ($chr > $max_qualscore) 
#			{
#				$max_qualscore = $chr;
#			}
#
#		}
#		$counter--;
#	}
#	# Phred+33 means quality score + 33 = ASCII
#	if ($min_qualscore >= 33 && $max_qualscore <= 126)
#	{
#		$fastq_format = "Sanger";
#	}
#
#	# Solexa+64 means quality score + 64 = ASCII
#	if ($min_qualscore >= 59 && $max_qualscore <= 126)
#	{
#		$fastq_format = "Solexa";
#	}
#
#	# Phred+64 means quality score + 64 = ASCII
#	if ($min_qualscore >= 64 && $max_qualscore <= 126)
#	{
#		$fastq_format = "Illumina 1.3+";
#	}
#
#	# Phred+64 in both both Illumina 1.3+ and Illumina 1.5+
#	#if ($min_qualscore >= 66 && $max_qualscore <= 126)
#	#{
#	#	$fastq_format = "Illumina 1.5+";
#	#}
#
#	# Illumina 1.8+ returned to the use of the Sanger format (Phred+33)
#	print "$fastq_format fastq format: ASCII(min, max) = ($min_qualscore, $max_qualscore)\n";
#	$this->{fastq_format} = $fastq_format;
#	if ($this->{fastq_format} eq "Sanger")
#	{
#		$this->{sanger_fq1} = $this->{input_fq1};
#		$this->{sanger_fq2} = $this->{input_fq2};
#	}
#	return $fastq_format;
#}

########################################################################################################
# 1. QCGroomer step: illumina2sanger
# Alternative converter "fastq_quality_converter -i $this->{input_fq1} -o $this->{sanger_fq1}"
#sub illumina2sanger
#{
#	my $this = shift;  # Figure out who I am
#	my $fastq_converter;
#	if ($this->{fastq_format} eq "Sanger")
#	{
#		print "The input fastq is already in Sanger format, no need for conversion!!\n";
#		return;
#	}
#	elsif ($this->{fastq_format} eq "Solexa")
#	{
#		$fastq_converter = "maq sol2sanger";
#	}
#	elsif ($this->{fastq_format} eq "Illumina 1.3+")
#	{
#		$fastq_converter = "maq ill2sanger";
#	}
#	else
#	{
#		print "No Sanger or Illumina fastq format found, exit now!!\n";
#		exit(0);
#	}
#	$this->run_cmd("$fastq_converter $this->{input_fq1} $this->{sanger_fq1}");
#
#	if($this->{pe}) 
#	{
#		$this->run_cmd("$fastq_converter $this->{input_fq2} $this->{sanger_fq2}");
#	}
#
#}
#
#





# Change hash reference for key with supplied value, or assign default
#	$this->{sample_name}			= $parms{sample_name}	|| 'undef_sample_name';
#	$this->{input_fq1}				= $parms{input_fq1}		|| 0;
#	$this->{input_fq2}				= $parms{input_fq2}		|| 0;
#
#	$this->{reference_fa}			= $parms{reference_fa}	|| 'hg18.fa';
#	$this->{pe}						= $parms{pe}			|| 1;
#	$this->{fastq_format}			= $parms{fastq_format}	|| 'Illumina';
#	$this->{sanger_fq1}				= $parms{sample_name}.'_read1_sanger.fq';
#	$this->{sanger_fq2}				= $parms{sample_name}.'_read2_sanger.fq';
#	$this->{dbsnp_rod}				= $parms{dbsnp_rod}		|| '/data/nextgen1/pipeline/dbsnp_130_hg18.rod';
#	$this->{target_bed}				= $parms{target_bed}	|| '/data/nextgen1/pipeline/targets/SureSelect_All_Exon_G3362_with_names.v2.bed';
#
#	$this->{java_mem}				= $parms{java_mem}		|| '-Xmx4g';
#	$this->{cpu_threads}				= $parms{cpu_threads}	|| 8;
#	$this->{gatk_dir}				= $parms{gatk_dir}		|| '/data/public/tools/gatk-v4333/Sting';
#	$this->{picard_dir}				= $parms{picard_dir}	|| '/home/gmslz/tools/picard-tools-1.44';
#	$this->{tools_dir}				= $parms{tools_dir}		|| '/data/public/tools/gatk_pipeline/src-pipeline';
#	$this->{tmp_dir}				= $parms{tmp_dir}		|| '/home/tmp';

#	print $this->{output_project_dir}."-- base dir***\n";
#	print $this->{gatk_dir}."-- gatk dir***\n";
#	print $this->{tmp_dir}."-- tmp_dir dir***\n";

