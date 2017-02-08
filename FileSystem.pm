package FileSystem;
################################################################
# FileSystem Package                                           #
# Author:   Zhengdeng Lei 09/20/2004                           #
# E-mail:   zlei2 at uic.edu                                   #
# Homepage: http://array.bioengr.uic.edu/~zlei2/               #
################################################################ 
use strict;
use warnings;

my %Filename2Handle=();

################################################################
# Function: Read the whole file to an array                    #
# Input:    File path                                          #
# Output:   An array of file content                           #
# Example:  Input -- "./test.txt"                              #
#			Output -- an array with the content of ./test.txt  #
################################################################
sub GetCmdRes
{
	my $cmd = shift;
	open (CMD,"$cmd |") || die("$!: $cmd failed! Make sure use full path if you call a command!\n");
	my $res = <CMD>;
	close (CMD);
	chop($res);
	return $res;
}

sub NoSlashPath
{
	my $path = shift;
	if ($path =~ /\/$/)
	{
		chop($path)
	}
	return $path;
}

sub GetFileNameFromPath
{
	my ($path) = @_;
	if ($path =~/(.*)\/(.*)/ || $path =~/(.*)\\(.*)/)
	{
		return $2;
	}


}

#sub GetFileNameFromPath
#{
#	my ($path) = @_;
#	my $reverse_path = scalar reverse($path);
#	print "\n****** $reverse_path  ****\n\n";
#	my $filename = "";
#	if ($reverse_path =~/(.*)\/\w/ || $reverse_path =~/(.*)\\\w/)
#	{
#		$filename = $1;
#		print "\n****** $filename  ****\n\n";
#	}
#	return scalar reverse($filename);
#}



sub ReadFile
{
    my ($filename) = @_;
    my @filedata = ();
    open(READ_FILE_DATA, $filename) || die("$!: $filename\nMake sure the file path is correct!!\n");
    @filedata = <READ_FILE_DATA>;
    close READ_FILE_DATA;

    return @filedata;
}

##############################################################
# NOTE: Remember to call FileSystem::CloseFileHandles() when 
#       You don't need the file handles anymore.
# Example:
# FileSystem::WriteFile("E:/a.txt", "a");
# FileSystem::WriteFile("E:/b.txt", "b");
# FileSystem::WriteFile("E:/a.txt", "aa");
# FileSystem::CloseFileHandles();
###############################################################
sub WriteFile
{
	my ($filename, @Buffer) = @_;
	my $fh; ## $fh is a file handle (TYPEGLOB) assigned by the file function: open() 
	if (!exists($Filename2Handle{$filename}))
	{
		###################
		open($fh, ">$filename") || die("$!: $filename\nMake sure the file path is correct!!\n");
		$Filename2Handle{$filename} = $fh;
	}
	else
	{
		$fh = $Filename2Handle{$filename};
	}
	if ($#Buffer>0)
	{
		foreach my $_line (@Buffer)
		{
		   chomp($_line);
		   print $fh $_line, "\n";
		}
	}
	else
	{
		print $fh $Buffer[0];
	}
}

sub AppendFile
{
	my ($filename, @Buffer) = @_;
	my $fh; ## $fh is a file handle (TYPEGLOB) assigned by the file function: open() 
	if (!exists($Filename2Handle{$filename}))
	{
		###################
		open($fh, ">>$filename") || die("$!: $filename\nMake sure the file path is correct!!\n");
		$Filename2Handle{$filename} = $fh;
	}
	else
	{
		$fh = $Filename2Handle{$filename};
	}
	if ($#Buffer>0)
	{
		foreach my $_line (@Buffer)
		{
		   chomp($_line);
		   print $fh $_line, "\n";
		}
	}
	else
	{
		print $fh $Buffer[0];
	}
}


sub CloseFileHandles
{
	if ($#_ >=0) #close one specified file handle
	{
		my $filename = $_[0];
		if (!exists($Filename2Handle{$filename}))
		{
			die("$filename: Fail to close file handle!\n");
		}
		
		close $Filename2Handle{$filename};
		delete $Filename2Handle{$filename};
	}
	else #close all file handles
	{
		foreach my $fh (values %Filename2Handle)
		{

			close $fh;
		}
	   %Filename2Handle = ();
	}
}

sub Close
{
	if ($#_ >=0) #close one specified file handle
	{
		my $filename = $_[0];
		if (!exists($Filename2Handle{$filename}))
		{
			die("$filename: Fail to close file handle!\n");
		}
		
		close $Filename2Handle{$filename};
		delete $Filename2Handle{$filename};
	}
	else #close all file handles
	{
		foreach my $fh (values %Filename2Handle)
		{

			close $fh;
		}
	   %Filename2Handle = ();
	}
}

sub MakeDir
{
	my ($Dir) = @_;
	if (-e $Dir)
	{
		return;
	}
	else
	{
		mkdir($Dir);
	}
}


sub GetSubfolders
{
	my ($Dir, $FullPath) = @_;
	my @subfolders = ();
	print "Reading the folder ...\n";
	opendir(HDIR, $Dir)||die("$!: $Dir\nMake sure it exists!!\n");
	my @FilePaths = readdir(HDIR);
	print "Reading the folder done!!!\n";

	closedir(HDIR);
	for(my $i=2;$i <=$#FilePaths;$i++)
	{
		my $dir_or_file = "";
		if ($FullPath)
		{
			$dir_or_file = $Dir.'\\'.$FilePaths[$i];

		}
		else
		{

			$dir_or_file = $FilePaths[$i];

		}
		if(-d $dir_or_file)
		{
		#	print "$dir_or_file\n";
			push(@subfolders, $dir_or_file);
		}
	}
	return @subfolders;
}


sub GetSubfiles
{
	my ($Dir, $FullPath) = @_;
	my @subfiles = ();
	opendir(HDIR, $Dir)||die("$!: $Dir\nMake sure it exists!!\n");
	my @FilePaths = readdir(HDIR);
	closedir(HDIR);
	for(my $i=2;$i < @FilePaths;$i++)
	{
		my $dir_or_file = "";
		if ($FullPath)
		{
			$dir_or_file = $Dir."\\".$FilePaths[$i];

		}
		else
		{

			$dir_or_file = $FilePaths[$i];
			#print "$dir_or_file\n"

		}
		if(-f $dir_or_file||-f $Dir."\\".$dir_or_file) 
		{
			push(@subfiles, $dir_or_file);
		}
	}
	return @subfiles;
}



#my $src_dir = 'J:\HTS Members\TTakagi\Projects\Ewing\'s sarcoma\ESAR_DRP_23cp_08182008';
#my @files = ();
#@files = GetFileByPattern($src_dir, 'TXT', @files);
sub GetFileByPattern
{
	my ($dir, $pattern, @matchedFiles)=@_;
	$pattern = lc($pattern);
	my $DIR_HANDLE;
	opendir ($DIR_HANDLE, $dir) || die print 'error';
	my @array=readdir($DIR_HANDLE);
	for(my $i=2;$i < @array;$i++)
	{
		my $dir_or_file = $dir.'\\'.$array[$i];
		if (lc($array[$i]) =~ /$pattern/)
		#if((-f $dir_or_file))
		{
			push(@matchedFiles, $dir_or_file);
		}
	}

	########################
	# Go through the sub-dirs
	for(my $i=2;$i < @array;$i++)
	{
		my $dir_or_file = $dir.'\\'.$array[$i];
		if(-d $dir.'\\'.$array[$i])
		{
			@matchedFiles = GetFileByPattern($dir.'\\'.$array[$i], $pattern, @matchedFiles);
		}
	}
    closedir($DIR_HANDLE);
	return @matchedFiles;
}


##############################################################
# Input: param1 -- SD file
#        param2 -- number of entry for each file
# Example:
# FileSystem::SplitSDF("c:\\Progenics.sdf", 50000);
###############################################################
sub SplitSDF
{
	my $sdfile = shift;
	my $num_entry_each_file = shift;

	open(READ_SDF, $sdfile) || die("$!: Make sure the file path is correct!!\n");
	my $prefix = "";
	if ($sdfile =~ /(.*)\.sdf/)
	{
		$prefix = $1;
	}
	my $entry_num = 0;
	my $file_num = 1;
	my $entry_start = ($file_num - 1)*$num_entry_each_file + 1;
	my $entry_end = $file_num*$num_entry_each_file;
	open(WRITE_SDF, ">$prefix.$file_num.$entry_start-$entry_end.sdf") || die("$!: Make sure the file path is correct!!\n");
	while (my $line = <READ_SDF>)
	{
		if ($line =~ /(\s+[\d.-]+\s+[\d.-]+\s+[\d.-]+\s+)(R)(\s+\d.*)/)
		{
			print $line;
			chomp($line);
			$line = $1."H".$3."\n";
		}


		print WRITE_SDF $line;
		if ($line =~ /\$\$\$\$/)
		{
			$entry_num++;
			my $mod_num = $entry_num % $num_entry_each_file;
			if ($mod_num == 0)
			{
				close WRITE_SDF;
				print "$entry_num compounds processed\n";
				$file_num++;
				$entry_start = ($file_num - 1)*$num_entry_each_file + 1;
				$entry_end = $file_num*$num_entry_each_file;
				open(WRITE_SDF, ">$prefix.$file_num.$entry_start-$entry_end.sdf") || die("$!: Make sure the file path is correct!!\n");
			}
		}
	}
	close READ_SDF;
	close WRITE_SDF;
	print "$entry_num compounds processed\n";
	rename("$prefix.$file_num.$entry_start-$entry_end.sdf","$prefix.$file_num.$entry_start-$entry_num.sdf") || die "Can't rename fred to barney: $!";
}

sub SplitSDFaddField
{
}

1;
