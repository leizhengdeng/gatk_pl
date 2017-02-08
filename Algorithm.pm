package Algorithm;
################################################################
# Algorithm Library                                            #
# Author:	Zhengdeng Lei 09/20/2004                           #
# E-mail:   zlei2 at uic.edu                                   #
# Homepage: http://array.bioengr.uic.edu/~zlei2/               #
################################################################ 
use strict;
use warnings;

my $EPSLON = 1E-10;


sub FormatNum
{
	my $num  = $_[0];
	if ($num < 10)
	{
		$num = '0'.$num;
	}
	return $num;
}

sub FormatWell
{
	my $well_name  = $_[0];
	if ($well_name =~ /([A-Z]+)(\d+)/)
	{
		$well_name = sprintf("%s%02d", $1, $2);
	}
	return $well_name;
}

sub Well_letter2num
{
	my $Well_letter = $_[0];

	my @letters = split(//, $Well_letter);
	my $num = 0;
	for (my $pos=0; $pos<=$#letters; $pos++)
	{
		my $pow = $#letters - $pos;
		$num += ord($letters[$pos]) - 65 + 26**$pow;

	}
	return $num;
}

sub Well_num2letter
{
	my $Well_num = $_[0];
	my $letters = '';
	while ($Well_num>26)
	{
		my $letter = int(($Well_num-1)/26);
		$letters .= chr($letter+64);
		$Well_num = $Well_num - $letter*26;

	}
	$letters .= chr($Well_num+64);
	return $letters;
}

sub trim_space
{
	my $string = shift;
	if ($string ne "")
	{
		$string =~ s/^\s+//;
		$string =~ s/\s+$//;
	}
	return $string;
}



# Left trim function to remove leading whitespace
sub ltrim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	return $string;
}
# Right trim function to remove trailing whitespace
sub rtrim($)
{
	my $string = shift;
	$string =~ s/\s+$//;
	return $string;
}



sub CheckDuplicate
{
	my @request_in = @_;
	my %no_duplicates = ();
	my @no_duplicates = ();

	my %duplicates = ();
	my @duplicates = ();
	for (my $i=0; $i<=$#request_in; $i++)
	{
		$request_in[$i] =~ s/^\s+//;
		$request_in[$i] =~ s/\s+$//;
		if (exists($no_duplicates{$request_in[$i]}))
		{
			$duplicates{$request_in[$i]} = 1;
		}
		else
		{
			$no_duplicates{$request_in[$i]} = 1;
			push(@no_duplicates, $request_in[$i]);
		}
	}
	@duplicates = sort(keys(%duplicates));
	return @duplicates;
}






################################################################
# Function:	partition of the data                              #
# Input:	$K - K-fold cross valiation                        #
#           $DATAT_NUM - number of data in the training set    # 
# Output:	partition of the data for each fold                #
# Example:	Algorithm::CVPartition(5, 32)                      #
#			Output -- (7, 7, 6, 6, 6)                          #
# Author:	Zhengdeng Lei 09/28/2004                           #
################################################################
sub CVPartition
{
	#$K-fold cross validation
	my ($K, $DATAT_NUM) = @_;
	my $ONE_FOLD = int($DATAT_NUM/$K);
	my $MOD = $DATAT_NUM - $ONE_FOLD*$K;
	my @PARTITION = ();
	for (my $i=0; $i<$K; $i++)
	{
		if ($i < $MOD)
		{
			$PARTITION[$i] = $ONE_FOLD + 1;
		}
		else
		{
			$PARTITION[$i] = $ONE_FOLD;
		}
	}
	return @PARTITION;
}



##############################################################################################
#  @votes;
#  for ($cl=0; $cl<=$CLASS_NUM; $cl++)
#  {
#	$votes[0][$cl] = 0; #votes
#	$votes[1][$cl] = 0; #use function margin score when ties arise
#  }
#  for (my $i=0; $i<=$#AA_sequences;$i++)
#  {
#	# Now dealing with one seq.
#	# Clear votes
#	for ($cl=0; $cl<=$CLASS_NUM; $cl++)
#	{
#		$votes[0][$cl] = 0;
#		$votes[1][$cl] = 0;
#	}
#
#	for ($model1=1; $model1<= $CLASS_NUM-1; $model1++)
#	{
#		for ($model2=$model1+1; $model2<= $CLASS_NUM; $model2++)
#		{
#			if ($D1X1_ALLTestRes{"$model1$model2"}[$i] >= 0)
#			{
#				$votes[0][$model1] += 1;                 #vote
#				$votes[1][$model1] += $D1X1_ALLTestRes{"$model1$model2"}[$i];      #function margin
#			}
#			else
#			{
#				$votes[0][$model2] += 1;                 #vote
#				$votes[1][$model2] += abs($D1X1_ALLTestRes{"$model1$model2"}[$i]); #function margin
#			}
#      }
#  }
#  $prediction =  Algorithm::JuryVoting(@votes);
sub JuryVoting
{
	my @VOTES = @_;
	my @index = Algorithm::GetIndexOfMaxVal(@{$VOTES[0]});
	my $PredictedClass = 0;
	if ($#index > 0) # there is a tie, more than one value has max val
	{
		my $max = -9E10;
		for (my $ind=0; $ind<=$#index; $ind++)
		{
			if ($VOTES[1][$index[$ind]] > $max)
			{
				$max = $VOTES[1][$index[$ind]];
				$PredictedClass = $index[$ind];
			}
		}
	}
	else # There is no tie
	{
		$PredictedClass = $index[0];
	}
	return $PredictedClass;
}



################################################################
# Function:	Get the indices of the maximum values for an array #
# Input:	A numerical array                                  #
# Output:	Index array                                        #
# Example:	Input -- (1.1, 2, 3.3, 3.3, 0)                     #
#			Output -- (2, 3)                                   #
# called by Algorithm::JuryVoting                              #
# Author:	Zhengdeng Lei 09/28/2004                           #
################################################################
sub GetIndexOfMaxVal
{
	my $EPSLON = 1E-6;
    my @array = @_;
	my @index = ();
	my $max = -9E10;

	for (my $i=0; $i<=$#array; $i++)
	{
		if (abs($array[$i]-$max) < $EPSLON) #$array[$i] == $max
		{
			$max = $array[$i];
			push(@index, $i);
		}
		elsif ($array[$i] > $max)
		{
			$max = $array[$i];
			@index = ();
			push(@index, $i);
		}
	}
    return @index;
}

################################################################
# Function:	Get the descending sorted indices of an array      #
# Input:	A numerical array                                  #
# Output:	Index array                                        #
# Example:	Input -- @a = (1.1, 2, 3.3, 3.3, 0);               #
#			Output -- (2, 3, 1, 0, 4)                          #
#	        #@b = sort{$a <=> $b}(@a);                         #
#           #@b = (3.3, 3.3, 2, 1.1, 0)                        #
#           #The elements in @b correspond to the indices in @a#
# Author:	Zhengdeng Lei 09/28/2004                           #
################################################################
sub SortIndex
{
    my @array = @_;
	my @index = (0);

	if ($#array < 0)
	{
		print ("Algorithm::SortIndex input array is empty!\n");
		return;
	}

	for (my $i=1; $i<=$#array; $i++)
	{
		my $inserted = 0;
		for (my $j=0; $j<=$#index; $j++)
		{
			if($array[$i] > $array[$index[$j]])
			{
				splice(@index, $j, 0, $i);
				$inserted = 1;
				last;
			}
		}
		if (!$inserted) #
		{
			push(@index, $i);
		}
	}
    return @index;
}



##################################################################################
# Function:	Get the MCC values of prediction                                     # 
# Input: 1. $K -- is the number of class                                         #
#        2. Array of data points	                                             #
#        @DATA_POINTS = ([predicted class, true class],                          #
#                        [predicted class, true class],                          #
#                                     ... ...                                    #
#                        [predicted class, true class]);                         #
#        NOTE: predicted class, true class rangs 1 ... $K                        #
#        $DATA_POINTS[$i][0] is the predicted class for $i_th data point         #
#        $DATA_POINTS[$i][1] is the true class for $i_th data point              #
# Output:	MCC: $MCC[0] is the overall MCC                                      #
#                $MCC[1] is the MCC of 1st class                                 #
#                $MCC[$K] is the MCC of $K_th class                              #
# Example:	@mcc = Algorithm::CalMCC($CLASS_NUM, @datapoints);                   #
# Author:	Zhengdeng Lei 11/01/2006                                             #
##################################################################################
#suppose there are K classes,
#for one specific class k, MCCk = (pk*nk - uk*ok)/sqrt((pk+uk)*(pk+ok)*(nk+uk)*(nk+ok))
#Mst  where s = predicted class, t = true class
#pk = Mkk                         correctly predicted classes belonging in class k
#nk = sum[s!=k]sum[t!=k]Mst       correctly predicted classes not belonging to class k
#uk = sum[s!=k]Msk                under-predicted points
#ok = sum[t!=k]mkt                over-predicted points
sub CalMCC
{
	#$K -- number of classes K
	#@DATA_POINTS = ([predicted class, true class],
	#                [predicted class, true class],
	#                      ... ...
	#                [predicted class, true class]);
	#$DATA_POINTS[$i][0] is the predicted class for $i_th data point
	#$DATA_POINTS[$i][1] is the true class for $i_th data point
	my ($K, @DATA_POINTS) = @_;
	my @MCC = ();
	my $mcc = 0;
	for (my $k=1; $k<=$K; $k++)
	{
		my $pk = 0;
		my $nk = 0;
		my $uk = 0;
		my $ok = 0;
		for (my $i=0; $i<=$#DATA_POINTS; $i++)
		{
			if (($DATA_POINTS[$i][0]==$k) && ($DATA_POINTS[$i][1]==$k))
			{
				$pk++;
			}
			if (($DATA_POINTS[$i][0]!=$k) && ($DATA_POINTS[$i][1]!=$k))
			{
				$nk++;
			}
			if (($DATA_POINTS[$i][0]!=$k) && ($DATA_POINTS[$i][1]==$k))
			{
				$uk++;
			}
			if (($DATA_POINTS[$i][0]==$k) && ($DATA_POINTS[$i][1]!=$k))
			{
				$ok++;
			}
		}
		my $denominator = sqrt(($pk+$uk)*($pk+$ok)*($nk+$uk)*($nk+$ok));
		if ($denominator<1E-10)
		{
			$denominator = 1E-10;
		}
		$MCC[$k] = ($pk*$nk - $uk*$ok)/$denominator;
		$mcc += $MCC[$k];
	}
	$MCC[0] = $mcc/$K;
	return @MCC; # $MCC[0] is the overall mcc
}

##################################################################################
# Function:	Get the Qacc values of prediction                                    # 
# Input: 1. $K -- is the number of class                                         #
#        2. Array of data points	                                             #
#        @DATA_POINTS = ([predicted class, true class],                          #
#                        [predicted class, true class],                          #
#                                     ... ...                                    #
#                        [predicted class, true class]);                         #
#        NOTE: predicted class, true class rangs 1 ... $K                        #
#        $DATA_POINTS[$i][0] is the predicted class for $i_th data point         #
#        $DATA_POINTS[$i][1] is the true class for $i_th data point              #
# Output:	Qacc: $Qacc[0] is the overall accuracy                               #
#                 $Qacc[1] is the accuracy of 1st class                          #
#                 $Qacc[$K] is the accuracy of $K_th class                       #
# Example:	@qacc = Algorithm::CalQacc($CLASS_NUM, @datapoints);                 #
# Author:	Zhengdeng Lei 11/01/2006                                             #
##################################################################################
sub CalQacc
{
	my ($K, @DATA_POINTS) = @_;
	my @Qacc = ();
	my @mk = ();
	my @pk = ();
	my $m = 0;
	my $p = 0;
	for (my $k=1; $k<=$K; $k++)
	{
		$mk[$k] = 0;
		$pk[$k] = 0;
		for (my $i=0; $i<=$#DATA_POINTS; $i++)
		{
			if ($DATA_POINTS[$i][1]==$k) #true class
			{
				$mk[$k]++;
				$m++;
				if ($DATA_POINTS[$i][0] == $k)
				{
					$pk[$k]++;
					$p++;
				}
			}
		}
		if ($mk[$k] !=0)
		{
			$Qacc[$k] = $pk[$k]/$mk[$k];
		}
	}
	if ($m !=0)
	{
		$Qacc[0] = $p/$m;
	}
	return @Qacc; # $Qacc[0] is the overall accuracy
}

##################################################################################
# Function:	Get the evalution values of prediction                               # 
# Input: Array of data points	                                                 #
#        @DATA_POINTS = ([predicted class, true class], postive = 1, negative =2 #
#                        [predicted class, true class],                          #
#                                     ... ...                                    #
#                        [predicted class, true class]);                         #
#        $DATA_POINTS[$i][0] is the predicted class for $i_th data point         #
#        $DATA_POINTS[$i][1] is the true class for $i_th data point              #
# Output:	TP, TN, FP, FN, Accuracy, Precision, Recall, F_score, MCC,           #
#           Specificity, Sensitivity                                             #  
# NOTE: Is it better to output a hash of the values                              #
# Example:	my @test_res = Algorithm::ConfusionMatrix(@datapoints);              #
# FileSystem::WriteFile("$dir/detail.out", "g=$gg c=$cc j=$jj:\tTP=$test_res[0] TN=$test_res[1] FP=$test_res[2] FN=$test_res[3] ");
# FileSystem::WriteFile("$dir/detail.out", "Accuracy=$test_res[4] Precision=$test_res[5] Recall=$test_res[6] F_score=$test_res[7] MCC=$test_res[8] Specificity=$test_res[9] Sensitivity=$test_res[10] ");
# Author:	Zhengdeng Lei 11/01/2006                                             #
##################################################################################
sub ConfusionMatrix
{
	#@DATA_POINTS = ([predicted class, true class], postive = 1, negative =2
	#                [predicted class, true class],
	#                      ... ...
	#                [predicted class, true class]);
	#$DATA_POINTS[$i][0] is the predicted class for $i_th data point
	#$DATA_POINTS[$i][1] is the true class for $i_th data point
	my @DATA_POINTS = @_;
	my $TP=0;
	my $TN=0;
	my $FP=0;
	my $FN=0;

	for (my $i=0; $i<=$#DATA_POINTS; $i++)
	{
		if ($DATA_POINTS[$i][0]==1 && $DATA_POINTS[$i][1]==1) {$TP++;}
		if ($DATA_POINTS[$i][0]==2 && $DATA_POINTS[$i][1]==2) {$TN++;}
		if ($DATA_POINTS[$i][0]==1 && $DATA_POINTS[$i][1]==2) {$FP++;}
		if ($DATA_POINTS[$i][0]==2 && $DATA_POINTS[$i][1]==1) {$FN++;}
	}
	my $Accuracy = ($TP+$TN)/($TP+$TN+$FP+$FN);
	my $Precision = 0;
	my $Recall = 0;
	my $F_score = 0;
	if ($TP > 0)
	{
		$Precision = $TP/($TP+$FP);
		$Recall = $TP/($TP+$FN);
		$F_score = 2*$Precision*$Recall/($Precision + $Recall+$EPSLON);
	}
	my $MCC = 0;
	if (sqrt(($TP+$FN)*($TP+$FP)*($TN+$FP)*($TN+$FN)) > 1E-9)
	{
		$MCC= ($TP*$TN - $FP*$FN)/sqrt(($TP+$FN)*($TP+$FP)*($TN+$FP)*($TN+$FN));
	}
	my $Specificity = $TN/($TN+$FP);
	my $Sensitivity = $Recall;

	my @res = ($TP, $TN, $FP, $FN, $Accuracy, $Precision, $Recall, $F_score, $MCC, $Specificity, $Sensitivity);
	for (my $r=4; $r<=$#res; $r++)
	{
		$res[$r] = sprintf("%0.4f", $res[$r]);
	}
	return @res;
}


##################################################################################
# Function:	Get the AUC of ROC                                                   # 
# Input:	1. the REFERENCE of the array of true positive rate (sensitivity)    #
#           2. the REFERENCE of the array of false positive  (1-specficity)      #
#           3. (optional) the csv filename of ROC TPR vs. FPR for output         #
# Output:	AUC                                                                  #
# Example:	Algorithm::calAUC(\@TPR, \@FPR, "./roc.csv")                         #
#           Algorithm::calAUC(\@TPR, \@FPR)                                      #
# Author:	Zhengdeng Lei 03/26/2007                                             #
##################################################################################
sub calAUC
{
	my @TPR =  @{$_[0]};
	my @FPR = @{$_[1]};
	my $ROC_filename = $_[2];
	if ($#TPR != $#FPR)
	{
		print "Error: The dimesions of @TPR and @FPR are different!\n";
	}
	push (@TPR, (1, 0));
	push (@FPR, (1, 0));
	@TPR = sort {$b <=> $a} @TPR; #sort in descending order
	@FPR = sort {$b <=> $a} @FPR; #sort in descending order
	my $AUC = 0;
	for (my $i=1; $i<=$#TPR; $i++)
	{
		$AUC += ($TPR[$i-1]+$TPR[$i])/2*($FPR[$i-1]-$FPR[$i]);

	}

	if ($ROC_filename)
	{
	    open(ROC, ">$ROC_filename");
		print ROC "1-Specificity(FPR), Sensitivity(TPR), , AUC=$AUC\n";
		for (my $i=0; $i<=$#TPR; $i++)
		{

			print ROC sprintf("%.4f,%.4f\n", $FPR[$i], $TPR[$i]);
		}
		close(ROC);
	}

	return sprintf("%.4f", $AUC);
}	


#########################################################################################
# Function:	Get the AUC of ROC for SVM                                                  # 
# Input:	1. the REFERENCE of the array of SVM sample test (pos=1 neg=2 or -1)        #
#           2. the REFERENCE of the array svm output, i.e. the array of function scores #
#           3. (optional) the csv filename of ROC TPR vs. FPR for output                #
# The input obtained from "./svm_classify test svm_model svm_output >scree.out"         # 
# 					    	if ($test[$t] =~ /^-/){	$SVM_test[$t] = 2;#negative	}       #
#							else {$SVM_test[$t] = 1}                                    #
# Output:	AUC                                                                         #
# Example:	Algorithm::cal_SVM_AUC(\@SVM_test, \@svm_output, "./roc.csv")               #
#           Algorithm::calAUC(\@SVM_test, \@svm_output)                                 #
# Author:	Zhengdeng Lei 03/26/2007                                                    #
#########################################################################################
sub cal_SVM_AUC
{
	my @svm_test = @{$_[0]};
	my @svm_output = @{$_[1]};
	my $ROC_file = $_[2];
	my @svm_output_sorted = sort{$a <=> $b}(@svm_output);
	my @threshold = ();
	push(@threshold, -1E20);
	for (my $s=0; $s<$#svm_output_sorted; $s++)
	{
		push(@threshold, (($svm_output_sorted[$s] + $svm_output_sorted[$s+1])/2));
	}
	push(@threshold, 1E20 );
	my @TPR = ();
	my @FPR = ();
	for (my $t=0; $t<=$#threshold; $t++)
	#for (my $t=0; $t<=0; $t++)
	{
		my $TP=0;
		my $TN=0;
		my $FP=0;
		my $FN=0;
		for (my $s=0; $s<=$#svm_output; $s++)
		{
			if ($svm_output[$s]  > $threshold[$t])
			{
				#predicted =1 print "1";
				if ($svm_test[$s] == 1)	
				{
					$TP++;
				}
				else
				{
					$FP++
				}
			}
			else
			{
				#predicted = 2 print "2";
				if ($svm_test[$s] == 1)	
				{
					$FN++;
				}
				else
				{
					$TN++
				}
			}
		}
		$TPR[$t] = $TP/($TP+$FN); #or Sensitivity
		$FPR[$t] = $FP/($FP+$TN); #or 1-Specificity
	#	FileSystem::WriteFile("./roc.csv", "$FPR[$t],$TPR[$t]\n");
	}
	
	return calAUC(\@TPR, \@FPR, $ROC_file);
}



#########################################################################################
# Function:	Get inner product of two vectors (strings)                                  # 
# Input:	1. vector 1 in string format                                                #
#           2. vector 2 in string format                                                #
# Output:	inner product                                                               #
# Example:	$a = "1:1 9:2 12:1";                                                        #
#           $b = "3:1 9:1 12:1 14:1 60:1";                                              #
#           print inner_product($a, $b);                                                #
# Author:	Zhengdeng Lei 06/20/2007                                                    #
#########################################################################################
sub inner_product
{
	my ($a, $b) = @_;
	my ($a_index, $a_score, $b_index, $b_score);
	print "$a\n$b\n";
	my @a_index_score = split(/\s/, $a);
	my @b_index_score = split(/\s/, $b);
	my $inner = 0;
	my $i=0;
	my $j=0;
	while (1)
	{
		#print "i=$i\tj=$j\n";
	   if ($a_index_score[$i] =~ /([0-9]*):([.0-9]*)/)
		{
			$a_index = $1;
			$a_score = $2;
			#print "a_index = $1\ta_score = $2\n";
		}
		if ($b_index_score[$j] =~ /([0-9]*):([.0-9]*)/)
		{
			$b_index = $1;
			$b_score = $2;
			#print "b_index = $1\tb_score = $2\n";
		}

		if ($a_index==$b_index)
		{
			#print "a_index == b_index\n";
			$inner += $a_score*$b_score;
			$i++;
			$j++;
			if ($i>$#a_index_score && $j>$#b_index_score)
			{
				last;
			}
		}
		elsif ($a_index>$b_index)
		{
			#print "a_index > b_index\n";
			$j++;
			if ($j>$#b_index_score)
			{
				last;
			}
		}
		else
		{
			#print "a_index < b_index\n";
			$i++;
			if ($i>$#a_index_score)
			{
				last;
			}
		}
	}
	return $inner;
} 

sub max
{
	my @v = @_;
	my $max = $v[0];
	foreach my $value (@v)
	{
		if ($max<$value)
		{
			$max = $value;
		}
	}
	return $max;

}

sub min
{
	my @v = @_;
	my $min = $v[0];
	foreach my $value (@v)
	{
		if ($min>$value)
		{
			$min = $value;
		}
	}
	return $min;

}

sub median
{
	my @v = @_;
	my @sorted_v = sort{$a <=> $b}(@v);
	my $max_idx = $#v;
	if (($max_idx%2)==0)
	{
		return $sorted_v[$max_idx/2];
	} 
	else
	{
		return ($sorted_v[($max_idx-1)/2]+$sorted_v[($max_idx+1)/2])/2;

	}

}

sub scale2zero_one
{
	my @v = @_;
	my @scaled_v = ();
	my $min = min(@v);
	my $max = max(@v);
	my $range = $max - $min;
	for (my $i=0; $i<=$#v; $i++)
	{
		$scaled_v[$i] = ($v[$i]-$min)/$range;
	}
	return @scaled_v;
}



sub MakeNotNegative
{
    my ($src, $dest) = @_;

    open(SRC_FILE_DATA, $src) || die("$!: $src\nMake sure the file path is correct!!\n");
	open (DEST_FILE_DATA, ">$dest");
	my $i=0;
	my $num_negative = 0;
    while (my $filedata = <SRC_FILE_DATA>)
	{
		$i++;
		chomp($filedata);
		my @values = split(/\t/, $filedata);


		for (my $v=0; $v<$#values; $v++)
		{
			if ($values[$v] =~ /^-[\d.+-eE]+/ && $values[$v] <=0)
			{
				$values[$v] = 0.000001;
				$num_negative++;
			}
			 print DEST_FILE_DATA "$values[$v]\t";
		}
		my $last_value = $values[$#values];
		if ($last_value =~ /^-[\d.+-eE]+/ && $last_value <=0)
		{
			$last_value = 0.000001;
			$num_negative++;

		}
		 print DEST_FILE_DATA "$last_value\n";
    }
    close SRC_FILE_DATA;
	close DEST_FILE_DATA; 
	print "Number of negatives found: $num_negative\n";
}




1;