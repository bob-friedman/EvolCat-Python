#run perl dsdn_dist.pl <input file> <R> (default is .5, enter ratio or type R to estimate at 3rd pos)
#print "gene1\tgene2\tratio\tps\tSEps\tpn\tSEpn\tds\tSEds\tdn\tSEdn\n";

#codon_number to aa
@aa_array = qw/F F L L L L L L I I I M V V V V S S S S P P P P T T T T A A A A Y Y Z Z H H Q Q N N K K D D E E C C Z W R R R R S S R R G G G G/;


$file = shift @ARGV;
$ratio = shift @ARGV;
if($ratio eq ""){$ratio=.5}

open(FILE,"$file") || die "Cant open file: ($file)";

while($line=<FILE>){
	chomp $line;
	if ($line =~ /^\s*$/){
		next;
	}elsif($line =~ /^\s*#/){
		next;
	}elsif($line =~ /^>/){
		if($switch==1){
			$seq=~s/\s//g;
			push(@seqs,$seq);
			#if ($seq !~ /[^-ACGT]/){die "You should have only the characters -,A,C,G,T in $line, try again.\n";}	
			$seq='';
		}
		$switch=1;
		$line =~ s/>//g;
		push(@names,$line);
	 	next;
	}else{$seq.=$line}
}
push(@seqs,$seq);
close FILE;

#print "$file\n";
$longest=length($seqs[0]);
$num_seqs = scalar(@names);



for ($na=0; $na < $num_seqs; ++$na) {
	for ($n=$na+1; $n < $num_seqs; ++$n) {
	$P1=$Q=$P2=$syn_codons=$nonsyn_codons=$SA_Nei=$SB_Nei=0;
	
	# compute R on 3rd positions
	if($ratio_argv eq 'R'){
		@aln1 = split (//,$seqs[$na]);
		@aln2 = split (//,$seqs[$n]);
		$len=$longest/3;
		for ($x=0;$x<$longest;$x+=3){
			$f=$aln1[$x+2];$s=$aln2[$x+2];
			if($f ne $s){
				if(($f=~/[AG]/)and($s=~/[AG]/)){$P1++}
				elsif(($f=~/[AG]/)and($s=~/[CT]/)){$Q++}
				elsif(($f=~/[CT]/)and($s=~/[AG]/)){$Q++}
				elsif(($f=~/[CT]/)and($s=~/[CT]/)){$P2++}
			}
		}
		
		$P=$P1+$P2;
		if($len==0){print "length 0\n";exit}
		$P=$P/$len;
		$Q=$Q/$len;
		$D=$P+$Q;
		$w1=1/(1-2*$P-$Q);
		$w2=1/(1-(2*$Q));
	
		if(($w1<0)or($w2<0)){
			print "$names[$na]\t$names[$n]\tundef\n";
			next;
		}else{
			$s=.5*log($w1)-.25*log($w2);
			$v=.5*log($w2);
			if($v!=0){
				$ratio=$s/$v;
			}else{
				print "$names[$na]\t$names[$n]\tundef\n";
				next;
			}
		}
	}
	

	print "$names[$na]\t$names[$n]\t";
	
	$ratio_rounded=sprintf("%1.2f",$ratio); 
	print "$ratio_rounded\t";
	
	$x=$v=$index=0;
	$R=$ratio_rounded*2;
	$T=$R/($R+2);

	$Astring = \$seqs[$na];  
	@Aarray = split(//,"$$Astring");
	$A = \@Aarray;

	$Bstring = \$seqs[$n];
	@Barray = split(//,"$$Bstring"); #splits the sequence string into an array
	$B = \@Barray;

	&countsubstitutions; 

	}
}



sub countsubstitutions{
	$count_codons = 0;
	for ($i=1; $i <= $longest; $i += 3) {
		@codonA = ("$$A[$i-1]", "$$A[$i]", "$$A[$i+1]");
		@codonB = ("$$B[$i-1]", "$$B[$i]", "$$B[$i+1]");
		$codA = \@codonA;
		$codB = \@codonB;

		++$count_codons;
		# potential number of synonymous changes
		$syn_siteA[$i] = &syn_site(@$codA);
		$syn_siteB[$i] = &syn_site(@$codB);

		$SA_Nei += $syn_siteA[$i];
		$SB_Nei += $syn_siteB[$i];
				
		#codon_conversion to number 0-63 
		$codon_numberA = &codon_conversion(@$codA);
		$codon_numberB = &codon_conversion(@$codB);
		
		#consider only codons with changes for subsequent steps
		if ($$codA[0] eq $$codB[0] && $$codA[1] eq $$codB[1] && $$codA[2] eq $$codB[2]){next}

		# syn changes of 1 base	
		elsif (($aa_array[$codon_numberA] eq $aa_array[$codon_numberB]) && (($$codA[0] ne $$codB[0] && 
			$$codA[1] eq $$codB[1] && $$codA[2] eq $$codB[2]) || ($$codA[0] eq $$codB[0] && $$codA[1] ne
		       	$$codB[1] && $$codA[2] eq $$codB[2]) || ($$codA[0] eq $$codB[0] && $$codA[1] eq $$codB[1] &&
		       	$$codA[2] ne $$codB[2]))) {$syn_codons++}

		# nonsyn & conserv changes of 1 base
		elsif(($aa_array[$codon_numberA] ne $aa_array[$codon_numberB]) &&
		       	(($$codA[0] ne $$codB[0] && $$codA[1] eq $$codB[1] && $$codA[2] eq $$codB[2]) ||
		       	($$codA[0] eq $$codB[0] && $$codA[1] ne $$codB[1] && $$codA[2] eq $$codB[2]) ||
		       	($$codA[0] eq $$codB[0] && $$codA[1] eq $$codB[1] && $$codA[2] ne $$codB[2]))) 
			{$nonsyn_codons++}

		# 3 elsifs tally 2 base
		#Two base change, example: AAA -> TTA
		elsif($$codA[0] ne $$codB[0] && $$codA[1] ne $$codB[1] && $$codA[2] eq $$codB[2]){
			$x = $$codA[0];
			$y = $$codA[1];
			$$codA[0] = $$codB[0];
			$codon_numberC = &codon_conversion(@$codA);
			$$codA[0] = $x;
			$$codA[1] = $$codB[1];
			$codon_numberD = &codon_conversion(@$codA);
			$$codA[1] = $y;
			$tmp_syn = $temp =0;

			if($aa_array[$codon_numberC] eq "Z"){$temp++}
			else{
				if($aa_array[$codon_numberA] eq $aa_array[$codon_numberC]){$tmp_syn++}
				if($aa_array[$codon_numberB] eq $aa_array[$codon_numberC]){$tmp_syn++}
			}

			if($aa_array[$codon_numberD] eq "Z"){$temp++}
			else{
				if($aa_array[$codon_numberA] eq $aa_array[$codon_numberD]){$tmp_syn++}
				if($aa_array[$codon_numberB] eq $aa_array[$codon_numberD]){$tmp_syn++}
			}
			
			$tmp_syn = $tmp_syn/(2-$temp);
			$syn_codons += $tmp_syn;
			$nonsyn_codons += (2 - $tmp_syn);
		}
	
		#Two base change, example: AAA -> TAT 
		elsif($$codA[0] ne $$codB[0] && $$codA[1] eq $$codB[1] && $$codA[2] ne $$codB[2]){
			$x = $$codA[0];
			$y = $$codA[2];
			$$codA[0] = $$codB[0];
			$codon_numberC = &codon_conversion(@$codA);
			$$codA[0] = $x;
			$$codA[2] = $$codB[2];
			$codon_numberD = &codon_conversion(@$codA);
			$$codA[2] = $y;
			$tmp_syn = $temp =0;

			if($aa_array[$codon_numberC] eq "Z"){$temp++}
			else{
				if($aa_array[$codon_numberA] eq $aa_array[$codon_numberC]){$tmp_syn++}
				if($aa_array[$codon_numberB] eq $aa_array[$codon_numberC]){$tmp_syn++}
			}

			if($aa_array[$codon_numberD] eq "Z"){$temp++}
			else{
				if($aa_array[$codon_numberA] eq $aa_array[$codon_numberD]){$tmp_syn++}
				if($aa_array[$codon_numberB] eq $aa_array[$codon_numberD]){$tmp_syn++}
			}
			
			$tmp_syn = $tmp_syn/(2-$temp);
			$syn_codons += $tmp_syn;
			$nonsyn_codons += (2 - $tmp_syn);
		}

		#Two base change, example: AAA -> ATT 
      		elsif($$codA[0] eq $$codB[0] && $$codA[1] ne $$codB[1] && $$codA[2] ne $$codB[2]){
			$x = $$codA[1];
			$y = $$codA[2];
			$$codA[1] = $$codB[1];
			$codon_numberC = &codon_conversion(@$codA);
			$$codA[1] = $x;
			$$codA[2] = $$codB[2];
			$codon_numberD = &codon_conversion(@$codA);
			$$codA[2] = $y;
			$tmp_syn = $temp =0;

			if($aa_array[$codon_numberC] eq "Z"){$temp++}
			else{
				if($aa_array[$codon_numberA] eq $aa_array[$codon_numberC]){$tmp_syn++}
				if($aa_array[$codon_numberB] eq $aa_array[$codon_numberC]){$tmp_syn++}
			}

			if($aa_array[$codon_numberD] eq "Z"){$temp++}
			else{
				if($aa_array[$codon_numberA] eq $aa_array[$codon_numberD]){$tmp_syn++}
				if($aa_array[$codon_numberB] eq $aa_array[$codon_numberD]){$tmp_syn++}
			}
			
			$tmp_syn = $tmp_syn/(2-$temp);
			$syn_codons += $tmp_syn;
			$nonsyn_codons += (2 - $tmp_syn);
		}


		#all three bases changed
	        # For example AAA -> TTT
        	elsif($$codA[0] ne $$codB[0] && $$codA[1] ne $$codB[1] && $$codA[2] ne $$codB[2]){
			$x = $$codA[0];
			$y = $$codA[1];
			$z = $$codA[2];
			$$codA[0] = $$codB[0];
			$codon_numberC = &codon_conversion(@$codA);
			$$codA[1] = $$codB[1];
			$codon_numberF = &codon_conversion(@$codA);
			$$codA[1] = $y;	
			$$codA[2] = $$codB[2];
			$codon_numberG = &codon_conversion(@$codA);
			$$codA[0] = $x;
			$codon_numberE = &codon_conversion(@$codA);
			$$codA[1] = $$codB[1];
			$codon_numberH = &codon_conversion(@$codA);
			$$codA[2] = $z;
			$codon_numberD = &codon_conversion(@$codA);
			$temp = $tmp_syn =0;

			if(($aa_array[$codon_numberC] eq "Z") or ($aa_array[$codon_numberF] eq "Z")){$temp++}
			else{
				if($aa_array[$codon_numberA] eq $aa_array[$codon_numberC]){$tmp_syn++;}
				if($aa_array[$codon_numberC] eq $aa_array[$codon_numberF]){$tmp_syn++;}
				if($aa_array[$codon_numberB] eq $aa_array[$codon_numberF]){$tmp_syn++;}
			}

			if(($aa_array[$codon_numberC] eq "Z") or ($aa_array[$codon_numberG] eq "Z")){$temp++}
			else{
				if($aa_array[$codon_numberA] eq $aa_array[$codon_numberC]){$tmp_syn++;}
				if($aa_array[$codon_numberC] eq $aa_array[$codon_numberG]){$tmp_syn++;}
				if($aa_array[$codon_numberB] eq $aa_array[$codon_numberG]){$tmp_syn++;}
			}

			if(($aa_array[$codon_numberD] eq "Z") or ($aa_array[$codon_numberF] eq "Z")){$temp++}
			else{
				if($aa_array[$codon_numberA] eq $aa_array[$codon_numberD]){$tmp_syn++;}
				if($aa_array[$codon_numberD] eq $aa_array[$codon_numberF]){$tmp_syn++;}
				if($aa_array[$codon_numberB] eq $aa_array[$codon_numberF]){$tmp_syn++;}
			}
	
			if(($aa_array[$codon_numberD] eq "Z") or ($aa_array[$codon_numberH] eq "Z")){$temp++}
			else{
				if($aa_array[$codon_numberA] eq $aa_array[$codon_numberD]){$tmp_syn++;}
				if($aa_array[$codon_numberD] eq $aa_array[$codon_numberH]){$tmp_syn++;}
				if($aa_array[$codon_numberB] eq $aa_array[$codon_numberH]){$tmp_syn++;}

			}
			if(($aa_array[$codon_numberE] eq "Z") or ($aa_array[$codon_numberG] eq "Z")){$temp++}
			else{
				if($aa_array[$codon_numberA] eq $aa_array[$codon_numberE]){$tmp_syn++;}
				if($aa_array[$codon_numberE] eq $aa_array[$codon_numberG]){$tmp_syn++;}
				if($aa_array[$codon_numberB] eq $aa_array[$codon_numberG]){$tmp_syn++;}
					
			}

			if(($aa_array[$codon_numberE] eq "Z") or ($aa_array[$codon_numberH] eq "Z")){$temp++}
			else{
				if($aa_array[$codon_numberA] eq $aa_array[$codon_numberE]){$tmp_syn++;}
				if($aa_array[$codon_numberE] eq $aa_array[$codon_numberH]){$tmp_syn++;}
				if($aa_array[$codon_numberB] eq $aa_array[$codon_numberH]){$tmp_syn++;}
					
			}
			
			$syn_codons += $tmp_syn/(6-$temp);
			$nonsyn_codons += 3.0 - ($tmp_syn/(6-$temp));
			
		}else{
			#print "Error\t$$codA[0]$$codA[1]$$codA[2]\t$aa_array[&codon_conversion(@$codA)]\t$$codB[0]$$codB[1]$$codB[2]\t$aa_array[&codon_conversion(@$codB)]\n";
		}
	}

	$potential_syn = ($SA_Nei/3 + $SB_Nei/3)/2;
	$potential_nonsyn = (3*$count_codons - $potential_syn);
	
	if($potential_syn == 0 || $potential_nonsyn == 0){print "denom 0\n";exit}
	$ps = $syn_codons/$potential_syn;
	$pn = $nonsyn_codons/$potential_nonsyn;
	print "$syn_codons\t$nonsyn_codons\t$potential_syn\t$potential_nonsyn\n";
}



#syn-site is a function, array is in the same order as the other amino acid translation arrays in this set.
sub syn_site {
	local(@codon)=@_;
	
	$codon_syn_sites[0]=$T*3;
	$codon_syn_sites[1]=$T*3;
	$codon_syn_sites[2]=(2*$T)*3;
	$codon_syn_sites[3]=(2*$T)*3;
	$codon_syn_sites[4]=3;
	$codon_syn_sites[5]=3;
	$codon_syn_sites[6]=(1+$T)*3;
	$codon_syn_sites[7]=(1+$T)*3;
	$codon_syn_sites[8]=(($R+1)/($R+2))*3;
	$codon_syn_sites[9]=(($R+1)/($R+2))*3;
	$codon_syn_sites[10]=(2/($R+2))*3;
	$codon_syn_sites[11]=0;
	$codon_syn_sites[12]=3;
	$codon_syn_sites[13]=3;
	$codon_syn_sites[14]=3;
	$codon_syn_sites[15]=3;
	$codon_syn_sites[16]=3;
	$codon_syn_sites[17]=3;
	$codon_syn_sites[18]=3;
	$codon_syn_sites[19]=3;
	$codon_syn_sites[20]=3;
	$codon_syn_sites[21]=3;
	$codon_syn_sites[22]=3;
	$codon_syn_sites[23]=3;
	$codon_syn_sites[24]=3;
	$codon_syn_sites[25]=3;
	$codon_syn_sites[26]=3;
	$codon_syn_sites[27]=3;
	$codon_syn_sites[28]=3;
	$codon_syn_sites[29]=3;
	$codon_syn_sites[30]=3;
	$codon_syn_sites[31]=3;
	$codon_syn_sites[32]=3;
	$codon_syn_sites[33]=3;
	$codon_syn_sites[34]=0;
	$codon_syn_sites[35]=0;
	$codon_syn_sites[36]=$T*3;
	$codon_syn_sites[37]=$T*3;
	$codon_syn_sites[38]=$T*3;
	$codon_syn_sites[39]=$T*3;
	$codon_syn_sites[40]=$T*3;
	$codon_syn_sites[41]=$T*3;
	$codon_syn_sites[42]=$T*3;
	$codon_syn_sites[43]=$T*3;
	$codon_syn_sites[44]=$T*3;
	$codon_syn_sites[45]=$T*3;
	$codon_syn_sites[46]=$T*3;
	$codon_syn_sites[47]=$T*3;
	$codon_syn_sites[48]=($R/($R+1))*3;
	$codon_syn_sites[49]=($R/($R+1))*3;
	$codon_syn_sites[50]=0;
	$codon_syn_sites[51]=0;
	$codon_syn_sites[52]=3;
	$codon_syn_sites[53]=3;
	$codon_syn_sites[54]=1.5*3;
	$codon_syn_sites[55]=(1+1/($R+2))*3;
	$codon_syn_sites[56]=$T*3;
	$codon_syn_sites[57]=$T*3;
	$codon_syn_sites[58]=($T+1/($R+1))*3;
	$codon_syn_sites[59]=($T+1/($R+2))*3;
	$codon_syn_sites[60]=3;
	$codon_syn_sites[61]=3;
	$codon_syn_sites[62]=3;
	$codon_syn_sites[63]=3;

	$codon_number = &codon_conversion(@codon);
	return $codon_syn_sites[$codon_number];
}

sub codon_conversion {
	local(@codon)=@_; 
	#hash base to number
	%baseNumber = ("T" => 0, "C" => 1, "A" => 2, "G" => 3);

	my ($xleft,$xmid,$xright);
	$xleft = $baseNumber{$codon[0]};
	$xmid = $baseNumber{$codon[1]};
	$xright = $baseNumber{$codon[2]};
	return ($xmid*16) + ($xleft*4) + $xright;
}
