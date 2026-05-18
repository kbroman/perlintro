#!/usr/local/bin/perl

# create data for a perl example
#
# simulate five nuclear families with no. sibs = (3, 2, 8, 4, 5) 
# at the markers in a ~/Projects/WGP/Data/wgp20.frq

srand();

@nsib = (3, 2, 8, 4, 5);

# read allele frequencies
$ifile = "/home1/biostats/kbroman/Projects/WGP/Data/wgp20.frq";
open(IN, $ifile) or die("Cannot read from $ifile");
$nmar = <IN>;
foreach $i (0..($nmar-1)) {
    $line = <IN>; chomp($line);
    ($mar, $nalle) = split(/\s+/, $line);
    $mar[$i] = $mar;
    $wh{$mar} = $i;
    $nalle{$mar} = $nalle;
    foreach $j (1..$nalle) {
	$line = <IN>; chomp($line);
	@v = split(/\s+/, $line);
	$freq{$mar[$i]}{$v[0]} = $v[1];
    }
}
close(IN);

# print list of markers
$ofile = "markers.txt";
open(OUT, ">$ofile") or die("Cannot write to $ofile");
foreach $i (0..($nmar-1)) {
    print OUT ("$mar[$i]\n");
}
close(OUT);

# print family information
$ofile = "families.txt";
open(OUT, ">$ofile") or die("Cannot write to $ofile");
print OUT ("Family Individual Father Mother Sex\n");
foreach $i (0..(@nsib-1)) {
    printf OUT ("%6d %10d %6d %6d %3d\n", $i+1, 1, 0, 0, 1);
    printf OUT ("%6d %10d %6d %6d %3d\n", $i+1, 2, 0, 0, 2);
    foreach $j (1..$nsib[$i]) {
	printf OUT ("%6d %10d %6d %6d %3d\n", $i+1, $j+2, 1, 2, rand()+1.5);
    }
}
close(OUT);

# simulate genotype data
foreach $i (0..($nmar-1)) {

    foreach $j (1..@nsib) {
	
	# simulate dad
	$gen{$j}{1}[$i*2] = sim_allele($mar[$i]);
	$gen{$j}{1}[$i*2+1] = sim_allele($mar[$i]);
	
	# simulate mom
	$gen{$j}{2}[$i*2] = sim_allele($mar[$i]);
	$gen{$j}{2}[$i*2+1] = sim_allele($mar[$i]);

	# simulate kids 
	foreach $k (1..$nsib[$j-1]) {
	    if(rand() < 0.5) {
		$gen{$j}{2+$k}[$i*2] = $gen{$j}{1}[$i*2];
	    }
	    else {
		$gen{$j}{2+$k}[$i*2] = $gen{$j}{1}[$i*2+1];
	    }
	    if(rand() < 0.5) {
		$gen{$j}{2+$k}[$i*2+1] = $gen{$j}{2}[$i*2];
	    }
	    else {
		$gen{$j}{2+$k}[$i*2+1] = $gen{$j}{2}[$i*2+1];
	    }
	}

	foreach $k (1..($nsib[$j-1]+2)) {
	    if(rand() < 0.1) {
		$gen{$j}{$k}[$i*2] = $gen{$j}{$k}[$i*2+1] = 0;
	    }
	}

    }
}

# write genotypes
$ofile = "genotypes.txt";
open(OUT, ">$ofile") or die("Cannot write to $ofile");

# header line
printf OUT ("%-9s", "Marker");
foreach $j (1..@nsib) {
    foreach $k (1..($nsib[$j-1]+2)) {
	printf OUT ("%-7s", $j . "-" . $k);
    }
}
print OUT ("\n");

foreach $mar (sort @mar) {
    printf OUT ("%-9s", $mar);

    $i = $wh{$mar};

    foreach $j (1..@nsib) {
	foreach $k (1..($nsib[$j-1]+2)) {
	    if($gen{$j}{$k}[$i*2] == 0 or $gen{$j}{$k}[$i*2+1]==0) {
		printf OUT ("%-7s", "");
	    }
	    else {
		printf OUT ("%-7s", $gen{$j}{$k}[$i*2] . "/" . $gen{$j}{$k}[$i*2+1]);
	    }
	}
    }

    print OUT ("\n");
}
close(OUT);


sub sim_allele 
{
    ($mar) = @_;
    
    $r = rand();
    
    foreach $alle (sort numerically keys %{$freq{$mar}}) {
	if($r < $freq{$mar}{$alle}) {
	    return($alle);
	}
	else {
	    $r -= $freq{$mar}{$alle};
	}
    }
    # this shouldn't happen
    return(-1);
}

sub numerically { $a <=> $b; }

	    
    
