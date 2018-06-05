use strict;
package DHspFormat;


# loads a file containing the high scoring pair information into a hash
# input
# target mrna                   start   end     hspID           binding motif id(hash)
# NM_175636|2045..3362|Homo	576	607	HSP9247007108	AAAAAAA
# NM_018191|1858..4042|Homo	1901	1932	HSP9247007304	AAAAAAA
# NM_018364|2446..6604|Homo	4091	4122	HSP9247007400	AAAAAAA
# NM_001429|7648..8765|Homo	47	78	HSP9247007442	AAAAAAA
# NM_014912|2271..5840|Homo	2281	2312	HSP9247007453	AAAAAAA
#
sub loadFile{
    my ($filename)=@_;

    open(FHANDLE,$filename) or die "Cannot open $filename \n";

    my %hspHash;
    while(<FHANDLE>){
	if(/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/){
	    if (defined $hspHash{$5}){
		push(@{$hspHash{$5}},[$1,$2,$3,$4]);

	    }else{
		$hspHash{$5}=[[$1,$2,$3,$4]];
	    }
	}else{die "error in HHspFormat: input line $_ cannot be parsed\n";}
    }

    return \%hspHash;
}

1;
