CREATE OR REPLACE FUNCTION genomic_sequence(INTEGER) RETURNS TEXT AS $$
BEGIN { 
    strict->import(); 
}

    my ($feature_id) = @_;
    my $locs = spi_exec_query("SELECT l.start, l.stop, l.strand,
             l.chromosome, f.dataset_id FROM location l, feature f WHERE
            f.feature_id = l.feature_id AND f.feature_id =" . $feature_id ."
            ORDER BY l.start;");
    my $nrows = $locs->{processed};
    if(!$nrows){ return "";}

    my $rows = $locs->{rows};
    my $start = $rows->[0]{'start'};
    my $did = $rows->[0]{'dataset_id'};

    # convert the actual start (e.g. 183891) to the 10Kmer
    # that contains it (e.g. 180001)
    my $fstart = $start - ( $start % 10000 ) + 1;
    my $stop = $rows->[-1]{'stop'};
    my $res;

    # an optimization to do the fastest query.
    # only use the IN() clause if it spans multiple 10Kmers
    if ($stop - $fstart > 10000){
        my @starts;
        for(my $i=$fstart;$i<$stop;$i+=10000){
            push(@starts,$i);
        }

        $res = spi_exec_query("SELECT g.sequence_data as sseq
        FROM genomic_sequence g WHERE g.dataset_id = ". $did ." AND g.start IN
        (". join(",",@starts)  .") ORDER BY g.start;");
    } # otherwise, just get the single 10Kmer.
    else{
        $res = spi_exec_query("SELECT g.sequence_data as sseq
        FROM genomic_sequence g WHERE g.dataset_id = ". $did ." AND 
        g.start = ". $fstart . " ORDER BY g.start;");
    }

    # convert from feature space to this sequence which
    # starts at basepair of 1.
    my $s0 = $start - ($rows->[0]{start} % 10000)+ 1;

    my $full_seq = "";
    map { $full_seq .= $_->{sseq} } @{$res->{rows}};
    my $seq = "";
    for(my $i=0;$i<$nrows;$i++){
        $seq .= substr($full_seq 
                      ,($rows->[$i]{start} - $s0)  
                      ,($rows->[$i]{stop}  - $rows->[$i]{start} + 1 ));
    }
    if($rows->[0]{strand} == -1){
        $seq = reverse $seq;
        $seq =~ tr/atcgATCG/tagcTAGC/;
    }
    return  $seq;
$$ LANGUAGE plperl;
