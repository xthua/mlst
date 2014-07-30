use Bio::SearchIO; 
use Bio::SeqIO; 
my $refresh_database=0;
my $mlst_file="ab.fas";

sub gene_length{
	my $length;
	my @list=@_; # input gene_name
	my $inseq = Bio::SeqIO->new(-file => "$mlst_file", format => 'fasta');
	while(my $seq = $inseq->next_seq){
		if($seq->display_id eq $list[0]){
			$length=length($seq->seq);
			#print $list[0]."$length\n";

		}
	}
	return($length);

}
opendir DIR, "." or die $!;
while((my $filename=readdir(DIR))){
	my $n=0;
	if($filename=~/fna$/){
		#print $filename;
		my @line=split(/\./,$filename);
		my $database=$line[0];
		#print $filename,$database;
		if($refresh_database){
			system("makeblastdb -in $filename -dbtype nucl -out $database");
		}
		# blast
		system(" blastn -db $database  -query $mlst_file -out $database.out");
		# read blast result and write seq to $database.fas
		open FH, ">$database.fa";
		my $in = new Bio::SearchIO(-format => 'blast', 
                           -file   => "$database.out");
		   while( my $result = $in->next_result ) {
  ## $result is a Bio::Search::Result::ResultI compliant object
  while( my $hit = $result->next_hit ) {
    ## $hit is a Bio::Search::Hit::HitI compliant object
    while( my $hsp = $hit->next_hsp ) {
      ## $hsp is a Bio::Search::HSP::HSPI compliant object
      if( $hsp->length('total') > 200 ) {
        if ( $hsp->percent_identity >= 75 && $n <= 7) {
		my $new_string=length($hsp->hit_string);
		my $query_name=$result->query_name;
		#print "$new_string\t",gene_length($query_name),"$query_name\n";
		my $gene_length=gene_length($query_name);
		if($gene_length <= $new_string){
          print FH  ">",   $result->query_name,"\n",
	  $hsp->hit_string,"\n";
	  $n++;
		}else{
			print "not eq\n";
			my $start=$hsp->start('hit');
			$start=$start-10;
			my $end=$start+$gene_length+10;
			#print "$start\t$end\n";
			my $file="$database.fna";
			#print $file;
			my $string;
			my $inseq = Bio::SeqIO->new(-file => "$file", format => 'fasta');
			while(my $seq = $inseq->next_seq){
				$string=$seq->subseq($start,$end);
				print length($string),"$new_string\t";
			}
			  print FH  ">",   $result->query_name,length($string),"\n",
			  $string,"\n";

		}
#            " Hit=",        $hit->name,
#            " Length=",     $hsp->length('total'),
#            " Percent_id=", $hsp->percent_identity, "\n";
        }
      }
    }  
  }
}
close(FH);

	}
}
close(DIR);
