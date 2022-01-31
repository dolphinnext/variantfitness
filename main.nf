$HOSTNAME = ""
params.outdir = 'results'  


if (!params.reads){params.reads = ""} 
if (!params.mate){params.mate = ""} 
if (!params.wtseq){params.wtseq = ""} 
if (!params.metadata){params.metadata = ""} 
if (!params.startpos){params.startpos = ""} 
// Stage empty file to be used as an optional input where required
ch_empty_file_1 = file("$baseDir/.emptyfiles/NO_FILE_1", hidden:true)
ch_empty_file_2 = file("$baseDir/.emptyfiles/NO_FILE_2", hidden:true)
ch_empty_file_3 = file("$baseDir/.emptyfiles/NO_FILE_3", hidden:true)
ch_empty_file_4 = file("$baseDir/.emptyfiles/NO_FILE_4", hidden:true)
ch_empty_file_5 = file("$baseDir/.emptyfiles/NO_FILE_5", hidden:true)
ch_empty_file_6 = file("$baseDir/.emptyfiles/NO_FILE_6", hidden:true)
ch_empty_file_7 = file("$baseDir/.emptyfiles/NO_FILE_7", hidden:true)
ch_empty_file_8 = file("$baseDir/.emptyfiles/NO_FILE_8", hidden:true)
ch_empty_file_9 = file("$baseDir/.emptyfiles/NO_FILE_9", hidden:true)
ch_empty_file_10 = file("$baseDir/.emptyfiles/NO_FILE_10", hidden:true)
ch_empty_file_11 = file("$baseDir/.emptyfiles/NO_FILE_11", hidden:true)

Channel
	.fromFilePairs( params.reads , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.into{g_2_reads_g15_3;g_2_reads_g15_18}

Channel.value(params.mate).into{g_3_mate_g_13;g_3_mate_g15_3;g_3_mate_g15_11;g_3_mate_g15_16;g_3_mate_g15_18;g_3_mate_g15_19;g_3_mate_g15_20;g_3_mate_g15_21}
Channel.value(params.wtseq).into{g_18_barcode_g_17;g_18_barcode_g_24;g_18_barcode_g_25}
g_26_txtFile_g_25 = file(params.metadata, type: 'any')
Channel.value(params.startpos).set{g_27_value_tuple_g_24}

//* params.run_Adapter_Removal =   "no"   //* @dropdown @options:"yes","no" @show_settings:"Adapter_Removal"
//* @style @multicolumn:{seed_mismatches, palindrome_clip_threshold, simple_clip_threshold} @condition:{Tool_for_Adapter_Removal="trimmomatic", seed_mismatches, palindrome_clip_threshold, simple_clip_threshold}, {Tool_for_Adapter_Removal="fastx_clipper", discard_non_clipped}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 1000
    $CPU  = 1
    $MEMORY = 24
    $QUEUE = "long"
}
//* platform
//* autofill
if (!((params.run_Adapter_Removal && (params.run_Adapter_Removal == "yes")) || !params.run_Adapter_Removal)){
g_2_reads_g15_18.set{g15_18_reads_g15_19}
g15_18_log_file_g15_11 = Channel.empty()
} else {


process Adapter_Trimmer_Quality_Module_Adapter_Removal {

input:
 set val(name), file(reads) from g_2_reads_g15_18
 val mate from g_3_mate_g15_18

output:
 set val(name), file("reads/*.fastq")  into g15_18_reads_g15_19
 file "*.{fastx,trimmomatic}.log"  into g15_18_log_file_g15_11

errorStrategy 'retry'

when:
(params.run_Adapter_Removal && (params.run_Adapter_Removal == "yes")) || !params.run_Adapter_Removal

shell:
phred = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.phred
Tool_for_Adapter_Removal = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.Tool_for_Adapter_Removal
Adapter_Sequence = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.Adapter_Sequence
//trimmomatic_inputs
min_length = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.min_length
seed_mismatches = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.seed_mismatches
palindrome_clip_threshold = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.palindrome_clip_threshold
simple_clip_threshold = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.simple_clip_threshold

//fastx_clipper_inputs
discard_non_clipped = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.discard_non_clipped
    
remove_previous_reads = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.remove_previous_reads
discard_non_clipped_text = ""
if (discard_non_clipped == "yes") {discard_non_clipped_text = "-c"}
nameAll = reads.toString()
nameArray = nameAll.split(' ')
file2 = ""
if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}
'''
#!/usr/bin/env perl
 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage;
 use Cwd qw();
 
runCmd("mkdir reads adapter unpaired");

open(OUT, ">adapter/adapter.fa");
my @adaps=split(/\n/,"!{Adapter_Sequence}");
my $i=1;
foreach my $adap (@adaps)
{
 print OUT ">adapter$i\\n$adap\\n";
 $i++;
}
close(OUT);

system("!{runGzip}");
my $quality="!{phred}";
print "fastq quality: $quality\\n";
print "tool: !{Tool_for_Adapter_Removal}\\n";

if ("!{mate}" eq "pair") {
    if ("!{Tool_for_Adapter_Removal}" eq "trimmomatic") {
        runCmd("trimmomatic PE -threads 1 -phred${quality} !{file1} !{file2} reads/!{name}.1.fastq unpaired/!{name}.1.fastq.unpaired reads/!{name}.2.fastq unpaired/!{name}.2.fastq.unpaired ILLUMINACLIP:adapter/adapter.fa:!{seed_mismatches}:!{palindrome_clip_threshold}:!{simple_clip_threshold} MINLEN:!{min_length} 2> !{name}.trimmomatic.log");
    } elsif ("!{Tool_for_Adapter_Removal}" eq "fastx_clipper") {
        print "Fastx_clipper is not suitable for paired reads.";
    }
} else {
    if ("!{Tool_for_Adapter_Removal}" eq "trimmomatic") {
        runCmd("trimmomatic SE -threads 1  -phred${quality} !{file1} reads/!{name}.fastq ILLUMINACLIP:adapter/adapter.fa:!{seed_mismatches}:!{palindrome_clip_threshold}:!{simple_clip_threshold} MINLEN:!{min_length} 2> !{name}.trimmomatic.log");
    } elsif ("!{Tool_for_Adapter_Removal}" eq "fastx_clipper") {
        runCmd("fastx_clipper  -Q $quality -a !{Adapter_Sequence} -l !{min_length} !{discard_non_clipped_text} -v -i !{file1} -o reads/!{name}.fastq > !{name}.fastx.log");
    }
}
if ("!{remove_previous_reads}" eq "true") {
    my $currpath = Cwd::cwd();
    my @paths = (split '/', $currpath);
    splice(@paths, -2);
    my $workdir= join '/', @paths;
    splice(@paths, -1);
    my $inputsdir = join '/', @paths;
    $inputsdir .= "/work";
    print "INFO: inputs reads will be removed if they are located in the $workdir $inputsdir\\n";
    my @listOfFiles = `readlink -e !{file1} !{file2}`;
    foreach my $targetFile (@listOfFiles){
        if (index($targetFile, $workdir) != -1 || index($targetFile, $inputsdir) != -1) {
            system("rm -f $targetFile");
            print "INFO: $targetFile deleted.\\n";
        }
    }
}


##Subroutines
sub runCmd {
    my ($com) = @_;
    my $error = system($com);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}
'''

}
}


//* params.run_Trimmer =   "no"   //* @dropdown @options:"yes","no" @show_settings:"Trimmer"
//* @style @multicolumn:{trim_length_5prime,trim_length_3prime}, {trim_length_5prime_R1,trim_length_3prime_R1}, {trim_length_5prime_R2,trim_length_3prime_R2} @condition:{single_or_paired_end_reads="single", trim_length_5prime,trim_length_3prime}, {single_or_paired_end_reads="pair", trim_length_5prime_R1,trim_length_3prime_R1,trim_length_5prime_R2,trim_length_3prime_R2}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 500
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "long"
}
//* platform
//* autofill
if (!((params.run_Trimmer && (params.run_Trimmer == "yes")) || !params.run_Trimmer)){
g15_18_reads_g15_19.set{g15_19_reads_g15_20}
g15_19_log_file_g15_21 = Channel.empty()
} else {


process Adapter_Trimmer_Quality_Module_Trimmer {

input:
 set val(name), file(reads) from g15_18_reads_g15_19
 val mate from g_3_mate_g15_19

output:
 set val(name), file("reads/*q")  into g15_19_reads_g15_20
 file "*.log" optional true  into g15_19_log_file_g15_21

errorStrategy 'retry'

when:
(params.run_Trimmer && (params.run_Trimmer == "yes")) || !params.run_Trimmer

shell:
phred = params.Adapter_Trimmer_Quality_Module_Trimmer.phred
single_or_paired_end_reads = params.Adapter_Trimmer_Quality_Module_Trimmer.single_or_paired_end_reads
trim_length_5prime = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime
trim_length_3prime = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime
trim_length_5prime_R1 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime_R1
trim_length_3prime_R1 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime_R1
trim_length_5prime_R2 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime_R2
trim_length_3prime_R2 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime_R2
remove_previous_reads = params.Adapter_Trimmer_Quality_Module_Trimmer.remove_previous_reads

nameAll = reads.toString()
nameArray = nameAll.split(' ')
file2 = ""
if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}
'''
#!/usr/bin/env perl
 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage; 
 use Cwd qw();
 
system("mkdir reads");
system("!{runGzip}");
my $file1 = "";
my $file2 = "";
if ("!{mate}" eq "pair") {
    $file1 = "!{file1}";
    $file2 = "!{file2}";
    my $trim1 = "!{trim_length_5prime_R1}:!{trim_length_3prime_R1}";
    my $trim2 = "!{trim_length_5prime_R2}:!{trim_length_3prime_R2}";
    my $len=getLength($file1);
    print "length of $file1: $len\\n";
    trimFiles($file1, $trim1, $len);
    my $len=getLength($file2);
    print "INFO: length of $file2: $len\\n";
    trimFiles($file2, $trim2, $len);
} else {
    $file1 = "!{file1}";
    my $trim1 = "!{trim_length_5prime}:!{trim_length_3prime}";
    my $len=getLength($file1);
    print "INFO: length of file1: $len\\n";
    trimFiles($file1, $trim1, $len);
}
if ("!{remove_previous_reads}" eq "true") {
    my $currpath = Cwd::cwd();
    my @paths = (split '/', $currpath);
    splice(@paths, -2);
    my $workdir= join '/', @paths;
    splice(@paths, -1);
    my $inputsdir= join '/', @paths;
    $inputsdir .= "/inputs";
    print "INFO: inputs reads will be removed if they are located in the workdir inputsdir\\n";
    my @listOfFiles = `readlink -e !{file1} !{file2}`;
    foreach my $targetFile (@listOfFiles){
        if (index($targetFile, $workdir) != -1 || index($targetFile, $inputsdir) != -1) {
            system("rm -f $targetFile");
            print "INFO: $targetFile deleted.\\n";
        }
    }
}



sub trimFiles
{
  my ($file, $trim, $len)=@_;
    my @nts=split(/[,:\\s\\t]+/,$trim);
    my $inpfile="";
    my $com="";
    my $i=1;
    my $outfile="";
    my $param="";
    my $quality="-Q!{phred}";

    if (scalar(@nts)==2)
    {
      $param = "-f ".($nts[0]+1) if (exists($nts[0]) && $nts[0] >= 0 );
      $param .= " -l ".($len-$nts[1]) if (exists($nts[0]) && $nts[1] > 0 );
      $outfile="reads/$file";  
      $com="fastx_trimmer $quality -v $param -o $outfile -i $file > !{name}.fastx_trimmer.log" if ((exists($nts[0]) && $nts[0] > 0) || (exists($nts[0]) && $nts[1] > 0 ));
      print "INFO: $com\\n";
      if ($com eq ""){
          print "INFO: Trimmer skipped for $file \\n";
          system("mv $file reads/.");
      } else {
          runCmd("$com");
          print "INFO: Trimmer executed for $file \\n";
      }
    }

    
}


sub getLength
{
   my ($filename)=@_;
   open (IN, $filename);
   my $j=1;
   my $len=0;
   while(my $line=<IN>)
   {
     chomp($line);
     if ($j >50) { last;}
     if ($j%4==0)
     {
        $len=length($line);
     }
     $j++;
   }
   close(IN);
   return $len;
}

sub runCmd {
    my ($com) = @_;
    my $error = system($com);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}

'''

}
}



process Adapter_Trimmer_Quality_Module_Trimmer_Removal_Summary {

input:
 file logfile from g15_19_log_file_g15_21.collect()
 val mate from g_3_mate_g15_21

output:
 file "trimmer_summary.tsv"  into g15_21_outputFileTSV_g_16

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %tsvDetail;
my %headerHash;
my %headerText;
my %headerTextDetail;

my $i = 0;
chomp( my $contents = `ls *.log` );

my @files = split( /[\\n]+/, $contents );
foreach my $file (@files) {
    $i++;
    my $mapOrder = "1";
    if ($file =~ /(.*)\\.fastx_trimmer\\.log/){
        $file =~ /(.*)\\.fastx_trimmer\\.log/;
        my $mapper   = "fastx_trimmer";
        my $name = $1;    ##sample name
        push( @header, $mapper );
        my $in;
        my $out;
        chomp( $in =`cat $file | grep 'Input:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $out =`cat $file | grep 'Output:' | awk '{sum+=\\$2} END {print sum}'` );

        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Trimmer" ];
    }
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary          = "trimmer_summary.tsv";
writeFile( $summary,          \\%headerText,       \\%tsv );

sub writeFile {
    my $summary    = $_[0];
    my %headerText = %{ $_[1] };
    my %tsv        = %{ $_[2] };
    open( OUT, ">$summary" );
    print OUT "Sample\\t";
    my @headArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push( @headArr, @{ $headerText{$mapOrder} } );
    }
    my $headArrAll = join( "\\t", @headArr );
    print OUT "$headArrAll\\n";

    foreach my $name ( keys %tsv ) {
        my @rowArr = ();
        for my $mapOrder (@sortedOrderArray) {
            push( @rowArr, @{ $tsv{$name}{ $headerHash{$mapOrder} } } );
        }
        my $rowArrAll = join( "\\t", @rowArr );
        print OUT "$name\\t$rowArrAll\\n";
    }
    close(OUT);
}

'''
}

//* params.run_Quality_Filtering =   "no"   //* @dropdown @options:"yes","no" @show_settings:"Quality_Filtering"
//* @style @multicolumn:{window_size,required_quality}, {leading,trailing,minlen}, {minQuality,minPercent} @condition:{tool="trimmomatic", minlen, trailing, leading, required_quality_for_window_trimming, window_size}, {tool="fastx", minQuality, minPercent}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "short"
}
//* platform
//* autofill
if (!((params.run_Quality_Filtering && (params.run_Quality_Filtering == "yes")) || !params.run_Quality_Filtering)){
g15_19_reads_g15_20.set{g15_20_reads_g_13}
g15_20_log_file_g15_16 = Channel.empty()
} else {


process Adapter_Trimmer_Quality_Module_Quality_Filtering {

input:
 set val(name), file(reads) from g15_19_reads_g15_20
 val mate from g_3_mate_g15_20

output:
 set val(name), file("reads/*q")  into g15_20_reads_g_13
 file "*.{fastx,trimmomatic}_quality.log" optional true  into g15_20_log_file_g15_16

errorStrategy 'retry'

when:
(params.run_Quality_Filtering && (params.run_Quality_Filtering == "yes")) || !params.run_Quality_Filtering    

shell:
tool = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.tool
phred = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.phred
window_size = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.window_size
required_quality_for_window_trimming = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.required_quality_for_window_trimming
leading = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.leading
trailing = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.trailing
minlen = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.minlen


// fastx parameters
minQuality = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.minQuality
minPercent = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.minPercent

remove_previous_reads = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.remove_previous_reads

nameAll = reads.toString()
nameArray = nameAll.split(' ')
file2 ="";
if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}
'''
#!/usr/bin/env perl
 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage; 
 use Cwd qw();
 
system("mkdir reads unpaired");
system("!{runGzip}");
my $param = "SLIDINGWINDOW:"."!{window_size}".":"."!{required_quality_for_window_trimming}";
$param.=" LEADING:"."!{leading}";
$param.=" TRAILING:"."!{trailing}";
$param.=" MINLEN:"."!{minlen}";

my $quality="!{phred}";

print "INFO: fastq quality: $quality\\n";
     
if ("!{tool}" eq "trimmomatic") {
    if ("!{mate}" eq "pair") {
        runCmd("trimmomatic PE -phred${quality} !{file1} !{file2} reads/!{name}.1.fastq unpaired/!{name}.1.fastq.unpaired reads/!{name}.2.fastq unpaired/!{name}.1.fastq.unpaired $param 2> !{name}.trimmomatic_quality.log");
    } else {
        runCmd("trimmomatic SE -phred${quality} !{file1} reads/!{name}.fastq $param 2> !{name}.trimmomatic_quality.log");
    }
} elsif ("!{tool}" eq "fastx") {
    if ("!{mate}" eq "pair") {
        print("WARNING: Fastx option is not suitable for paired reads. This step will be skipped.");
        system("mv !{file1} !{file2} reads/.");
    } else {
        runCmd("fastq_quality_filter  -Q $quality -q !{minQuality} -p !{minPercent} -v -i !{file1} -o reads/!{name}.fastq > !{name}.fastx_quality.log");
    }
}
if ("!{remove_previous_reads}" eq "true") {
    my $currpath = Cwd::cwd();
    my @paths = (split '/', $currpath);
    splice(@paths, -2);
    my $workdir= join '/', @paths;
    splice(@paths, -1);
    my $inputsdir= join '/', @paths;
    $inputsdir .= "/inputs";
    print "INFO: inputs reads will be removed if they are located in the workdir inputsdir\\n";
    my @listOfFiles = `readlink -e !{file1} !{file2}`;
    foreach my $targetFile (@listOfFiles){
        if (index($targetFile, $workdir) != -1 || index($targetFile, $inputsdir) != -1) {
            system("rm -f $targetFile");
            print "INFO: $targetFile deleted.\\n";
        }
    }
}

##Subroutines
sub runCmd {
    my ($com) = @_;
    my $error = system($com);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}


'''

}
}



process merge_pairs {

input:
 set val(name),file(reads) from g15_20_reads_g_13
 val mate from g_3_mate_g_13

output:
 set val(name), file("reads/*.fastq")  into g_13_reads_g_14

errorStrategy 'retry'
maxRetries 2

label 'pandaseq'

script:
run_merge = params.run_Merge_Pairs
nameAll = reads.toString()
nameArray = nameAll.split(' ')
def file2;

if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}

"""
#!/bin/bash
mkdir reads
if [ "${mate}" == "pair" ]; then
    if [ "${run_merge}" == "yes" ]; then
        pandaseq -F -N -T 40 -l 50 -f ${file1} -r ${file2} -w reads/${name}.fastq > /dev/null 2>&1  
    else
        cat ${file1} ${file2} > reads/${name}.fastq
    fi
else
    mv $file1 ${name}.fastq > /dev/null 2>&1
    mv ${name}.fastq reads/.
fi

"""
 
}


if (!((params.run_Extract_Barcode && (params.run_Extract_Barcode == "yes")))){
g_13_reads_g_14.set{g_14_reads_g_17}
} else {

process Extract_Barcode {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /ext\/.*.fastq$/) "barcodes/$filename"}
input:
 set val(name),file(reads) from g_13_reads_g_14

output:
 set val(name),file("ext/*.fastq")  into g_14_reads_g_17

errorStrategy 'retry'
maxRetries 3

when:
(params.run_Extract_Barcode && (params.run_Extract_Barcode == "yes"))

shell:
count_reverse_complement = params.count_reverse_complement
five_prime = params.Extract_Barcode.five_prime
three_prime = params.Extract_Barcode.three_prime
min_len = params.Extract_Barcode.min_len
max_len = params.Extract_Barcode.max_len
mismatch = params.Extract_Barcode.mismatch

'''
#!/usr/bin/env perl

my $count_reverse_complement = "!{count_reverse_complement}";
my $five_prime = "!{five_prime}";
my $three_prime = "!{three_prime}";
my $min_len = "!{min_len}";
my $max_len = "!{max_len}";
my $mismatch = "!{mismatch}";

runCmd("mkdir -p ext");

if ($count_reverse_complement eq "no"){
	runCmd("cutadapt -a $three_prime -n $mismatch -o ext/!{name}.tmp !{reads}");
	runCmd("cutadapt -g $five_prime -m $min_len -M $max_len -n $mismatch -o ext/!{name}.fastq ext/!{name}.tmp");
}
else {
	my $rev_threeprime=$three_prime;
	my $rev_fiveprime=$five_prime;
	$rev_threeprime=~tr/ACGT/TGCA/;
	$rev_fiveprime=~tr/ACGT/TGCA/;
	$rev_threeprime = reverse($rev_threeprime);
	$rev_fiveprime = reverse($rev_fiveprime);
	runCmd("cutadapt -a $three_prime -n $mismatch -o ext/!{name}.tmp.1 !{reads}");
	runCmd("cutadapt -g $five_prime -m $min_len -M $max_len -n $mismatch -o ext/!{name}.tmp.2 ext/!{name}.tmp.1");
	runCmd("cutadapt -a $rev_fiveprime -n $mismatch -o ext/!{name}.tmp.3 !{reads}");
	runCmd("cutadapt -g $rev_threeprime -m $min_len -M $max_len -n $mismatch -o ext/!{name}.tmp.4 ext/!{name}.tmp.3");
	runCmd("cat ext/!{name}.tmp.2 ext/!{name}.tmp.4 > ext/!{name}.fastq");
}
##Subroutines
sub runCmd {
    my ($com) = @_;
    my $error = system($com);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}
'''

}
}



process countA {

input:
 set val(name),file(reads) from g_14_reads_g_17
 val wtseq from g_18_barcode_g_17

output:
 set "*.tsv"  into g_17_outputFileTSV_g_20

errorStrategy 'retry'
maxRetries 0
when:
(params.run_countA && (params.run_countA == "yes"))

shell:
'''
#!/usr/bin/env perl

open(READS, "!{reads}");
my $wtseq= "!{wtseq}";
my $name ="!{name}";
my @aalist= qw(A R N D C Q E G H I L K M F P S T W Y V *);
my %table = (
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '*',    # Stop
    'TAG' => '*',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '*',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );

my %count=();
my $wtcount=0;
my $wtcount_aa=0;
my %countaa=();
my %annot=();
my %annotaa=();
my $len=length($wtseq)/3;

# prepare has variables to count R abundance values 
for(my $i=1; $i<=$len; $i++)
{
   $count{$i}=() if (!exists($count{$i}));
   $countaa{$i}=() if (!exists($countaa{$i}));
   $annot{$i}=() if (!exists($annot{$i}));
   $annotaa{$i}=() if (!exists($annotaa{$i}));
   my $wtnt = substr($wtseq, ($i-1)*3, 3);
   my $wtaa= $table{$wtnt};
   foreach my $codon (keys %table){
      $count{$i}{$codon}=0 if (!exists($count{$i}{$codon}));
      $annot{$i}{$codon}=0 if (!exists($annot{$i}{$codon}));
      my $aa = $table{$codon};
      $countaa{$i}{$aa}=0 if (!exists($countaa{$i}{$aa}));
      $annotaa{$i}{$aa}="NA" if (!exists($annotaa{$i}{$aa}));
      my $ann="mut";
      $ann = "syn" if($wtaa eq $aa);
      $ann = "wt" if($wtnt eq $codon);
      $ann = "stop" if("*" eq $aa);
      $annot{$i}{$codon}=$ann;  
      $ann = "wt" if ($ann eq "syn");
      $annotaa{$i}{$aa}=$ann;
  }
}
# count R abundance values per position per codon and aminocacid 
my $j=0;
while($read=<READS>)
{
  chomp($read);
  $j++;
  if ($j % 100000 == 0){
     print "$j reads processed!\n"; 
  }
  if ($read=~/^([ACGT]*)$/) {
     # Find which position has a mutation
     # if there is no mutation icrease WT
     # If there is more than one mutation discard the read
     my $mutpos=0;
     my $mutnt="";
     my $num_of_mut=0;
     for(my $i=1; $i<=$len; $i++)
     {
          my $codon = substr($read, ($i-1)*3, 3);
          my $wtcodon = substr($wtseq, ($i-1)*3, 3);
          if ($codon ne $wtcodon){
          	$mutpos=$i;
          	$num_of_mut++;
          	$mutnt=$codon;
          }
     }
     my $aa = $table{$mutnt};
     if ($num_of_mut == 1){ 
		$count{$mutpos}{$mutnt}++;
		$countaa{$mutpos}{$aa}++;
	 } 
	 if ($num_of_mut == 0) {
		$wtcount++;
		$wtcount_aa++;
	}
  }
}

#Add wild type counts to each position
for(my $i=1; $i<=$len; $i++)
{
	my $wtcodon = substr($wtseq, ($i-1)*3, 3);
	$count{$i}{$wtcodon}=$wtcount;
	my $aa = $table{$wtcodon};
	$countaa{$i}{$aa}=$wtcount_aa;
}

# Get WT counts per position to be able to calculate A relative abundance values using log2(Rmut/Rwt) 
# if any R value is 0, a pseduo count 1 used to calculate log2 for that position


open(OUTNT_R, ">", $name."_nt_R_counts.tsv");
print OUTNT_R "codon\\taa\\tpos\\tannot\\tRcount\n";
open(OUTAA_R, ">", $name."_aa_R_counts.tsv");
print OUTAA_R "aa\\tpos\\tannot\\tRcount\n";
open(OUTNT_A, ">", $name."_nt_A_counts.tsv");
print OUTNT_A "codon\\taa\\tpos\\tannot\\tAcount\n";
open(OUTAA_A, ">", $name."_aa_A_counts.tsv");
print OUTAA_A "aa\\tpos\\tannot\\tAcount\n";
for(my $i=1; $i<=$len; $i++)
{
   $wtcount{$i}=1 if ($wtcount{$i} == 0);
   foreach my $codon (sort keys %table){
      print OUTNT_R $codon."\\t".$table{$codon}."\\t".$i."\\t".$annot{$i}{$codon}."\\t".$count{$i}{$codon}."\n";
      $count{$i}{$codon}=1 if ($count{$i}{$codon}==0);
      print OUTNT_A $codon."\\t".$table{$codon}."\\t".$i."\\t".$annot{$i}{$codon}."\\t".log($count{$i}{$codon}/$wtcount)/log(2)."\n";
   }
   $wtcount_aa{$i}=1 if ($wtcount_aa{$i} == 0);
   foreach my $aa (@aalist){
      print OUTAA_R $aa."\\t".$i."\\t".$annotaa{$i}{$aa}."\\t".$countaa{$i}{$aa}."\n";
      $countaa{$i}{$aa} = 1 if ($countaa{$i}{$aa}==0);
      print OUTAA_A $aa."\\t".$i."\\t".$annotaa{$i}{$aa}."\\t".log($countaa{$i}{$aa}/$wtcount_aa)/log(2)."\n";
   }
}
close(OUTNN_R);
close(OUTAA_R);
close(OUTNN_A);
close(OUTAA_A);
'''
}


process mergeACounts {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /out\/.*.tsv$/) "ARcounts/$filename"}
input:
 file tsv from g_17_outputFileTSV_g_20.collect()

output:
 file "out/*.tsv"  into g_20_outFileTSV0
 file "out/nt_R_counts.tsv"  into g_20_outFileTSV_g_24
 file "out/nt_A_counts.tsv"  into g_20_outFileTSV_g_25

errorStrategy 'retry'
maxRetries 1
	
script:
'''
#!/usr/bin/env perl

use strict;

#################### VARIABLES ######################
 my %tf = (
	"aa_R_counts" => 4,
	"nt_R_counts" => 5,
	"aa_A_counts" => 4,
	"nt_A_counts" => 5
 );
my $indir=$ENV{'PWD'};
my $outdir=$ENV{'PWD'}."/out";
################### PARAMETER PARSING ####################
`mkdir -p $outdir`;
foreach my $type (keys %tf) {
	my $out = "out/".${type}.".tsv";
	opendir I, $indir or die "Could not open ./\\n";
	my @files = grep(/_${type}.tsv\$/, readdir(I));
	closedir I; 
	print "type ".$type." ... col: ".$tf{$type}."\\n";
	my @a=();
	my %b=();
	my $i=0;
	foreach my $f (@files) { 
	  my $ind=$indir."/".$f;
	  my $libname=$f;
	  $libname=~s/_${type}.tsv//;
	  print "$libname\\n";
	  $i++;
	  $a[$i]=$libname;
	  print "${ind}\\n";
	  open IN, $ind or die "Could not open ".$ind."\\n";
	  my $lineNum = 0;
	
	while(<IN>)
	  {
	      next if ($lineNum++ == 0);
	      my @v=split; 
	      #print "Using $type $tf{$type}\\n";
	      my $count = pop @v;
	      $b{join("-", @v)}{$i}= $count;
	  }
	  print "\\n";
	  close IN;
	}
	
	open OUT, ">".${out} or die "Could not open ".${out}." file for writting\\n";
	
	if ($tf{$type} == 4) {
	  print OUT "aa\\tpos\\tannot";
	} else {
	  print OUT "codon\\taa\\tpos\\tannot";
	}
	
	for(my $j=1;$j<=$i;$j++)
	{
	 print OUT "\\t$a[$j]";
	}
	print OUT "\\n";
	
	foreach my $key (keys %b)
	{
	 print OUT replace("-", "\\t", $key); #only print out the first columns
	 
	 for(my $j=1;$j<=$i;$j++)
	 {
	   if (exists($b{$key}) && exists($b{$key}{$j})){
	     print OUT "\\t".$b{$key}{$j};
	   } else {
	     print OUT "\\t0";
	   }
	 }
	 print OUT "\\n";
	}
	close OUT;
}
sub replace {
  my ($from,$to,$string) = @_;
  $string =~s/\$from/\$to/ig;
  return $string;
}
'''
}


process calcFS {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tsv$/) "SFVals/$filename"}
input:
 file nt_A from g_20_outFileTSV_g_25
 file meta from g_26_txtFile_g_25
 val wtseq from g_18_barcode_g_25

output:
 file "*.tsv"  into g_25_outFileTSV0
 file "*AA_SF_vals.tsv"  into g_25_outFileTSV_g_29

errorStrategy 'retry'
maxRetries 0
	
when:
(params.run_calcFS && (params.run_calcFS == "yes"))

shell:
'''
#!/usr/bin/env perl

use strict;
use POSIX qw/ceil/;

#################### VARIABLES ######################
my $indir=$ENV{'PWD'};
my $outdir=$ENV{'PWD'}."/out";
 my %Afile = (
	"codon" => 0,
	"aa" => 1,
	"pos" => 2,
	"annot" => 3
 ); 

 sub median {( sort { $a <=> $b } @_ )[ ceil( $#_/2 ) ] }
 my @aalist= qw(A R N D C Q E G H I L K M F P S T W Y V *);
 my %table = (
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '*',    # Stop
    'TAG' => '*',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '*',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );
################### PARAMETER PARSING ####################
`mkdir -p $outdir`;

open(IN, "!{nt_A}");
open(META, "!{meta}");
my $wtseq= "!{wtseq}";
my $len=length($wtseq)/3;
my @comps = [];

# Read A values and annotations for each libray, position and codon.
# The hash will be Avals{lib}{pos}{codon}
# While reading, read syn and stop for each lib as well into other hash variables per library
# Astop and Asyn will be array into hash per position %Asyn{lib}{pos} = \\@AllSyn and same for stop codon
my %Avals = {};
my %Avalsaa = {};
my %Annots = {};
my %Annotsaa = {};
my %Astop = {};
my %Asyn = {};
my @Aheader = [];
my $i=0;
while(my $row=<IN>){
	my @line = split(/[\\t,\\s]/,$row);
	my $pos = $line[$Afile{"pos"}];
	my $aa = $line[$Afile{"aa"}];
	my $codon = $line[$Afile{"codon"}];
	my $ann = $line[$Afile{"annot"}];
	# line = codon,aa,pos,annot,lib{1},lib{2} ... 
	if ($i==0){
		@Aheader = @line;
		for(my $j=4; $j<=$#Aheader; $j++){
			my $lib = $Aheader[$j];
			$Avals{$lib}=() if (!exists($Avals{$lib}));
			$Avalsaa{$lib}=() if (!exists($Avalsaa{$lib}));
			$Annots{$lib}=() if (!exists($Annots{$lib}));
			$Annotsaa{$lib}=() if (!exists($Annotsaa{$lib}));
			$Astop{$lib}=() if (!exists($Astop{$lib}));
			$Asyn{$lib}=() if (!exists($Asyn{$lib}));
		}
	}else{
		for(my $j=4; $j<=$#Aheader; $j++){
			my $lib = $Aheader[$j];
			
			$Avals{$lib}{$pos}=() if (!exists($Avals{$lib}{$pos}));
			$Avalsaa{$lib}{$pos}=() if (!exists($Avalsaa{$lib}{$pos}));
			$Annots{$lib}{$pos}=() if (!exists($Annots{$lib}{$pos}));
			$Annotsaa{$lib}{$pos}=() if (!exists($Annotsaa{$lib}{$pos}));
			$Astop{$lib}{$pos}=[] if (!exists($Astop{$lib}{$pos}));
			$Asyn{$lib}{$pos}=[] if (!exists($Asyn{$lib}{$pos}));
	
			$Avals{$lib}{$pos}{$codon} = $line[$j];
			$Annots{$lib}{$pos}{$codon} = "$codon\\t$aa\\t$pos\\t$ann";
			my $aa = $table{$codon};
			if (!exists($Avalsaa{$lib}{$pos}{$aa})){
			  $Avalsaa{$lib}{$pos}{$aa} = 2**$line[$j];
			  $ann="wt" if($ann=="syn");
			  $Annotsaa{$lib}{$pos}{$aa} = "$aa\\t$pos\\t$ann";
			}else{
			  $Avalsaa{$lib}{$pos}{$aa} += 2**$line[$j];
			}
			if ($ann eq "syn"){
				push(@{$Asyn{$lib}{$pos}}, $line[$j]);
				print $lib."\\t".$Annots{$lib}{$pos}{$codon}."\\t".$line[$j]."\\n";
			}elsif ($ann eq "stop"){
				push(@{$Astop{$lib}{$pos}}, $line[$j]);
				print $lib."\\t".$Annots{$lib}{$pos}{$codon}."\\t".$line[$j]."\\n";
			}
		}
	}
	$i++;
}

my %Fvals = {};
my %Fvalsaa = {};
my %Fwtsyn = {};
my %Fstop = {};
my %Svals = {};
my %Svalsaa = {};

$i=0;
# Read metadata file
while(my $meta=<META>){
        if($i>0){
        chomp($meta);
        my @line = split(/[\\t,\\s]/,$meta);
		my $P1=$line[0];
		my $P0=$line[1];
		my $comp=$P1."_".$P0;
		open(OUTNT, ">".$comp."_NT_SF_vals.tsv");
		print OUTNT "codon\\taa\\tpos\\tannot\\tFvals\\tSvals\n";
		open(OUTAA, ">".$comp."_AA_SF_vals.tsv");
		print OUTAA "aa\\tpos\\tannot\\tFvals\\tSvals\n";
		
		$Fvals{$comp}=() if (!exists( $Fvals{$comp}));
		$Fvalsaa{$comp}=() if (!exists( $Fvals{$comp}));
		$Fwtsyn{$comp}=() if (!exists( $Fwtsyn{$comp}));
		$Fstop{$comp}=() if (!exists( $Fstop{$comp}));
		
		$Svals{$comp}=()  if (!exists( $Svals{$comp}));
		$Svalsaa{$comp}=()  if (!exists( $Svals{$comp}));
		
		for(my $pos=1; $pos<=$len; $pos++) {
			# Calc median Fwtsyn and Fstop for each position
			$Fwtsyn{$comp}{$pos} = median(@{$Asyn{$P1}{$pos}}) - median(@{$Asyn{$P0}{$pos}});
			$Fstop{$comp}{$pos} = median(@{$Astop{$P1}{$pos}}) - median(@{$Astop{$P0}{$pos}});
			print "$comp\\tFwtsyn\\t".$comp."\\t".$pos."\\t".$Fwtsyn{$comp}{$pos}."\\t". median(@{$Asyn{$P1}{$pos}})." - ".median(@{$Asyn{$P0}{$pos}})."\\n";
			print "$comp\\tFstop\\t".$comp."\\t".$pos."\\t".$Fstop{$comp}{$pos}."\\t". median(@{$Astop{$P1}{$pos}})." - ".median(@{$Astop{$P0}{$pos}})."\\n";
			
			# Calc F and S for all nts
			
			$Fvals{$comp}{$pos}=() if (!exists( $Fvals{$comp}{$pos}));
			$Svals{$comp}{$pos}=() if (!exists( $Svals{$comp}{$pos}));
			foreach my $codon (sort keys %table){
				# F = Ap1 - Ap0
				$Fvals{$comp}{$pos}{$codon} = $Avals{$P1}{$pos}{$codon} - $Avals{$P0}{$pos}{$codon};
				print "$comp\\tAvals\\t".$Annots{$P1}{$pos}{$codon}."\\t".$Avals{$P1}{$pos}{$codon}."\\t".$pos."\\t".$Avals{$P0}{$pos}{$codon}."\\n";
				# S = Fmut - Fwtsyn / Fwtsyn - Fstop
				$Svals{$comp}{$pos}{$codon} = ($Fvals{$comp}{$pos}{$codon} - $Fwtsyn{$comp}{$pos}) / ($Fwtsyn{$comp}{$pos} - $Fstop{$comp}{$pos});
			    print "$comp\\tFvals\\t".$Annots{$P1}{$pos}{$codon}."\\t(".$Fvals{$comp}{$pos}{$codon}." - ".$Fwtsyn{$comp}{$pos}.") / (".$Fwtsyn{$comp}{$pos}." - ".$Fstop{$comp}{$pos}.")\\n";
			    print "$comp\\tFvals\\t".$Annots{$P1}{$pos}{$codon}."\\t".$Fvals{$comp}{$pos}{$codon}."\\t".$pos."\\t".$Svals{$comp}{$pos}{$codon}."\\n";
			    print OUTNT $Annots{$P1}{$pos}{$codon}."\\t".$Fvals{$comp}{$pos}{$codon}."\\t".$Svals{$comp}{$pos}{$codon}."\\n";
			}
			
			# Calc F and S for all aa
			$Fvalsaa{$comp}{$pos}=() if (!exists( $Fvalsaa{$comp}{$pos}));
			$Svalsaa{$comp}{$pos}=() if (!exists( $Svalsaa{$comp}{$pos}));
			foreach my $aa (@aalist){
			    $Fvalsaa{$comp}{$pos}{$aa} = log($Avalsaa{$P1}{$pos}{$aa})/log(2) - log($Avalsaa{$P0}{$pos}{$aa})/log(2);
			    $Svalsaa{$comp}{$pos}{$aa} = ($Fvalsaa{$comp}{$pos}{$aa} - $Fwtsyn{$comp}{$pos}) / ($Fwtsyn{$comp}{$pos} - $Fstop{$comp}{$pos});
			    print OUTAA $Annotsaa{$P1}{$pos}{$aa}."\\t".$Fvalsaa{$comp}{$pos}{$aa}."\\t".$Svalsaa{$comp}{$pos}{$aa}."\\n";
			}
	    }
	close(OUTNT);
    close(OUTAA);

    }
    $i++;
}
'''


}


process fitnessHeatmap {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.pdf$/) "heatmaps/$filename"}
input:
 file SFvals from g_25_outFileTSV_g_29

output:
 file '*.pdf'  into g_29_outputFilePdf0

errorStrategy 'retry'
maxRetries 1

shell:
'''
#!/usr/bin/env Rscript

library(ComplexHeatmap)
library(reshape2)
library(circlize)

# make matrix of values
filename <- "!{SFvals}"
comp <- str_split(filename, "_AA_SF_vals.tsv")
comp <- comp[[1]][1]
datax <- read.table(filename, sep = "\t", header = T)
datax$pos <- datax$pos + startpos - 1
datax_unmelted = dcast(datax, aa~pos)
rownames(datax_unmelted) <- datax_unmelted[,1]
datax_unmelted <- datax_unmelted[,-1]
rownames_datax_unmelted <- rownames(datax_unmelted)

# make matrix of annotations
datax_label = dcast(datax, aa~pos, value.var = "annot")
rownames(datax_label) <- datax_label[,1]
datax_label <- datax_label[,-1]
datax_label <- datax_label[rownames_datax_unmelted,]

# set legend 
at <- seq(min(datax_unmelted), max(datax_unmelted), by = round((max(datax_unmelted)-min(datax_unmelted))/6))
at <- c(at, max(datax_unmelted))
col_fun2 = colorRamp2(breaks = c(-50, 0, 50), colors = c("blue" ,"black", "red"))

# Heatmap
x <- ComplexHeatmap::Heatmap(datax_unmelted, 
     heatmap_legend_param = list(title = paste("Fitness Effect ", comp), legend_direction = "horizontal", color_bar = "continuous",
             at = c(at[1], at[4], at[7]), labels = c("Deleterious","WT-like","Beneficial"), col = col_fun2), 
     show_column_dend = FALSE, show_row_dend = FALSE, 
     cluster_rows = FALSE, cluster_columns = FALSE,
     column_names_rot = 90,
     cell_fun = function(j, i, x, y, width, height, fill) {
         if(datax_label[i,j] == "wt") {
             grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "black", fill = "transparent", lwd=unit(2, "cm")))
         }
     }
)

pdf(paste0(comp, "_heatmap.pdf"), height=8, width=3,  paper = "a4r", compress=TRUE)
draw(x, heatmap_legend_side = "top")
graphics.off()

'''
}


process codonPlot {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.pdf$/) "codonPlots/$filename"}
input:
 file tsvfile from g_20_outFileTSV_g_24
 val wtseq from g_18_barcode_g_24
 val startpos from g_27_value_tuple_g_24

output:
 file "*.pdf"  into g_24_outputFilePdf0

errorStrategy 'retry'
maxRetries 1

shell:
'''
#!/usr/bin/env Rscript
library(plyr)
library(dplyr)
library(data.table)
library(reshape)
library(ggplot2)

Rcounts <- read.delim("!{tsvfile}")
startpos <- !{startpos}
wtseq <- "!{wtseq}"
Rc_sorted <- Rcounts[with(Rcounts, order(pos, codon)),]

libs <- colnames(Rcounts)[5:ncol(Rcounts)]

for (lib in libs) {

	p  <- list()
	pdf(paste0("barfreq",lib,".pdf"), height=20, width=20,  paper = "a4r", compress=TRUE)
	for (pos in seq(1:(nchar(wtseq)/3))) {
	
		Rc_sorted_pos <- Rc_sorted[Rc_sorted$pos==pos,c("pos", "aa", "codon", "annot", lib)]
		names(Rc_sorted_pos) <-  c("pos", "aa", "codon", "annot", "count")
		freq.by_aa <- data.table(ddply(Rc_sorted_pos, c("pos","aa"), summarize, Freq = sum(count)))
		
		freq.aa_table <- dcast(freq.by_aa , aa ~ pos, value.var="Freq")
		
		aa.totals <- freq.by_aa
		names(aa.totals) <- c("pos", "AA", "Freq.AA")
		codon.totals <- Rc_sorted_pos[, c("codon", "aa", "pos", "annot", "count")]
		names(codon.totals) <- c("Codon", "AA", "pos", "annot", "Freq")
		codon.totals$Codon <- paste0(codon.totals$Codon, "\\n(",codon.totals$annot,")")
		aa.codon.totals <- codon.totals %>% right_join(aa.totals, by=c("AA", "pos"))
		aa.codon.totals$Codon.Pct.Freq <- with(aa.codon.totals, Freq/Freq.AA)
		
		label_aa_count <- function(aa) paste(aa, aa.totals[aa.totals$AA == aa, "Freq.AA"], 
		        sep = ": ")
		
		aa.codon.totals$AA.Display <- unlist(lapply(aa.codon.totals$AA, label_aa_count))
		
		p[[pos]] <- ggplot(aa.codon.totals, 
		  aes(y = Codon.Pct.Freq, x = Codon, label = Freq)) + geom_bar(stat="identity") + 
		  ylim(0, 1.2) +
		  facet_wrap(~AA.Display, scale = "free_x") + 
		  geom_text(colour = "black", vjust = -0.5) + theme_bw() +
		  ggtitle(paste0(lib," pos:", pos + startpos - 1, " codon frequency plot")) +
		  theme(plot.title = element_text(hjust = 0.5)) +
		  xlab("Codon") + ylab("Codon % Freq.")
		
		Rc_p<- Rc_sorted[Rc_sorted$pos==pos & Rc_sorted$annot!="wt" ,]
		Rc_p_freq <- Rc_p
		Rc_p_freq[,lib] <- Rc_p[,lib]/sum( Rc_p[,lib])
        barplot(Rc_p_freq[,lib], names.arg = Rc_p_freq[,"codon"], 
		    main=paste0(lib," pos:", pos, " codon frequency plot"))

	}
	graphics.off()
	pdf(paste0("bar",lib,".pdf"), height=20, width=20,  paper = "a4r", compress=TRUE)
	for (pos in seq(1:(nchar(wtseq)/3))) {
	  print(p[[pos]])
	}
	graphics.off()
	
}

'''
}


process Adapter_Trimmer_Quality_Module_Quality_Filtering_Summary {

input:
 file logfile from g15_20_log_file_g15_16.collect()
 val mate from g_3_mate_g15_16

output:
 file "quality_filter_summary.tsv"  into g15_16_outputFileTSV_g_16

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %headerHash;
my %headerText;

my $i = 0;
chomp( my $contents = `ls *.log` );
my @files = split( /[\\n]+/, $contents );
foreach my $file (@files) {
    $i++;
    my $mapper   = "";
    my $mapOrder = "1";
    if ($file =~ /(.*)\\.fastx_quality\\.log/){
        $mapper   = "fastx";
        $file =~ /(.*)\\.fastx_quality\\.log/;
        my $name = $1;    ##sample name
        push( @header, $mapper );
        my $in;
        my $out;
        chomp( $in =`cat $file | grep 'Input:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $out =`cat $file | grep 'Output:' | awk '{sum+=\\$2} END {print sum}'` );
        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Quality Filtering" ];
    } elsif ($file =~ /(.*)\\.trimmomatic_quality\\.log/){
        $mapper   = "trimmomatic";
        $file =~ /(.*)\\.trimmomatic_quality\\.log/;
        my $name = $1;    ##sample name
        push( @header, $mapper );
        my $in;
        my $out;
        if ( "!{mate}" eq "pair"){
            chomp( $in =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$4} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$7} END {print sum}'` );
        } else {
            chomp( $in =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$3} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$5} END {print sum}'` );
        }
        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Quality Filtering" ];
    }
    
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary          = "quality_filter_summary.tsv";
writeFile( $summary,          \\%headerText,       \\%tsv );

sub writeFile {
    my $summary    = $_[0];
    my %headerText = %{ $_[1] };
    my %tsv        = %{ $_[2] };
    open( OUT, ">$summary" );
    print OUT "Sample\\t";
    my @headArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push( @headArr, @{ $headerText{$mapOrder} } );
    }
    my $headArrAll = join( "\\t", @headArr );
    print OUT "$headArrAll\\n";

    foreach my $name ( keys %tsv ) {
        my @rowArr = ();
        for my $mapOrder (@sortedOrderArray) {
            push( @rowArr, @{ $tsv{$name}{ $headerHash{$mapOrder} } } );
        }
        my $rowArrAll = join( "\\t", @rowArr );
        print OUT "$name\\t$rowArrAll\\n";
    }
    close(OUT);
}

'''
}


process Adapter_Trimmer_Quality_Module_Adapter_Removal_Summary {

input:
 file logfile from g15_18_log_file_g15_11.collect()
 val mate from g_3_mate_g15_11

output:
 file "adapter_removal_summary.tsv"  into g15_11_outputFileTSV_g_16
 file "adapter_removal_detailed_summary.tsv" optional true  into g15_11_outputFile1

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %tsvDetail;
my %headerHash;
my %headerText;
my %headerTextDetail;

my $i = 0;
chomp( my $contents = `ls *.log` );

my @files = split( /[\\n]+/, $contents );
foreach my $file (@files) {
    $i++;
    my $mapOrder = "1";
    if ($file =~ /(.*)\\.fastx\\.log/){
        $file =~ /(.*)\\.fastx\\.log/;
        my $mapper   = "fastx";
        my $name = $1;    ##sample name
        push( @header, $mapper );

        my $in;
        my $out;
        my $tooshort;
        my $adapteronly;
        my $noncliped;
        my $Nreads;

        chomp( $in =`cat $file | grep 'Input:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $out =`cat $file | grep 'Output:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $tooshort =`cat $file | grep 'too-short reads' | awk '{sum+=\\$2} END {print sum}'`);
        chomp( $adapteronly =`cat $file | grep 'adapter-only reads' | awk '{sum+=\\$2} END {print sum}'`);
        chomp( $noncliped =`cat $file | grep 'non-clipped reads.' | awk '{sum+=\\$2} END {print sum}'`);
        chomp( $Nreads =`cat $file | grep 'N reads.' | awk '{sum+=\\$2} END {print sum}'` );

        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Adapter Removal" ];
        $tsvDetail{$name}{$mapper} = [ $in, $tooshort, $adapteronly, $noncliped, $Nreads, $out ];
        $headerTextDetail{$mapOrder} = ["Total Reads","Too-short reads","Adapter-only reads","Non-clipped reads","N reads","Reads After Adapter Removal"];
    } elsif ($file =~ /(.*)\\.trimmomatic\\.log/){
        $file =~ /(.*)\\.trimmomatic\\.log/;
        my $mapper   = "trimmomatic";
        my $name = $1;    ##sample name
        push( @header, $mapper );
        
        my $in;
        my $out;

        if ( "!{mate}" eq "pair"){
            chomp( $in =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$4} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$7} END {print sum}'` );
        } else {
            chomp( $in =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$3} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$5} END {print sum}'` );
        }
        


        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Adapter Removal" ];
        
    }
    
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary          = "adapter_removal_summary.tsv";
my $detailed_summary = "adapter_removal_detailed_summary.tsv";
writeFile( $summary,          \\%headerText,       \\%tsv );
if (%headerTextDetail){
    writeFile( $detailed_summary, \\%headerTextDetail, \\%tsvDetail );  
}

sub writeFile {
    my $summary    = $_[0];
    my %headerText = %{ $_[1] };
    my %tsv        = %{ $_[2] };
    open( OUT, ">$summary" );
    print OUT "Sample\\t";
    my @headArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push( @headArr, @{ $headerText{$mapOrder} } );
    }
    my $headArrAll = join( "\\t", @headArr );
    print OUT "$headArrAll\\n";

    foreach my $name ( keys %tsv ) {
        my @rowArr = ();
        for my $mapOrder (@sortedOrderArray) {
            push( @rowArr, @{ $tsv{$name}{ $headerHash{$mapOrder} } } );
        }
        my $rowArrAll = join( "\\t", @rowArr );
        print OUT "$name\\t$rowArrAll\\n";
    }
    close(OUT);
}

'''
}

g15_11_outputFileTSV_g_16= g15_11_outputFileTSV_g_16.ifEmpty([""]) 
g15_21_outputFileTSV_g_16= g15_21_outputFileTSV_g_16.ifEmpty([""]) 
g15_16_outputFileTSV_g_16= g15_16_outputFileTSV_g_16.ifEmpty([""]) 

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 30
    $CPU  = 1
    $MEMORY = 10
    $QUEUE = "short"
}
//* platform
//* autofill

process Overall_Summary {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /overall_summary.tsv$/) "summary/$filename"}
input:
 file adapterSum from g15_11_outputFileTSV_g_16
 file trimmerSum from g15_21_outputFileTSV_g_16
 file qualitySum from g15_16_outputFileTSV_g_16

output:
 file "overall_summary.tsv"  into g_16_outputFileTSV0

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_rows;
my @seen_cols;
my $ID_header;

chomp(my $contents = `ls *.tsv`);
my @rawFiles = split(/[\\n]+/, $contents);
my @files = ();
# order must be in this order for chipseq pipeline: bowtie->dedup
# rsem bam pipeline: dedup->rsem, star->dedup
# riboseq ncRNA_removal->star
my @order = ("adapter_removal","trimmer","quality","extractUMI","extractValid","tRAX","sequential_mapping","ncRNA_removal","bowtie","star","hisat2","tophat2", "dedup","rsem","kallisto","esat","count");
for ( my $k = 0 ; $k <= $#order ; $k++ ) {
    for ( my $i = 0 ; $i <= $#rawFiles ; $i++ ) {
        if ( $rawFiles[$i] =~ /$order[$k]/ ) {
            push @files, $rawFiles[$i];
        }
    }
}

print Dumper \\@files;
##add rest of the files
for ( my $i = 0 ; $i <= $#rawFiles ; $i++ ) {
    push(@files, $rawFiles[$i]) unless grep{$_ == $rawFiles[$i]} @files;
}
print Dumper \\@files;

##Merge each file according to array order

foreach my $file (@files){
        open IN,"$file";
        my $line1 = <IN>;
        chomp($line1);
        ( $ID_header, my @header) = ( split("\\t", $line1) );
        push @seen_cols, @header;

        while (my $line=<IN>) {
        chomp($line);
        my ( $ID, @fields ) = ( split("\\t", $line) ); 
        my %this_row;
        @this_row{@header} = @fields;

        #print Dumper \\%this_row;

        foreach my $column (@header) {
            if (! exists $all_rows{$ID}{$column}) {
                $all_rows{$ID}{$column} = $this_row{$column}; 
            }
        }   
    }
    close IN;
}

#print for debugging
#print Dumper \\%all_rows;
#print Dumper \\%seen_cols;

#grab list of column headings we've seen, and order them. 
my @cols_to_print = uniq(@seen_cols);
my $summary = "overall_summary.tsv";
open OUT, ">$summary";
print OUT join ("\\t", $ID_header,@cols_to_print),"\\n";
foreach my $key ( keys %all_rows ) { 
    #map iterates all the columns, and gives the value or an empty string. if it's undefined. (prevents errors)
    print OUT join ("\\t", $key, (map { $all_rows{$key}{$_} // '' } @cols_to_print)),"\\n";
}
close OUT;

sub uniq {
    my %seen;
    grep ! $seen{$_}++, @_;
}

'''


}

//* params.run_FastQC =  "no"  //* @dropdown @options:"yes","no" @description:"FastQC provides quality control checks on raw sequence data."



process Adapter_Trimmer_Quality_Module_FastQC {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.(html|zip)$/) "fastQC/$filename"}
input:
 val mate from g_3_mate_g15_3
 set val(name), file(reads) from g_2_reads_g15_3

output:
 file '*.{html,zip}'  into g15_3_FastQCout0

errorStrategy 'retry'
maxRetries 3

when:
(params.run_FastQC && (params.run_FastQC == "yes"))

script:
nameAll = reads.toString()
if (nameAll.contains('.gz')) {
    file =  nameAll - '.gz' - '.gz'
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    file =  nameAll 
    runGzip = ''
}
"""
${runGzip}
fastqc ${file} 
"""
}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
