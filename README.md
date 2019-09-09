# nubeam-dedup
Removing exact PCR duplicates for sequencing reads. **Currently we do not offer license. Use the tool or copy source code by author's permission only**.
## Compiling:
### Dependency
`zlib` is required to compile. To install `zlib`, run the following commands:

`wget https://www.zlib.net/zlib1211.zip`

`unzip zlib1211.zip`

`cd zlib-1.2.11/`

`./configure`

`make`

`make install`
### Compile nubeam-dedup
Run the following commands:

`wget https://github.com/daihang16/nubeamdedup/archive/master.zip`

`unzip master.zip`

`cd nubeamdedup-master/Linux/` or `cd nubeamdedup-master/macOS/`

`make && make clean`

`./nubeam-dedup -i1 ../toydata/1.fq.gz -i2 ../toydata/2.fq.gz`

`wc -l *.fastq`

You should see both output files have 276884 lines.

We also offer pre-compliled executable file for Linux. The executable file was compliled on Red Hat Enterprise Linux Server 7.0 (Maipo) by compiler gcc with the version of 4.8.2 20140120 (Red Hat 4.8.2-16). C++11 was used.

## Usage:
`./nubeam-dedup -h` gives you the following messages:
```
./nubeam-dedup [-i -o -d    -i1 -i2 -o1 -o2 -d1 -d2    -s -r -z -h]

Remove exact PCR duplicates for sequencing reads in (gzipped) fastq format.
Produces de-duplicated reads in fastq files.

Parameters for single-end (SE) reads:
--in or -i: Input file name. The parameter is mandatory for SE reads.
--out or -o: Output file name for unique reads. The default is input file name prefix appended with '.uniq.fastq(.gz)', under the current directory.
--duplicate or -d: File name for removed duplicated reads. The parameter is only valid when --remove or -r is set as 1 (see below). The default is input file name prefix appended with '.removed.fastq(.gz)', under the current directory.

Parameters for paired-end (PE) reads:
--in1 or -i1: Input file name for read 1 file. The parameter is mandatory for PE reads.
--in2 or -i2: Input file name for read 2 file. The parameter is mandatory for PE reads.
--out1 or -o1: Output file name for unique read pairs read 1 file. The default is read 1 file name prefix appended with '.uniq.fastq(.gz)', under the current working directory.
--out2 or -o2: Output file name for unique read pairs read 2 file. The default is read 2 file name prefix appended with '.uniq.fastq(.gz)', under the current working directory.
--duplicate1 or -d1: File name for removed duplicated read pairs read 1 file. The parameter is only valid when --remove or -r is set as 1 (see below). The default is input file name prefix appended with '.removed.fastq(.gz)', under the current directory.
--duplicate2 or -d2: File name for removed duplicated read pairs read 2 file. The parameter is only valid when --remove or -r is set as 1 (see below). The default is input file name prefix appended with '.removed.fastq(.gz)', under the current directory.

Miscellaneous parameters:
--strand or -s: Whether take reads from complementary strand into account. Accept boolean 1 (default) or 0.
--remove or -r: Whether output removed duplicated reads. Accept boolean 0 (default) or 1.
--gz or -z: Compression level of output file. Accept integer 0 (default) to 9. If 0 (default), the output data will not be compressed and will be written to plain text file; otherwise, the output data will be written to gzip format file, with the compression level suggested by user. If compression is needed, a compression level of 6 is recommended as a compromise between speed and compression.

-h: print this help
```
## Examples:
- For single-end reads
  - Consider reads from complementary strand (default)
  
    `./nubeam-dedup -i read.fq`
    
    The command gives the following output on screen:
    
    `Output unique reads to /current/working/directory/read.uniq.fastq`
  - Do not consider reads from complementary strand
    
    `./nubeam-dedup -i read.fq -s 0`
  - Consider reads from complementary strand (default), output gzipped file with a compression level of 6, output removed duplicated reads 
  
    `./nubeam-dedup -i read.fq -z 6 -r 1`
    
    The command gives the following output on screen:
    
    `Output removed duplicated reads to /current/working/directory/read.removed.fastq.gz`
    
    `Output unique reads to /current/working/directory/read.uniq.fastq.gz`

- For paired-end reads
  - Consider reads from complementary strand (default)
  
    `./nubeam-dedup -i1 read1.fastq.gz -i2 read2.fastq.gz`
    
    The command gives the following output on screen:
    
    `Output unique read pairs read 1 to /current/working/directory/reads1.uniq.fastq`
    
    `Output unique read pairs read 2 to /current/working/directory/reads2.uniq.fastq`
  - Do not consider reads from complementary strand
  
    `./nubeam-dedup -i1 read1.fastq.gz -i2 read2.fastq.gz -s 0`
  - Consider reads from complementary strand (default), output gzipped file with a compression level of 6, output removed duplicated reads
  
    `./nubeam-dedup -i1 read1.fastq.gz -i2 read2.fastq.gz -z 6 -r 1`
    
    The command gives the following output on screen:
    
    `Output removed duplicated read pairs read 1 to /current/working/directory/read1.removed.fastq.gz`
    
    `Output removed duplicated read pairs read 2 to /current/working/directory/read2.removed.fastq.gz`
    
    `Output unique read pairs read 1 to /current/working/directory/reads1.uniq.fastq.gz`
    
    `Output unique read pairs read 2 to /current/working/directory/reads2.uniq.fastq.gz`
