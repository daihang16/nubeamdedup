# nubeamdedup
Removing PCR duplicates for sequencing reads. **Currently we do not offer license. Use the tool or copy source code by author's permission only**.
## Usage:
To complile, download the * *model.cpp* * and * *Makefile* * to the same directory and type `make`.

`./nubeamdedup -h` gives you the following messages:
```
./nubeamdedup [-i -o -i1 -i2 -o1 -o2 -s -h]

Remove exact PCR duplicates for sequencing reads in (gzipped) fastq format.
Produces de-duplicated reads in fastq files with user-given name.

--in or -i: input file name for SE reads
--out or -o: output file name for SE reads
--in1 or -i1: input file name for PE reads read 1 file
--in2 or -i2: input file name for PE reads read 2 file
--out1 or -o1: output file name for PE reads read 1 file
--out2 or -o2: output file name for PE reads read 2 file
--strand or -s: whether take reads from complementary strand into account. Accept boolean 0 (default) or 1.

-h: print this help
```
## Examples:
- For single-end reads
  - Do not consider reads from complementary strand
    
    `./nubeamdedup -i read.fq -o read.uniq.fq`
  - Consider reads from complementary strand
  
    `./nubeamdedup -i read.fq -o read.uniq.fq -s 1`

- For paired-end reads
  - Do not consider reads from complementary strand
  
    `./nubeamdedup -i1 read1.fq -i2 read2.fq -o1 read1.uniq.fq -o2 read2.uniq.fq`
  - Consider reads from complementary strand
  
    `./nubeamdedup -i1 read1.fq -i2 read2.fq -o1 read1.uniq.fq -o2 read2.uniq.fq -s 1`
