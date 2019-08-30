# nubeam-dedup
Removing PCR duplicates for sequencing reads. **Currently we do not offer license. Use the tool or copy source code by author's permission only**.
## Usage:
To complile, download the * *model.cpp* * and * *Makefile* * to the same directory and type `make`.

`./nubeam-dedup -h` gives you the following messages:
```
./nubeam-dedup [-i -o -i1 -i2 -o1 -o2 -s -h]

Remove exact PCR duplicates for sequencing reads in (gzipped) fastq format.
Produces de-duplicated reads in fastq files with user-given name.

--in or -i: input file name for SE reads
--out or -o: Output file name for SE reads. The default is input file name appended with '.uniq'.
--in1 or -i1: Input file name for PE reads read 1 file.
--in2 or -i2: Input file name for PE reads read 2 file.
--out1 or -o1: Output file name for PE reads read 1 file. The default is read 1 file name appended with '.uniq'.
--out2 or -o2: Output file name for PE reads read 2 file. The default is read 2 file name appended with '.uniq'.
--strand or -s: whether take reads from complementary strand into account. Accept boolean 1 (default) or 0.

-h: print this help
```
## Examples:
- For single-end reads
  - Consider reads from complementary strand (default)
  
    `./nubeam-dedup -i read.fq`
  - Do not consider reads from complementary strand
    
    `./nubeam-dedup -i read.fq -s 0`

- For paired-end reads
  - Consider reads from complementary strand (default)
  
    `./nubeam-dedup -i1 read1.fq -i2 read2.fq`
  - Do not consider reads from complementary strand
  
    `./nubeam-dedup -i1 read1.fq -i2 read2.fq -s 0`
  
