# nubeam-dedup
`nubeam-dedup` is a fast and easy-to-use bioinformatics tool removing exact PCR duplicates for sequencing reads, single-end or paired-end. We appreciate your interest in `nubeam-dedup`. If you use `nubeam-dedup`, please kindly cite:

Hang Dai and Yongtao Guan, **Nubeam-dedup: a fast and RAM-efficient tool to de-duplicate sequencing reads without mapping.** *Bioinformatics* 36(10), P.3254-3256, (2020) DOI: [10.1093/bioinformatics/btaa112](https://doi.org/10.1093/bioinformatics/btaa112)

## Compiling:
<!---
### Dependency
`zlib` is required to compile. For most computers it was already installed. If not, run the following commands to install `zlib`:

On `Linux`:
```Shell
wget https://www.zlib.net/zlib1211.zip
```
On `macOS`:
```Shell
curl -o zlib1211.zip https://www.zlib.net/zlib1211.zip
```
Then:
```Shell
unzip zlib1211.zip
cd zlib-1.2.11/
./configure
make
sudo make install
```
-->
### Compile nubeam-dedup
Run the following commands:

On `Linux`:
```console
foo@bar:~$ wget --no-check-certificate --content-disposition https://github.com/daihang16/nubeamdedup/archive/master.zip
foo@bar:~$ unzip nubeamdedup-master.zip
foo@bar:~$ cd nubeamdedup-master/Linux/
```
On `macOS`:
```console
foo@bar:~$ curl -LJO https://github.com/daihang16/nubeamdedup/archive/master.zip
foo@bar:~$ unzip nubeamdedup-master.zip
foo@bar:~$ cd nubeamdedup-master/macOS/
```
Then:
```console
foo@bar:Linux$ make && make clean
foo@bar:Linux$ ./nubeam-dedup -i1 ../toydata/1.fq.gz -i2 ../toydata/2.fq.gz 1> out.txt 2> log.txt
foo@bar:Linux$ cat log.txt
Output unique read pairs read 1 to nubeamdedup/Linux/1.uniq.fastq
Output unique read pairs read 2 to nubeamdedup/Linux/2.uniq.fastq
foo@bar:Linux$ cat out.txt
69221/142250 read pairs are unique.
foo@bar:Linux$ wc -l *.fastq
276884 1.uniq.fastq
276884 2.uniq.fastq
553768 total

```

You should see the expected output as above.

We also offer pre-compiled executable file for Linux. The executable file was compiled on Ubuntu 18.04.2 LTS by compiler gcc with the version of 7.4.0 (Ubuntu 7.4.0-1ubuntu1~18.04). C++11 was used.

[//]: # 'We also offer pre-compliled executable file for Linux. The executable file was compiled on Red Hat Enterprise Linux Server 7.0 (Maipo) by compiler gcc with the version of 4.8.2 20140120 (Red Hat 4.8.2-16). C++11 was used.'


## Usage:
`./nubeam-dedup -h` gives you the following messages:
```console
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
--gz or -z: Compression level of output file. Accept integer 0 (default) to 9. If 0 (default), the output data will not be compressed and will be written to plain text file; otherwise, the output data will be written to gzip format file, with the compression level suggested by user. If compression is needed, a compression level less than 3 is recommended as a compromise between speed and compression.

-h: print this help
```
## Examples:
- For single-end reads
  - Consider reads from complementary strand (default)
  
    ```Shell
    ./nubeam-dedup -i read.fq
    ```
    
    The command gives the following output on screen:
    
    ```console
    Output unique reads to /current/working/directory/read.uniq.fastq
    x/y reads are unique.
    ```
  - Do not consider reads from complementary strand
    
    ```Shell
    ./nubeam-dedup -i read.fq -s 0
    ```
  - Consider reads from complementary strand (default), output gzipped file with a compression level of 6, output removed duplicated reads 
  
    ```Shell
    ./nubeam-dedup -i read.fq -z 6 -r 1
    ```
    
    The command gives the following output on screen:
    
    ```console
    Output removed duplicated reads to /current/working/directory/read.removed.fastq.gz    
    Output unique reads to /current/working/directory/read.uniq.fastq.gz
    x/y reads are unique.
    ```

- For paired-end reads
  - Consider reads from complementary strand (default)
  
    ```Shell
    ./nubeam-dedup -i1 read1.fastq.gz -i2 read2.fastq.gz
    ```
    
    The command gives the following output on screen:
    
    ```console
    Output unique read pairs read 1 to /current/working/directory/reads1.uniq.fastq    
    Output unique read pairs read 2 to /current/working/directory/reads2.uniq.fastq
    x/y read pairs are unique.
    ```
  - Do not consider reads from complementary strand
  
    ```Shell
    ./nubeam-dedup -i1 read1.fastq.gz -i2 read2.fastq.gz -s 0
    ```
  - Consider reads from complementary strand (default), output gzipped file with a compression level of 2, output removed duplicated reads
  
    ```Shell
    ./nubeam-dedup -i1 read1.fastq.gz -i2 read2.fastq.gz -z 2 -r 1
    ```
    
    The command gives the following output on screen:
    
    ```console
    Output removed duplicated read pairs read 1 to /current/working/directory/read1.removed.fastq.gz
    Output removed duplicated read pairs read 2 to /current/working/directory/read2.removed.fastq.gz    
    Output unique read pairs read 1 to /current/working/directory/reads1.uniq.fastq.gz    
    Output unique read pairs read 2 to /current/working/directory/reads2.uniq.fastq.gz
    x/y read pairs are unique.
    ```
## Miscellaneous:    
- A large value (like 6) for `-z` tag might significantly increase the running time. From Figures 1, 2 and 7 in [this post](https://clearlinux.org/news-blogs/linux-os-data-compression-options-comparing-behavior), `-z 6` would increase the amount of time by a factor of 2.5-3 compared with `-z 1` (with a limited gain regarding to compression ratio); and `-z 1` would increase the amount of time by a factor of 2.5 compared with `-z 0`, which is the default setting of `nubeam-dedup`. The recommended practice is: either use a smaller compression level (1-3) or do not use the `-z` tag at all. For the latter choice, if compression was required, [pigz](https://zlib.net/pigz/) could be used after `nubeam-dedup` finishes---this can significantly accelerate the compression.
- For the convenience of users, `nubeam-dedup` outputs the number of unique reads and total reads to `stdout` and the output file information to `stderr`. Use `1>` and `2>` to [redirect the two streams](https://www.howtogeek.com/435903/what-are-stdin-stdout-and-stderr-on-linux/) respectively. 
