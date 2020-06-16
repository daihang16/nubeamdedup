#include <fstream> 
#include <iostream> 
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set> 
#include <utility>
#include <zlib.h>
#include <string>
#include <cmath> 
#include <stdio.h>
#include <iterator>
#include <unistd.h>
#include <sysexits.h>

using namespace std; 

// utility functions
void operation(double prod[][2], bool yes); 
void matrix_multiplication_helper(double prod[][2], char * raw_seq, char nucleotide);
char * reverse_complement(char * raw_seq);

// do not consider reads from complementary strand
void quantify_reads(string fin, string fout, string compression_level); // SE reads
void quantify_reads(string fin1, string fin2, string fout1, string fout2, string compression_level); // PE reads
void quantify_reads(string fin, string fout, string fremoved, string compression_level); // SE reads, output removed reads
void quantify_reads(string fin1, string fin2, string fout1, string fout2, string fremoved1, string fremoved2, string compression_level); // PE reads, output removed reads

// consider reads from complementary strand
void quantify_reads_rc(string fin, string fout, string compression_level); // SE reads
void quantify_reads_rc(string fin1, string fin2, string fout1, string fout2, string compression_level); // PE reads
void quantify_reads_rc(string fin, string fout, string fremoved, string compression_level); // SE reads, output removed reads
void quantify_reads_rc(string fin1, string fin2, string fout1, string fout2, string fremoved1, string fremoved2, string compression_level); // PE reads, output removed reads

