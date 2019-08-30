// /*************************************************************************
// * 
// * Duke University CONFIDENTIAL
// * __________________
// * 
// * [2019] - [20**] Duke University & Yongtao Guan & Hang Dai
// * All Rights Reserved.
// * 
// * NOTICE: All information contained herein is, and remains
// * the property of Duke University and the creator. 
// * The intellectual and technical concepts contained
// * herein are proprietary to Duke University and the creator
// * and may be covered by U.S. and Foreign Patents, patents in process, 
// * and are protected by trade secret or copyright law.
// * Dissemination of this information or reproduction of this material
// * is strictly forbidden unless prior written permission is obtained
// * from Duke University.

#include <fstream> 
#include <iostream> 
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <limits>
#include <string.h>
#include <zlib.h>
#include <string>
#include <cmath> 
#include <vector> 
#include <map>
#include <set> 
#include <stdio.h>
#include <sys/vfs.h>
#include <iterator>
	 
using namespace std; 
						  
#define dd 2 
void operation(double prod[][dd], bool yes); 
void matrix_multiplication_helper(double prod[][dd], char * raw_seq, char nucleotide);
char * reverse_complement(char * raw_seq);
void quantify_reads(string fin, string fout); // SE reads
void quantify_reads(string fin1, string fin2, string fout1, string fout2); // PE reads
void quantify_reads_rc(string fin, string fout); // SE reads
void quantify_reads_rc(string fin1, string fin2, string fout1, string fout2); // PE reads
void quantify_reads_approx(string fin1, string fin2, string fout1, string fout2); // PE reads
void quantify_reads_rc_approx(string fin1, string fin2, string fout1, string fout2); // PE reads

// SE reads
void quantify_reads(string fin, string fout) 
{
	gzFile output_file = gzopen(fout.c_str(), "wT"); 
	if(output_file == NULL) 
	{
		printf("can't open %s file to write\n", fout.c_str()); 
		exit(0); 
	}

	gzFile input_file = gzopen(fin.c_str(), "r"); 
	if(input_file == NULL) {
		printf("can't open %s file to read. \n", fin.c_str()); 
		exit(0); 
	}
	
	char * seq_id = new char[200];  
	char * raw_sequence = new char[200]; 
	char * seq_id_repeat = new char[200];
	char * tscore = new char[200];  

	unordered_set<double> vdat;

	double read_identifier;
	while (gzgets(input_file, seq_id, 200) &&
		gzgets(input_file, raw_sequence, 200) &&
		gzgets(input_file, seq_id_repeat, 200) &&
		gzgets(input_file, tscore, 200)) {
		double prod[dd][dd] = {{1,0},{0,1}};
		matrix_multiplication_helper(prod, raw_sequence, 'A');
		matrix_multiplication_helper(prod, raw_sequence, 'T');
		matrix_multiplication_helper(prod, raw_sequence, 'C');
		matrix_multiplication_helper(prod, raw_sequence, 'G');
		read_identifier = prod[0][0] + sqrt(3) * prod[0][1] + M_SQRT2 * prod[1][0] + sqrt(5) * prod[1][1];
		// check if the read has already appeared before
		// if unseen, keep it; else, skip it 
		if (vdat.insert(read_identifier).second == 1) {
			gzprintf(output_file, "%s%s+\n%s", seq_id, raw_sequence, tscore);
		}
	}
	gzclose(input_file);
	gzclose(output_file); 

	vdat.clear();
	delete[] seq_id;
	delete[] raw_sequence; 
	delete[] seq_id_repeat;
	delete[] tscore; 
}

// SE reads
void quantify_reads_rc(string fin, string fout) 
{
	gzFile output_file = gzopen(fout.c_str(), "wT"); 
	if(output_file == NULL) 
	{
		printf("can't open %s file to write\n", fout.c_str()); 
		exit(0); 
	}

	gzFile input_file = gzopen(fin.c_str(), "r"); 
	if(input_file == NULL) {
		printf("can't open %s file to read. \n", fin.c_str()); 
		exit(0); 
	}
	
	char * seq_id = new char[200];  
	char * raw_sequence = new char[200]; 
	char * seq_id_repeat = new char[200];
	char * tscore = new char[200];  

	unordered_set<double> vdat;

	double read_identifier;
	double read_identifier_rc;
	while (gzgets(input_file, seq_id, 200) &&
		gzgets(input_file, raw_sequence, 200) &&
		gzgets(input_file, seq_id_repeat, 200) &&
		gzgets(input_file, tscore, 200)) {
		double prod[dd][dd] = {{1,0},{0,1}};
		matrix_multiplication_helper(prod, raw_sequence, 'A');
		matrix_multiplication_helper(prod, raw_sequence, 'T');
		matrix_multiplication_helper(prod, raw_sequence, 'C');
		matrix_multiplication_helper(prod, raw_sequence, 'G');		
		read_identifier = prod[0][0] + sqrt(3) * prod[0][1] + M_SQRT2 * prod[1][0] + sqrt(5) * prod[1][1];		
		// check if the read has already appeared before
		// if unseen, keep it; else, skip it 
		if (vdat.insert(read_identifier).second == 1) {
			char * raw_sequence_rc = reverse_complement(raw_sequence);
			double prod_rc[dd][dd] = {{1,0},{0,1}};
			matrix_multiplication_helper(prod_rc, raw_sequence_rc, 'A');
			matrix_multiplication_helper(prod_rc, raw_sequence_rc, 'T');
			matrix_multiplication_helper(prod_rc, raw_sequence_rc, 'C');
			matrix_multiplication_helper(prod_rc, raw_sequence_rc, 'G');
			read_identifier_rc = prod_rc[0][0] + sqrt(3) * prod_rc[0][1] + M_SQRT2 * prod_rc[1][0] + sqrt(5) * prod_rc[1][1];
			vdat.insert(read_identifier_rc);
			gzprintf(output_file, "%s%s+\n%s", seq_id, raw_sequence, tscore);
			delete[] raw_sequence_rc;
		}
	}
	gzclose(input_file);
	gzclose(output_file); 

	vdat.clear();
	delete[] seq_id;
	delete[] raw_sequence; 
	delete[] seq_id_repeat;
	delete[] tscore; 
}

// PE reads
void quantify_reads(string fin1, string fin2, string fout1, string fout2) 
{
	gzFile output_file1 = gzopen(fout1.c_str(), "wT"); 
	if(output_file1 == NULL) 
	{
		printf("can't open %s file to write\n", fout1.c_str()); 
		exit(0); 
	}
	gzFile output_file2 = gzopen(fout2.c_str(), "wT"); 
	if(output_file2 == NULL) 
	{
		printf("can't open %s file to write\n", fout2.c_str()); 
		exit(0); 
	}

	gzFile input_file1 = gzopen(fin1.c_str(), "r"); 
	if(input_file1 == NULL) {
		printf("can't open %s file to read. \n", fin1.c_str()); 
		exit(0); 
	}
	gzFile input_file2 = gzopen(fin2.c_str(), "r"); 
	if(input_file2 == NULL) {
		printf("can't open %s file to read. \n", fin2.c_str()); 
		exit(0); 
	}
	
	char * seq_id1 = new char[200];  
	char * raw_sequence1 = new char[200]; 
	char * seq_id_repeat1 = new char[200];
	char * tscore1 = new char[200];  

	char * seq_id2 = new char[200];  
	char * raw_sequence2 = new char[200]; 
	char * seq_id_repeat2 = new char[200];
	char * tscore2 = new char[200];  

	unordered_multimap<double, unsigned int> vdat1;
	unordered_multimap<double, unsigned int> vdat2;

	set<unsigned int> temp_set;

	double read_identifier1;
	double read_identifier2;

	unsigned int line_num = 0;
	while (gzgets(input_file1, seq_id1, 200) &&
		gzgets(input_file1, raw_sequence1, 200) &&
		gzgets(input_file1, seq_id_repeat1, 200) &&
		gzgets(input_file1, tscore1, 200) &&
		gzgets(input_file2, seq_id2, 200) &&
		gzgets(input_file2, raw_sequence2, 200) &&
		gzgets(input_file2, seq_id_repeat2, 200) &&
		gzgets(input_file2, tscore2, 200)) {
		line_num++;
		double prod1[dd][dd] = {{1,0},{0,1}};
		double prod2[dd][dd] = {{1,0},{0,1}};
		matrix_multiplication_helper(prod1, raw_sequence1, 'A');
		matrix_multiplication_helper(prod1, raw_sequence1, 'T');
		matrix_multiplication_helper(prod1, raw_sequence1, 'C');
		matrix_multiplication_helper(prod1, raw_sequence1, 'G');
		matrix_multiplication_helper(prod2, raw_sequence2, 'A');
		matrix_multiplication_helper(prod2, raw_sequence2, 'T');
		matrix_multiplication_helper(prod2, raw_sequence2, 'C');
		matrix_multiplication_helper(prod2, raw_sequence2, 'G');
		read_identifier1 = prod1[0][0] + sqrt(3) * prod1[0][1] + M_SQRT2 * prod1[1][0] + sqrt(5) * prod1[1][1];
		read_identifier2 = prod2[0][0] + sqrt(3) * prod2[0][1] + M_SQRT2 * prod2[1][0] + sqrt(5) * prod2[1][1];
		// check if the read pair has already appeared before
		// if unseen, keep it; else, skip it 
		auto pos1 = vdat1.equal_range(read_identifier1);
		auto pos2 = vdat2.equal_range(read_identifier2);
		if (pos1.first == vdat1.end() || pos2.first == vdat2.end()) {
			vdat1.insert(pos1.second, pair<double, unsigned int>(read_identifier1, line_num)); // gives it hint
			vdat2.insert(pos2.second, pair<double, unsigned int>(read_identifier2, line_num));
			gzprintf(output_file1, "%s%s+\n%s", seq_id1, raw_sequence1, tscore1);
			gzprintf(output_file2, "%s%s+\n%s", seq_id2, raw_sequence2, tscore2);
		} else {
			while (pos1.first != pos1.second) {
				temp_set.insert(pos1.first->second);
				pos1.first++;
			}
			unsigned int count = 0;
			while (count == 0 && pos2.first != pos2.second) {
				count = temp_set.count(pos2.first->second);
				pos2.first++;
			}
			if (count == 0) {
				vdat1.insert(pos1.second, pair<double, unsigned int>(read_identifier1, line_num)); // gives it hint
				vdat2.insert(pos2.second, pair<double, unsigned int>(read_identifier2, line_num));
				gzprintf(output_file1, "%s%s+\n%s", seq_id1, raw_sequence1, tscore1);
				gzprintf(output_file2, "%s%s+\n%s", seq_id2, raw_sequence2, tscore2);
			}
			temp_set.clear();
		}
	}
	gzclose(input_file1);
	gzclose(input_file2);
	gzclose(output_file1);
	gzclose(output_file2); 

	vdat1.clear();
	vdat2.clear();
	delete[] seq_id1;
	delete[] raw_sequence1; 
	delete[] seq_id_repeat1;
	delete[] tscore1; 
	delete[] seq_id2;
	delete[] raw_sequence2; 
	delete[] seq_id_repeat2;
	delete[] tscore2; 
}

// PE reads
void quantify_reads_rc(string fin1, string fin2, string fout1, string fout2) 
{
	gzFile output_file1 = gzopen(fout1.c_str(), "wT"); 
	if(output_file1 == NULL) 
	{
		printf("can't open %s file to write\n", fout1.c_str()); 
		exit(0); 
	}
	gzFile output_file2 = gzopen(fout2.c_str(), "wT"); 
	if(output_file2 == NULL) 
	{
		printf("can't open %s file to write\n", fout2.c_str()); 
		exit(0); 
	}

	gzFile input_file1 = gzopen(fin1.c_str(), "r"); 
	if(input_file1 == NULL) {
		printf("can't open %s file to read. \n", fin1.c_str()); 
		exit(0); 
	}
	gzFile input_file2 = gzopen(fin2.c_str(), "r"); 
	if(input_file2 == NULL) {
		printf("can't open %s file to read. \n", fin2.c_str()); 
		exit(0); 
	}
	
	char * seq_id1 = new char[200];  
	char * raw_sequence1 = new char[200]; 
	char * seq_id_repeat1 = new char[200];
	char * tscore1 = new char[200];  

	char * seq_id2 = new char[200];  
	char * raw_sequence2 = new char[200]; 
	char * seq_id_repeat2 = new char[200];
	char * tscore2 = new char[200];  

	unordered_multimap<double, unsigned int> vdat1;
	unordered_multimap<double, unsigned int> vdat2;

	set<unsigned int> temp_set;

	double read_identifier1;
	double read_identifier2;

	unsigned int line_num = 0;
	while (gzgets(input_file1, seq_id1, 200) &&
		gzgets(input_file1, raw_sequence1, 200) &&
		gzgets(input_file1, seq_id_repeat1, 200) &&
		gzgets(input_file1, tscore1, 200) &&
		gzgets(input_file2, seq_id2, 200) &&
		gzgets(input_file2, raw_sequence2, 200) &&
		gzgets(input_file2, seq_id_repeat2, 200) &&
		gzgets(input_file2, tscore2, 200)) {
		line_num++;
		double prod1[dd][dd] = {{1,0},{0,1}};
		double prod2[dd][dd] = {{1,0},{0,1}};
		matrix_multiplication_helper(prod1, raw_sequence1, 'A');
		matrix_multiplication_helper(prod1, raw_sequence1, 'T');
		matrix_multiplication_helper(prod1, raw_sequence1, 'C');
		matrix_multiplication_helper(prod1, raw_sequence1, 'G');
		matrix_multiplication_helper(prod2, raw_sequence2, 'A');
		matrix_multiplication_helper(prod2, raw_sequence2, 'T');
		matrix_multiplication_helper(prod2, raw_sequence2, 'C');
		matrix_multiplication_helper(prod2, raw_sequence2, 'G');
		read_identifier1 = prod1[0][0] + sqrt(3) * prod1[0][1] + M_SQRT2 * prod1[1][0] + sqrt(5) * prod1[1][1];
		read_identifier2 = prod2[0][0] + sqrt(3) * prod2[0][1] + M_SQRT2 * prod2[1][0] + sqrt(5) * prod2[1][1];
		// check if the read pair has already appeared before
		// if unseen, keep it; else, skip it 
		auto pos1 = vdat1.equal_range(read_identifier1);
		auto pos2 = vdat2.equal_range(read_identifier2);
		if (pos1.first == vdat1.end() || pos2.first == vdat2.end()) {
			vdat1.insert(pos1.second, pair<double, unsigned int>(read_identifier1, line_num)); // gives it hint
			vdat2.insert(pos2.second, pair<double, unsigned int>(read_identifier2, line_num));
			vdat1.insert({read_identifier2, line_num});
			vdat2.insert({read_identifier1, line_num});
			gzprintf(output_file1, "%s%s+\n%s", seq_id1, raw_sequence1, tscore1);
			gzprintf(output_file2, "%s%s+\n%s", seq_id2, raw_sequence2, tscore2);
		} else {
			while (pos1.first != pos1.second) {
				temp_set.insert(pos1.first->second);
				pos1.first++;
			}
			unsigned int count = 0;
			while (count == 0 && pos2.first != pos2.second) {
				count = temp_set.count(pos2.first->second);
				pos2.first++;
			}
			if (count == 0) {
				vdat1.insert(pos1.second, pair<double, unsigned int>(read_identifier1, line_num)); // gives it hint
				vdat2.insert(pos2.second, pair<double, unsigned int>(read_identifier2, line_num));
				vdat1.insert({read_identifier2, line_num});
				vdat2.insert({read_identifier1, line_num});
				gzprintf(output_file1, "%s%s+\n%s", seq_id1, raw_sequence1, tscore1);
				gzprintf(output_file2, "%s%s+\n%s", seq_id2, raw_sequence2, tscore2);
			}
			temp_set.clear();
		}
	}
	gzclose(input_file1);
	gzclose(input_file2);
	gzclose(output_file1);
	gzclose(output_file2); 

	vdat1.clear();
	vdat2.clear();
	delete[] seq_id1;
	delete[] raw_sequence1; 
	delete[] seq_id_repeat1;
	delete[] tscore1; 
	delete[] seq_id2;
	delete[] raw_sequence2; 
	delete[] seq_id_repeat2;
	delete[] tscore2; 
}

// PE reads
void quantify_reads_approx(string fin1, string fin2, string fout1, string fout2) 
{
	gzFile output_file1 = gzopen(fout1.c_str(), "wT"); 
	if(output_file1 == NULL) 
	{
		printf("can't open %s file to write\n", fout1.c_str()); 
		exit(0); 
	}
	gzFile output_file2 = gzopen(fout2.c_str(), "wT"); 
	if(output_file2 == NULL) 
	{
		printf("can't open %s file to write\n", fout2.c_str()); 
		exit(0); 
	}

	gzFile input_file1 = gzopen(fin1.c_str(), "r"); 
	if(input_file1 == NULL) {
		printf("can't open %s file to read. \n", fin1.c_str()); 
		exit(0); 
	}
	gzFile input_file2 = gzopen(fin2.c_str(), "r"); 
	if(input_file2 == NULL) {
		printf("can't open %s file to read. \n", fin2.c_str()); 
		exit(0); 
	}
	
	char * seq_id1 = new char[200];  
	char * raw_sequence1 = new char[200]; 
	char * seq_id_repeat1 = new char[200];
	char * tscore1 = new char[200];  

	char * seq_id2 = new char[200];  
	char * raw_sequence2 = new char[200]; 
	char * seq_id_repeat2 = new char[200];
	char * tscore2 = new char[200]; 

	unordered_set<double> vdat;

	double read_identifier;

	while (gzgets(input_file1, seq_id1, 200) &&
		gzgets(input_file1, raw_sequence1, 200) &&
		gzgets(input_file1, seq_id_repeat1, 200) &&
		gzgets(input_file1, tscore1, 200) &&
		gzgets(input_file2, seq_id2, 200) &&
		gzgets(input_file2, raw_sequence2, 200) &&
		gzgets(input_file2, seq_id_repeat2, 200) &&
		gzgets(input_file2, tscore2, 200)) {
		double prod[dd][dd] = {{1,0},{0,1}};
		matrix_multiplication_helper(prod, raw_sequence1, 'A');
		matrix_multiplication_helper(prod, raw_sequence1, 'T');
		matrix_multiplication_helper(prod, raw_sequence1, 'C');
		matrix_multiplication_helper(prod, raw_sequence1, 'G');
		matrix_multiplication_helper(prod, raw_sequence2, 'A');
		matrix_multiplication_helper(prod, raw_sequence2, 'T');
		matrix_multiplication_helper(prod, raw_sequence2, 'C');
		matrix_multiplication_helper(prod, raw_sequence2, 'G');
		read_identifier = prod[0][0] + sqrt(3) * prod[0][1] + M_SQRT2 * prod[1][0] + sqrt(5) * prod[1][1];
		// check if the read pair has already appeared before
		// if unseen, keep it; else, skip it 
		if (vdat.insert(read_identifier).second == 1) {
			gzprintf(output_file1, "%s%s+\n%s", seq_id1, raw_sequence1, tscore1);
			gzprintf(output_file2, "%s%s+\n%s", seq_id2, raw_sequence2, tscore2);
		}
	}
	gzclose(input_file1);
	gzclose(input_file2);
	gzclose(output_file1);
	gzclose(output_file2); 

	vdat.clear();
	delete[] seq_id1;
	delete[] raw_sequence1; 
	delete[] seq_id_repeat1;
	delete[] tscore1; 
	delete[] seq_id2;
	delete[] raw_sequence2; 
	delete[] seq_id_repeat2;
	delete[] tscore2; 
}

// PE reads
void quantify_reads_rc_approx(string fin1, string fin2, string fout1, string fout2) 
{
	gzFile output_file1 = gzopen(fout1.c_str(), "wT"); 
	if(output_file1 == NULL) 
	{
		printf("can't open %s file to write\n", fout1.c_str()); 
		exit(0); 
	}
	gzFile output_file2 = gzopen(fout2.c_str(), "wT"); 
	if(output_file2 == NULL) 
	{
		printf("can't open %s file to write\n", fout2.c_str()); 
		exit(0); 
	}

	gzFile input_file1 = gzopen(fin1.c_str(), "r"); 
	if(input_file1 == NULL) {
		printf("can't open %s file to read. \n", fin1.c_str()); 
		exit(0); 
	}
	gzFile input_file2 = gzopen(fin2.c_str(), "r"); 
	if(input_file2 == NULL) {
		printf("can't open %s file to read. \n", fin2.c_str()); 
		exit(0); 
	}
	
	char * seq_id1 = new char[200];  
	char * raw_sequence1 = new char[200]; 
	char * seq_id_repeat1 = new char[200];
	char * tscore1 = new char[200];  

	char * seq_id2 = new char[200];  
	char * raw_sequence2 = new char[200]; 
	char * seq_id_repeat2 = new char[200];
	char * tscore2 = new char[200]; 

	unordered_set<double> vdat;

	double read_identifier;
	double read_identifier_rc;
	while (gzgets(input_file1, seq_id1, 200) &&
		gzgets(input_file1, raw_sequence1, 200) &&
		gzgets(input_file1, seq_id_repeat1, 200) &&
		gzgets(input_file1, tscore1, 200) &&
		gzgets(input_file2, seq_id2, 200) &&
		gzgets(input_file2, raw_sequence2, 200) &&
		gzgets(input_file2, seq_id_repeat2, 200) &&
		gzgets(input_file2, tscore2, 200)) {
		double prod[dd][dd] = {{1,0},{0,1}};
		matrix_multiplication_helper(prod, raw_sequence1, 'A');
		matrix_multiplication_helper(prod, raw_sequence1, 'T');
		matrix_multiplication_helper(prod, raw_sequence1, 'C');
		matrix_multiplication_helper(prod, raw_sequence1, 'G');
		matrix_multiplication_helper(prod, raw_sequence2, 'A');
		matrix_multiplication_helper(prod, raw_sequence2, 'T');
		matrix_multiplication_helper(prod, raw_sequence2, 'C');
		matrix_multiplication_helper(prod, raw_sequence2, 'G');
		read_identifier = prod[0][0] + sqrt(3) * prod[0][1] + M_SQRT2 * prod[1][0] + sqrt(5) * prod[1][1];
		// check if the read pair has already appeared before
		// if unseen, keep it; else, skip it 
		if (vdat.insert(read_identifier).second == 1) {
			double prod_rc[dd][dd] = {{1,0},{0,1}};
			matrix_multiplication_helper(prod_rc, raw_sequence2, 'A');
			matrix_multiplication_helper(prod_rc, raw_sequence2, 'T');
			matrix_multiplication_helper(prod_rc, raw_sequence2, 'C');
			matrix_multiplication_helper(prod_rc, raw_sequence2, 'G');
			matrix_multiplication_helper(prod_rc, raw_sequence1, 'A');
			matrix_multiplication_helper(prod_rc, raw_sequence1, 'T');
			matrix_multiplication_helper(prod_rc, raw_sequence1, 'C');
			matrix_multiplication_helper(prod_rc, raw_sequence1, 'G');
			read_identifier_rc = prod_rc[0][0] + sqrt(3) * prod_rc[0][1] + M_SQRT2 * prod_rc[1][0] + sqrt(5) * prod_rc[1][1];
			vdat.insert(read_identifier_rc);
			gzprintf(output_file1, "%s%s+\n%s", seq_id1, raw_sequence1, tscore1);
			gzprintf(output_file2, "%s%s+\n%s", seq_id2, raw_sequence2, tscore2);
		}
	}
	gzclose(input_file1);
	gzclose(input_file2);
	gzclose(output_file1);
	gzclose(output_file2); 

	vdat.clear();
	delete[] seq_id1;
	delete[] raw_sequence1; 
	delete[] seq_id_repeat1;
	delete[] tscore1; 
	delete[] seq_id2;
	delete[] raw_sequence2; 
	delete[] seq_id_repeat2;
	delete[] tscore2; 
}

void operation(double prod[][dd], bool yes) 
{
	if (yes) {
		prod[0][1] += prod[0][0]; 
		prod[1][1] += prod[1][0]; 
	}
	else {
		prod[0][0] += prod[0][1]; 
		prod[1][0] += prod[1][1]; 
	}
}

void matrix_multiplication_helper(double prod[][dd], char * raw_seq, char nucleotide) {
	while (*raw_seq != '\n') {
		operation(prod, *raw_seq == nucleotide);
		raw_seq++;
	}
}

char * reverse_complement(char * raw_seq) {
	unsigned short int read_length = 0;
	while (*raw_seq != '\n') {
		read_length++;
		raw_seq++;
	}
	// now raw_seq points to '\n', read_length is the length of read
	char * raw_seq_rc = new char [read_length + 1];
	// go backward
	for (unsigned short int i = 0; i < read_length; i++) {
		switch (*(--raw_seq)) {
			case 'A':
				raw_seq_rc[i] = 'T';
				break;
			case 'T':
				raw_seq_rc[i] = 'A';
				break;
			case 'C':
				raw_seq_rc[i] = 'G';
				break;
			case 'G':
				raw_seq_rc[i] = 'C';
				break;
			default:
				raw_seq_rc[i] = 'N';
				break;
		}
	}
	raw_seq_rc[read_length] = '\n';
	return(raw_seq_rc);
}

int main(int argc, char ** argv)
{
	string fin; 
	string fout; 
	string fin1; 
	string fin2; 
	string fout1; 
	string fout2;
	string fnlog = "\0"; 
	bool strand = 0; // whether consider reads from complementary strand or not
	bool qtf = 0; // let it run or not 

	for (int i = 1; i < argc; i++) {
		string str;
		
		if (i>1 && argv[i][0] != '-') 
			continue;
		str.assign(argv[i]);
		if (str.compare("-h") == 0) {
			printf("./nubeamdedup [-i -o -i1 -i2 -o1 -o2 -s -h]\n"); 
			printf("Remove exact PCR duplicates for sequencing reads in (gzipped) fastq format.\n");
			printf("Produces de-duplicated reads in fastq files with user-given name.\n"); 
			printf("--in or -i: input file name for SE reads\n"); 
			printf("--out or -o: output file name for SE reads\n"); 
			printf("--in1 or -i1: input file name for PE reads read 1 file\n"); 
			printf("--in2 or -i2: input file name for PE reads read 2 file\n"); 
			printf("--out1 or -o1: output file name for PE reads read 1 file\n");
			printf("--out2 or -o2: output file name for PE reads read 2 file\n");
			printf("--strand or -s: whether take reads from complementary strand into account. Accept boolean 0 (default) or 1.\n");
			printf("-h: print this help\n"); 
			exit(0); 
		}
		else if (str.compare("--in") == 0 || str.compare("-i") == 0) {
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				fin.clear();
				fin.assign(argv[i+1]);
				qtf = 1;
		}
		else if (str.compare("--out") == 0 || str.compare("-o") == 0) {
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				fout.clear();
				fout.assign(argv[i+1]);
				fnlog.assign(fout); 
		}
		else if (str.compare("--in1") == 0 || str.compare("-i1") == 0) {
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				fin1.clear();
				fin1.assign(argv[i+1]);
				qtf = 1;
		}
		else if (str.compare("--in2") == 0 || str.compare("-i2") == 0) {
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				fin2.clear();
				fin2.assign(argv[i+1]);
		}
		else if (str.compare("--out1") == 0 || str.compare("-o1") == 0) {
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				fout1.clear();
				fout1.assign(argv[i+1]);
				fnlog.assign(fout1); 
		}
		else if (str.compare("--out2") == 0 || str.compare("-o2") == 0) {
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				fout2.clear();
				fout2.assign(argv[i+1]);
		}
		else if (str.compare("--strand") == 0 || str.compare("-s") == 0) {
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				if(!isdigit(argv[i+1][0]))
				{
					printf("wrong augument after option.\n");
					exit(0);
				}
				strand = atoi(argv[i+1]);
		}
		else {
			fprintf(stderr,"Bad option %s\n", argv[i]);
			exit(0);
		}
	}
	
	if (qtf) {
		if (fin1 == "") { // SE reads
			if (strand) {
				quantify_reads_rc(fin, fout);
			} else {
				quantify_reads(fin, fout);
			}
		} else { // PE reads
			if (strand) {
				quantify_reads_rc_approx(fin1, fin2, fout1, fout2);
			} else {
				quantify_reads(fin1, fin2, fout1, fout2);
			}
		}
		return 1; 
	} else {
		printf("No input file(s). Use ./nubeamdedup -h for more options.\n"); 
		exit(0);
	}
	return 0; 
}
