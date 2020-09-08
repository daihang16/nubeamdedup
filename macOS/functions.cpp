#include "functions.h"

// SE reads
void quantify_reads(string fin, string fout, string compression_level) 
{
	gzFile output_file = gzopen(fout.c_str(), compression_level.c_str());
	if (output_file == NULL) {
		printf("can't open %s file to write. \n", fout.c_str()); 
		exit(EX_CANTCREAT); 
	}

	gzFile input_file = gzopen(fin.c_str(), "r"); 
	if (input_file == NULL) {
		printf("can't open %s file to read. \n", fin.c_str()); 
		exit(EX_NOINPUT); 
	}
	
	char * seq_id = new char[310];  
	char * raw_sequence = new char[310]; 
	char * seq_id_repeat = new char[310];
	char * tscore = new char[310];  

	unordered_set<double> vdat;

	double read_identifier = 0;

	unsigned int line_num = 0; // number of reads
	unsigned int n_uniq_reads = 0; // number of unique reads
	while (gzgets(input_file, seq_id, 310) &&
		gzgets(input_file, raw_sequence, 310) &&
		gzgets(input_file, seq_id_repeat, 310) &&
		gzgets(input_file, tscore, 310)) {
		line_num++;
		double prod[2][2] = {{1,0},{0,1}};
		matrix_multiplication_helper(prod, raw_sequence, 'A');
		matrix_multiplication_helper(prod, raw_sequence, 'T');
		matrix_multiplication_helper(prod, raw_sequence, 'C');
		matrix_multiplication_helper(prod, raw_sequence, 'G');
		read_identifier = prod[0][0] + sqrt(3) * prod[0][1] + M_SQRT2 * prod[1][0] + sqrt(5) * prod[1][1];
		// check if the read has already appeared before
		// if unseen, keep it; else, skip it 
		if (vdat.insert(read_identifier).second == 1) {
			gzprintf(output_file, "%s%s+\n%s", seq_id, raw_sequence, tscore);
			n_uniq_reads++;
		}
	}
	gzclose(input_file);
	gzclose(output_file); 
	cout << n_uniq_reads << "/" << line_num << " reads are unique." << endl; 

	vdat.clear();
	delete[] seq_id;
	delete[] raw_sequence; 
	delete[] seq_id_repeat;
	delete[] tscore; 
}

// SE reads, output removed reads
void quantify_reads(string fin, string fout, string fremoved, string compression_level) 
{
	gzFile output_file = gzopen(fout.c_str(), compression_level.c_str());
	if (output_file == NULL) {
		printf("can't open %s file to write\n", fout.c_str()); 
		exit(EX_CANTCREAT); 
	}

	gzFile removed_file = gzopen(fremoved.c_str(), compression_level.c_str());
	if (removed_file == NULL) {
		printf("can't open %s file to write\n", fremoved.c_str()); 
		exit(EX_CANTCREAT);
	}

	gzFile input_file = gzopen(fin.c_str(), "r"); 
	if (input_file == NULL) {
		printf("can't open %s file to read. \n", fin.c_str()); 
		exit(EX_NOINPUT); 
	}
	
	char * seq_id = new char[310];  
	char * raw_sequence = new char[310]; 
	char * seq_id_repeat = new char[310];
	char * tscore = new char[310];  

	unordered_set<double> vdat;

	double read_identifier = 0;

	unsigned int line_num = 0; // number of reads
	unsigned int n_uniq_reads = 0; // number of unique reads
	while (gzgets(input_file, seq_id, 310) &&
		gzgets(input_file, raw_sequence, 310) &&
		gzgets(input_file, seq_id_repeat, 310) &&
		gzgets(input_file, tscore, 310)) {
		line_num++;
		double prod[2][2] = {{1,0},{0,1}};
		matrix_multiplication_helper(prod, raw_sequence, 'A');
		matrix_multiplication_helper(prod, raw_sequence, 'T');
		matrix_multiplication_helper(prod, raw_sequence, 'C');
		matrix_multiplication_helper(prod, raw_sequence, 'G');
		read_identifier = prod[0][0] + sqrt(3) * prod[0][1] + M_SQRT2 * prod[1][0] + sqrt(5) * prod[1][1];
		// check if the read has already appeared before
		// if unseen, keep it; else, skip it 
		if (vdat.insert(read_identifier).second == 1) {
			gzprintf(output_file, "%s%s+\n%s", seq_id, raw_sequence, tscore);
			n_uniq_reads++;
		} else {
			gzprintf(removed_file, "%s%s+\n%s", seq_id, raw_sequence, tscore);
		}
	}
	gzclose(input_file);
	gzclose(removed_file);
	gzclose(output_file); 
	cout << n_uniq_reads << "/" << line_num << " reads are unique." << endl; 

	vdat.clear();
	delete[] seq_id;
	delete[] raw_sequence; 
	delete[] seq_id_repeat;
	delete[] tscore; 
}

// SE reads
void quantify_reads_rc(string fin, string fout, string compression_level) 
{
	gzFile output_file = gzopen(fout.c_str(), compression_level.c_str()); 
	if(output_file == NULL) 
	{
		printf("can't open %s file to write\n", fout.c_str()); 
		exit(EX_CANTCREAT); 
	}

	gzFile input_file = gzopen(fin.c_str(), "r"); 
	if(input_file == NULL) {
		printf("can't open %s file to read. \n", fin.c_str()); 
		exit(EX_NOINPUT); 
	}
	
	char * seq_id = new char[310];  
	char * raw_sequence = new char[310]; 
	char * seq_id_repeat = new char[310];
	char * tscore = new char[310];  

	unordered_set<double> vdat;

	double read_identifier = 0;
	double read_identifier_rc = 0;

	unsigned int line_num = 0; // number of reads
	unsigned int n_uniq_reads = 0; // number of unique reads
	while (gzgets(input_file, seq_id, 310) &&
		gzgets(input_file, raw_sequence, 310) &&
		gzgets(input_file, seq_id_repeat, 310) &&
		gzgets(input_file, tscore, 310)) {
		line_num++;
		double prod[2][2] = {{1,0},{0,1}};
		matrix_multiplication_helper(prod, raw_sequence, 'A');
		matrix_multiplication_helper(prod, raw_sequence, 'T');
		matrix_multiplication_helper(prod, raw_sequence, 'C');
		matrix_multiplication_helper(prod, raw_sequence, 'G');		
		read_identifier = prod[0][0] + sqrt(3) * prod[0][1] + M_SQRT2 * prod[1][0] + sqrt(5) * prod[1][1];		
		// check if the read has already appeared before
		// if unseen, keep it; else, skip it 
		if (vdat.insert(read_identifier).second == 1) {
			char * raw_sequence_rc = reverse_complement(raw_sequence);
			double prod_rc[2][2] = {{1,0},{0,1}};
			matrix_multiplication_helper(prod_rc, raw_sequence_rc, 'A');
			matrix_multiplication_helper(prod_rc, raw_sequence_rc, 'T');
			matrix_multiplication_helper(prod_rc, raw_sequence_rc, 'C');
			matrix_multiplication_helper(prod_rc, raw_sequence_rc, 'G');
			read_identifier_rc = prod_rc[0][0] + sqrt(3) * prod_rc[0][1] + M_SQRT2 * prod_rc[1][0] + sqrt(5) * prod_rc[1][1];
			vdat.insert(read_identifier_rc);
			gzprintf(output_file, "%s%s+\n%s", seq_id, raw_sequence, tscore);
			n_uniq_reads++;
			delete[] raw_sequence_rc;
		}
	}
	gzclose(input_file);
	gzclose(output_file); 
	cout << n_uniq_reads << "/" << line_num << " reads are unique." << endl; 

	vdat.clear();
	delete[] seq_id;
	delete[] raw_sequence; 
	delete[] seq_id_repeat;
	delete[] tscore; 
}

// SE reads, output removed reads
void quantify_reads_rc(string fin, string fout, string fremoved, string compression_level) 
{
	gzFile output_file = gzopen(fout.c_str(), compression_level.c_str()); 
	if(output_file == NULL) 
	{
		printf("can't open %s file to write\n", fout.c_str()); 
		exit(EX_CANTCREAT); 
	}

	gzFile removed_file = gzopen(fremoved.c_str(), compression_level.c_str());
	if (removed_file == NULL) {
		printf("can't open %s file to write\n", fremoved.c_str()); 
		exit(EX_CANTCREAT);
	}

	gzFile input_file = gzopen(fin.c_str(), "r"); 
	if(input_file == NULL) {
		printf("can't open %s file to read. \n", fin.c_str()); 
		exit(EX_NOINPUT); 
	}
	
	char * seq_id = new char[310];  
	char * raw_sequence = new char[310]; 
	char * seq_id_repeat = new char[310];
	char * tscore = new char[310];  

	unordered_set<double> vdat;

	double read_identifier = 0;
	double read_identifier_rc = 0;

	unsigned int line_num = 0; // number of reads
	unsigned int n_uniq_reads = 0; // number of unique reads
	while (gzgets(input_file, seq_id, 310) &&
		gzgets(input_file, raw_sequence, 310) &&
		gzgets(input_file, seq_id_repeat, 310) &&
		gzgets(input_file, tscore, 310)) {
		line_num++;
		double prod[2][2] = {{1,0},{0,1}};
		matrix_multiplication_helper(prod, raw_sequence, 'A');
		matrix_multiplication_helper(prod, raw_sequence, 'T');
		matrix_multiplication_helper(prod, raw_sequence, 'C');
		matrix_multiplication_helper(prod, raw_sequence, 'G');		
		read_identifier = prod[0][0] + sqrt(3) * prod[0][1] + M_SQRT2 * prod[1][0] + sqrt(5) * prod[1][1];		
		// check if the read has already appeared before
		// if unseen, keep it; else, skip it 
		if (vdat.insert(read_identifier).second == 1) {
			char * raw_sequence_rc = reverse_complement(raw_sequence);
			double prod_rc[2][2] = {{1,0},{0,1}};
			matrix_multiplication_helper(prod_rc, raw_sequence_rc, 'A');
			matrix_multiplication_helper(prod_rc, raw_sequence_rc, 'T');
			matrix_multiplication_helper(prod_rc, raw_sequence_rc, 'C');
			matrix_multiplication_helper(prod_rc, raw_sequence_rc, 'G');
			read_identifier_rc = prod_rc[0][0] + sqrt(3) * prod_rc[0][1] + M_SQRT2 * prod_rc[1][0] + sqrt(5) * prod_rc[1][1];
			vdat.insert(read_identifier_rc);
			gzprintf(output_file, "%s%s+\n%s", seq_id, raw_sequence, tscore);
			n_uniq_reads++;
			delete[] raw_sequence_rc;
		} else {
			gzprintf(removed_file, "%s%s+\n%s", seq_id, raw_sequence, tscore);
		}
	}
	gzclose(input_file);
	gzclose(removed_file);
	gzclose(output_file); 
	cout << n_uniq_reads << "/" << line_num << " reads are unique." << endl; 

	vdat.clear();
	delete[] seq_id;
	delete[] raw_sequence; 
	delete[] seq_id_repeat;
	delete[] tscore; 
}

// PE reads
void quantify_reads(string fin1, string fin2, string fout1, string fout2, string compression_level) 
{
	gzFile output_file1 = gzopen(fout1.c_str(), compression_level.c_str());	
	if(output_file1 == NULL) 
	{
		printf("can't open %s file to write\n", fout1.c_str()); 
		exit(EX_CANTCREAT); 
	}
	gzFile output_file2 = gzopen(fout2.c_str(), compression_level.c_str()); 
	if(output_file2 == NULL) 
	{
		printf("can't open %s file to write\n", fout2.c_str()); 
		exit(EX_CANTCREAT); 
	}

	gzFile input_file1 = gzopen(fin1.c_str(), "r"); 
	if(input_file1 == NULL) {
		printf("can't open %s file to read. \n", fin1.c_str()); 
		exit(EX_NOINPUT); 
	}
	gzFile input_file2 = gzopen(fin2.c_str(), "r"); 
	if(input_file2 == NULL) {
		printf("can't open %s file to read. \n", fin2.c_str()); 
		exit(EX_NOINPUT); 
	}
	
	char * seq_id1 = new char[310];  
	char * raw_sequence1 = new char[310]; 
	char * seq_id_repeat1 = new char[310];
	char * tscore1 = new char[310];  

	char * seq_id2 = new char[310];  
	char * raw_sequence2 = new char[310]; 
	char * seq_id_repeat2 = new char[310];
	char * tscore2 = new char[310];  

	unordered_multimap<double, unsigned int> vdat1;
	unordered_map<unsigned int, double> vdat2;

	double read_identifier1 = 0;
	double read_identifier2 = 0;

	unsigned int line_num = 0; // number of reads
	unsigned int n_uniq_reads = 0; // number of unique reads
	while (gzgets(input_file1, seq_id1, 310) &&
		gzgets(input_file1, raw_sequence1, 310) &&
		gzgets(input_file1, seq_id_repeat1, 310) &&
		gzgets(input_file1, tscore1, 310) &&
		gzgets(input_file2, seq_id2, 310) &&
		gzgets(input_file2, raw_sequence2, 310) &&
		gzgets(input_file2, seq_id_repeat2, 310) &&
		gzgets(input_file2, tscore2, 310)) {
		line_num++;
		double prod1[2][2] = {{1,0},{0,1}};
		double prod2[2][2] = {{1,0},{0,1}};
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
		if (pos1.first == vdat1.end()) {
			vdat1.insert({read_identifier1, line_num});
			vdat2.insert({line_num, read_identifier2});
			gzprintf(output_file1, "%s%s+\n%s", seq_id1, raw_sequence1, tscore1);
			gzprintf(output_file2, "%s%s+\n%s", seq_id2, raw_sequence2, tscore2);
			n_uniq_reads++;
		} else {
			bool seen = false;
			while (!seen && pos1.first != pos1.second) {
				if (vdat2[pos1.first->second] == read_identifier2) {
					seen = true;
				}
				pos1.first++;
			}
			if (!seen) {
				vdat1.insert(pos1.second, pair<double, unsigned int>(read_identifier1, line_num)); // gives it hint
				vdat2.insert({line_num, read_identifier2});
				gzprintf(output_file1, "%s%s+\n%s", seq_id1, raw_sequence1, tscore1);
				gzprintf(output_file2, "%s%s+\n%s", seq_id2, raw_sequence2, tscore2);
				n_uniq_reads++;
			}
		}
	}
	gzclose(input_file1);
	gzclose(input_file2);
	gzclose(output_file1);
	gzclose(output_file2); 
	cout << n_uniq_reads << "/" << line_num << " read pairs are unique." << endl; 

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

// PE reads, output removed reads
void quantify_reads(string fin1, string fin2, string fout1, string fout2, string fremoved1, string fremoved2, string compression_level) 
{
	gzFile output_file1 = gzopen(fout1.c_str(), compression_level.c_str());	
	if(output_file1 == NULL) 
	{
		printf("can't open %s file to write\n", fout1.c_str()); 
		exit(EX_CANTCREAT); 
	}
	gzFile output_file2 = gzopen(fout2.c_str(), compression_level.c_str()); 
	if(output_file2 == NULL) 
	{
		printf("can't open %s file to write\n", fout2.c_str()); 
		exit(EX_CANTCREAT); 
	}

	gzFile removed_file1 = gzopen(fremoved1.c_str(), compression_level.c_str());
	if (removed_file1 == NULL) {
		printf("can't open %s file to write\n", fremoved1.c_str()); 
		exit(EX_CANTCREAT);
	}
	gzFile removed_file2 = gzopen(fremoved2.c_str(), compression_level.c_str());
	if (removed_file2 == NULL) {
		printf("can't open %s file to write\n", fremoved2.c_str()); 
		exit(EX_CANTCREAT);
	}

	gzFile input_file1 = gzopen(fin1.c_str(), "r"); 
	if(input_file1 == NULL) {
		printf("can't open %s file to read. \n", fin1.c_str()); 
		exit(EX_NOINPUT); 
	}
	gzFile input_file2 = gzopen(fin2.c_str(), "r"); 
	if(input_file2 == NULL) {
		printf("can't open %s file to read. \n", fin2.c_str()); 
		exit(EX_NOINPUT); 
	}
	
	char * seq_id1 = new char[310];  
	char * raw_sequence1 = new char[310]; 
	char * seq_id_repeat1 = new char[310];
	char * tscore1 = new char[310];  

	char * seq_id2 = new char[310];  
	char * raw_sequence2 = new char[310]; 
	char * seq_id_repeat2 = new char[310];
	char * tscore2 = new char[310];  

	unordered_multimap<double, unsigned int> vdat1;
	unordered_map<unsigned int, double> vdat2;

	double read_identifier1 = 0;
	double read_identifier2 = 0;

	unsigned int line_num = 0; // number of reads
	unsigned int n_uniq_reads = 0; // number of unique reads
	while (gzgets(input_file1, seq_id1, 310) &&
		gzgets(input_file1, raw_sequence1, 310) &&
		gzgets(input_file1, seq_id_repeat1, 310) &&
		gzgets(input_file1, tscore1, 310) &&
		gzgets(input_file2, seq_id2, 310) &&
		gzgets(input_file2, raw_sequence2, 310) &&
		gzgets(input_file2, seq_id_repeat2, 310) &&
		gzgets(input_file2, tscore2, 310)) {
		line_num++;
		double prod1[2][2] = {{1,0},{0,1}};
		double prod2[2][2] = {{1,0},{0,1}};
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
		if (pos1.first == vdat1.end()) {
			vdat1.insert({read_identifier1, line_num});
			vdat2.insert({line_num, read_identifier2});
			gzprintf(output_file1, "%s%s+\n%s", seq_id1, raw_sequence1, tscore1);
			gzprintf(output_file2, "%s%s+\n%s", seq_id2, raw_sequence2, tscore2);
			n_uniq_reads++;
		} else {
			bool seen = false;
			while (!seen && pos1.first != pos1.second) {
				if (vdat2[pos1.first->second] == read_identifier2) {
					seen = true;
				}
				pos1.first++;
			}
			if (!seen) {
				vdat1.insert(pos1.second, pair<double, unsigned int>(read_identifier1, line_num)); // gives it hint
				vdat2.insert({line_num, read_identifier2});
				gzprintf(output_file1, "%s%s+\n%s", seq_id1, raw_sequence1, tscore1);
				gzprintf(output_file2, "%s%s+\n%s", seq_id2, raw_sequence2, tscore2);
				n_uniq_reads++;
			} else {
				gzprintf(removed_file1, "%s%s+\n%s", seq_id1, raw_sequence1, tscore1);
				gzprintf(removed_file2, "%s%s+\n%s", seq_id2, raw_sequence2, tscore2);
			}
		}
	}
	gzclose(input_file1);
	gzclose(input_file2);
	gzclose(output_file1);
	gzclose(output_file2); 
	gzclose(removed_file1);
	gzclose(removed_file2);
	cout << n_uniq_reads << "/" << line_num << " read pairs are unique." << endl; 

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
void quantify_reads_rc(string fin1, string fin2, string fout1, string fout2, string compression_level) 
{
	gzFile output_file1 = gzopen(fout1.c_str(), compression_level.c_str());	
	if(output_file1 == NULL) 
	{
		printf("can't open %s file to write\n", fout1.c_str()); 
		exit(EX_CANTCREAT); 
	}
	gzFile output_file2 = gzopen(fout2.c_str(), compression_level.c_str()); 
	if(output_file2 == NULL) 
	{
		printf("can't open %s file to write\n", fout2.c_str()); 
		exit(EX_CANTCREAT); 
	}

	gzFile input_file1 = gzopen(fin1.c_str(), "r"); 
	if(input_file1 == NULL) {
		printf("can't open %s file to read. \n", fin1.c_str()); 
		exit(EX_NOINPUT); 
	}
	gzFile input_file2 = gzopen(fin2.c_str(), "r"); 
	if(input_file2 == NULL) {
		printf("can't open %s file to read. \n", fin2.c_str()); 
		exit(EX_NOINPUT); 
	}
	
	char * seq_id1 = new char[310];  
	char * raw_sequence1 = new char[310]; 
	char * seq_id_repeat1 = new char[310];
	char * tscore1 = new char[310];  

	char * seq_id2 = new char[310];  
	char * raw_sequence2 = new char[310]; 
	char * seq_id_repeat2 = new char[310];
	char * tscore2 = new char[310];  

	unordered_multimap<double, unsigned int> vdat1;
	unordered_map<unsigned int, double> vdat2;

	double read_identifier1 = 0;
	double read_identifier2 = 0;

	unsigned int line_num = 0; // number of reads
	unsigned int n_uniq_reads = 0; // number of unique reads
	while (gzgets(input_file1, seq_id1, 310) &&
		gzgets(input_file1, raw_sequence1, 310) &&
		gzgets(input_file1, seq_id_repeat1, 310) &&
		gzgets(input_file1, tscore1, 310) &&
		gzgets(input_file2, seq_id2, 310) &&
		gzgets(input_file2, raw_sequence2, 310) &&
		gzgets(input_file2, seq_id_repeat2, 310) &&
		gzgets(input_file2, tscore2, 310)) {
		line_num++;
		double prod1[2][2] = {{1,0},{0,1}};
		double prod2[2][2] = {{1,0},{0,1}};
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
		bool seen = false;
		if (pos1.first != vdat1.end()) {
			while (!seen && pos1.first != pos1.second) {
				if (vdat2[pos1.first->second] == read_identifier2) {
					seen = true;
				}
				pos1.first++;
			}
		}
		// seen == 1 denotes duplicate
		// seen == 0 denotes map1 doesn't contain r1 
		// or map1 contains r1 but map2 no r2 in corresponding lines
		if (!seen) {
			auto pos2 = vdat1.equal_range(read_identifier2);
			if (pos2.first == vdat1.end()) { // r2 not in map1
				vdat1.insert(pos1.second, pair<double, unsigned int>(read_identifier1, line_num)); // gives it hint
				vdat2.insert({line_num, read_identifier2});
				gzprintf(output_file1, "%s%s+\n%s", seq_id1, raw_sequence1, tscore1);
				gzprintf(output_file2, "%s%s+\n%s", seq_id2, raw_sequence2, tscore2);
				n_uniq_reads++;
			} else { // r2 in map1
				while (!seen && pos2.first != pos2.second) {
					if (vdat2[pos2.first->second] == read_identifier1) {
						seen = true;
					}
					pos2.first++;
				}
				if (!seen) {
					vdat1.insert(pos1.second, pair<double, unsigned int>(read_identifier1, line_num));
					vdat2.insert({line_num, read_identifier2});
					gzprintf(output_file1, "%s%s+\n%s", seq_id1, raw_sequence1, tscore1);
					gzprintf(output_file2, "%s%s+\n%s", seq_id2, raw_sequence2, tscore2);
					n_uniq_reads++;
				}
			}
		}
	}
	gzclose(input_file1);
	gzclose(input_file2);
	gzclose(output_file1);
	gzclose(output_file2); 
	cout << n_uniq_reads << "/" << line_num << " read pairs are unique." << endl; 

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

// PE reads, output removed reads
void quantify_reads_rc(string fin1, string fin2, string fout1, string fout2, string fremoved1, string fremoved2, string compression_level) 
{
	gzFile output_file1 = gzopen(fout1.c_str(), compression_level.c_str());	
	if(output_file1 == NULL) 
	{
		printf("can't open %s file to write\n", fout1.c_str()); 
		exit(EX_CANTCREAT); 
	}
	gzFile output_file2 = gzopen(fout2.c_str(), compression_level.c_str()); 
	if(output_file2 == NULL) 
	{
		printf("can't open %s file to write\n", fout2.c_str()); 
		exit(EX_CANTCREAT); 
	}

	gzFile removed_file1 = gzopen(fremoved1.c_str(), compression_level.c_str());
	if (removed_file1 == NULL) {
		printf("can't open %s file to write\n", fremoved1.c_str()); 
		exit(EX_CANTCREAT);
	}
	gzFile removed_file2 = gzopen(fremoved2.c_str(), compression_level.c_str());
	if (removed_file2 == NULL) {
		printf("can't open %s file to write\n", fremoved2.c_str()); 
		exit(EX_CANTCREAT);
	}

	gzFile input_file1 = gzopen(fin1.c_str(), "r"); 
	if(input_file1 == NULL) {
		printf("can't open %s file to read. \n", fin1.c_str()); 
		exit(EX_NOINPUT); 
	}
	gzFile input_file2 = gzopen(fin2.c_str(), "r"); 
	if(input_file2 == NULL) {
		printf("can't open %s file to read. \n", fin2.c_str()); 
		exit(EX_NOINPUT); 
	}
	
	char * seq_id1 = new char[310];  
	char * raw_sequence1 = new char[310]; 
	char * seq_id_repeat1 = new char[310];
	char * tscore1 = new char[310];  

	char * seq_id2 = new char[310];  
	char * raw_sequence2 = new char[310]; 
	char * seq_id_repeat2 = new char[310];
	char * tscore2 = new char[310];  

	unordered_multimap<double, unsigned int> vdat1;
	unordered_map<unsigned int, double> vdat2;

	double read_identifier1 = 0;
	double read_identifier2 = 0;

	unsigned int line_num = 0; // number of reads
	unsigned int n_uniq_reads = 0; // number of unique reads
	while (gzgets(input_file1, seq_id1, 310) &&
		gzgets(input_file1, raw_sequence1, 310) &&
		gzgets(input_file1, seq_id_repeat1, 310) &&
		gzgets(input_file1, tscore1, 310) &&
		gzgets(input_file2, seq_id2, 310) &&
		gzgets(input_file2, raw_sequence2, 310) &&
		gzgets(input_file2, seq_id_repeat2, 310) &&
		gzgets(input_file2, tscore2, 310)) {
		line_num++;
		double prod1[2][2] = {{1,0},{0,1}};
		double prod2[2][2] = {{1,0},{0,1}};
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
		bool seen = false;
		if (pos1.first != vdat1.end()) {
			while (!seen && pos1.first != pos1.second) {
				if (vdat2[pos1.first->second] == read_identifier2) {
					seen = true;
				}
				pos1.first++;
			}
		}
		// seen == 1 denotes duplicate
		// seen == 0 denotes map1 doesn't contain r1 
		// or map1 contains r1 but map2 no r2 in corresponding lines
		if (!seen) {
			auto pos2 = vdat1.equal_range(read_identifier2);
			if (pos2.first == vdat1.end()) { // r2 not in map1
				vdat1.insert(pos1.second, pair<double, unsigned int>(read_identifier1, line_num)); // gives it hint
				vdat2.insert({line_num, read_identifier2});
				gzprintf(output_file1, "%s%s+\n%s", seq_id1, raw_sequence1, tscore1);
				gzprintf(output_file2, "%s%s+\n%s", seq_id2, raw_sequence2, tscore2);
				n_uniq_reads++;
			} else { // r2 in map1
				while (!seen && pos2.first != pos2.second) {
					if (vdat2[pos2.first->second] == read_identifier1) {
						seen = true;
					}
					pos2.first++;
				}
				if (!seen) {
					vdat1.insert(pos1.second, pair<double, unsigned int>(read_identifier1, line_num));
					vdat2.insert({line_num, read_identifier2});
					gzprintf(output_file1, "%s%s+\n%s", seq_id1, raw_sequence1, tscore1);
					gzprintf(output_file2, "%s%s+\n%s", seq_id2, raw_sequence2, tscore2);
					n_uniq_reads++;
				} else {
					gzprintf(removed_file1, "%s%s+\n%s", seq_id1, raw_sequence1, tscore1);
					gzprintf(removed_file2, "%s%s+\n%s", seq_id2, raw_sequence2, tscore2);
				}
			}
		} else {
			gzprintf(removed_file1, "%s%s+\n%s", seq_id1, raw_sequence1, tscore1);
			gzprintf(removed_file2, "%s%s+\n%s", seq_id2, raw_sequence2, tscore2);
		}
	}
	gzclose(input_file1);
	gzclose(input_file2);
	gzclose(output_file1);
	gzclose(output_file2); 
	gzclose(removed_file1);
	gzclose(removed_file2);
	cout << n_uniq_reads << "/" << line_num << " read pairs are unique." << endl; 

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

void operation(double prod[][2], bool yes) 
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

void matrix_multiplication_helper(double prod[][2], char * raw_seq, char nucleotide) {
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
