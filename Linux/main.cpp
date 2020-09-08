#include "functions.h"

int main(int argc, char ** argv)
{
	string fin; 
	string fout; 
	string fin1; 
	string fin2; 
	string fout1; 
	string fout2;
	string fremoved; 
	string fremoved1; 
	string fremoved2; 
	bool strand = true; // whether consider reads from complementary strand or not
	bool write_remove = false; // whether output removed reads or not
	int gz = 0; // indicating compression level when writing gzipped file, 0-9
	string compression_level = "wT"; // default, change if gz > 0
	string file_name_suffix = ".uniq.fastq"; // default, change if gz > 0
	string removed_name_suffix = ".removed.fastq"; // default, change if gz > 0
	bool qtf = false; // let it run or not 

	for (int i = 1; i < argc; i++) {
		string str;
		
		if (i>1 && argv[i][0] != '-') 
			continue;
		str.assign(argv[i]);
		if (str.compare("-h") == 0) {
			cout << "./nubeam-dedup [-i -o -d    -i1 -i2 -o1 -o2 -d1 -d2    -s -r -z -h]\n\n"; 
			cout << "Remove exact PCR duplicates for sequencing reads in (gzipped) fastq format.\n";
			cout << "Produces de-duplicated reads in fastq files.\n\n"; 

			cout << "Parameters for single-end (SE) reads:\n";
			cout << "--in or -i: Input file name. The parameter is mandatory for SE reads.\n"; 
			cout << "--out or -o: Output file name for unique reads. The default is input file name prefix appended with '.uniq.fastq(.gz)', under the current directory.\n"; 
			cout << "--duplicate or -d: File name for removed duplicated reads. The parameter is only valid when --remove or -r is set as 1 (see below). The default is input file name prefix appended with '.removed.fastq(.gz)', under the current directory.\n\n"; 
			
			cout << "Parameters for paired-end (PE) reads:\n";
			cout << "--in1 or -i1: Input file name for read 1 file. The parameter is mandatory for PE reads.\n"; 
			cout << "--in2 or -i2: Input file name for read 2 file. The parameter is mandatory for PE reads.\n"; 
			cout << "--out1 or -o1: Output file name for unique read pairs read 1 file. The default is read 1 file name prefix appended with '.uniq.fastq(.gz)', under the current working directory.\n";
			cout << "--out2 or -o2: Output file name for unique read pairs read 2 file. The default is read 2 file name prefix appended with '.uniq.fastq(.gz)', under the current working directory.\n";
			cout << "--duplicate1 or -d1: File name for removed duplicated read pairs read 1 file. The parameter is only valid when --remove or -r is set as 1 (see below). The default is input file name prefix appended with '.removed.fastq(.gz)', under the current directory.\n"; 
			cout << "--duplicate2 or -d2: File name for removed duplicated read pairs read 2 file. The parameter is only valid when --remove or -r is set as 1 (see below). The default is input file name prefix appended with '.removed.fastq(.gz)', under the current directory.\n\n"; 
			
			cout << "Miscellaneous parameters:\n";
			cout << "--strand or -s: Whether take reads from complementary strand into account. Accept boolean 1 (default) or 0.\n";
			cout << "--remove or -r: Whether output removed duplicated reads. Accept boolean 0 (default) or 1.\n";
			cout << "--gz or -z: Compression level of output file. Accept integer 0 (default) to 9. " << 
			"If 0 (default), the output data will not be compressed and will be written to plain text file; " << 
			"otherwise, the output data will be written to gzip format file, with the compression level suggested by user. " << 
			"If compression is needed, a compression level less than 3 is recommended as a compromise between speed and compression.\n\n";
			cout << "-h: print this help\n"; 
			exit(0); 
		}
		else if (str.compare("--in") == 0 || str.compare("-i") == 0) {
			if (argv[i+1] == NULL || argv[i+1][0] == '-') continue;
			fin.assign(argv[i+1]);
			qtf = true;
		}
		else if (str.compare("--out") == 0 || str.compare("-o") == 0) {
			if (argv[i+1] == NULL || argv[i+1][0] == '-') continue;
			fout.assign(argv[i+1]);
		}
		else if (str.compare("--duplicate") == 0 || str.compare("-d") == 0) {
			if (argv[i+1] == NULL || argv[i+1][0] == '-') continue;
			fremoved.assign(argv[i+1]);
		}
		else if (str.compare("--in1") == 0 || str.compare("-i1") == 0) {
			if (argv[i+1] == NULL || argv[i+1][0] == '-') continue;
			fin1.assign(argv[i+1]);
			qtf = true;
		}
		else if (str.compare("--in2") == 0 || str.compare("-i2") == 0) {
			if (argv[i+1] == NULL || argv[i+1][0] == '-') continue;
			fin2.assign(argv[i+1]);
		}
		else if (str.compare("--out1") == 0 || str.compare("-o1") == 0) {
			if (argv[i+1] == NULL || argv[i+1][0] == '-') continue;
			fout1.assign(argv[i+1]);
		}
		else if (str.compare("--out2") == 0 || str.compare("-o2") == 0) {
			if (argv[i+1] == NULL || argv[i+1][0] == '-') continue;
			fout2.assign(argv[i+1]);
		}
		else if (str.compare("--duplicate1") == 0 || str.compare("-d1") == 0) {
			if (argv[i+1] == NULL || argv[i+1][0] == '-') continue;
			fremoved1.assign(argv[i+1]);
		}
		else if (str.compare("--duplicate2") == 0 || str.compare("-d2") == 0) {
			if (argv[i+1] == NULL || argv[i+1][0] == '-') continue;
			fremoved2.assign(argv[i+1]);
		}
		else if (str.compare("--strand") == 0 || str.compare("-s") == 0) {
			if (argv[i+1] == NULL || argv[i+1][0] == '-') continue;
			if (!isdigit(argv[i+1][0])) {
				printf("wrong augument after option --strand.\n");
				exit(EX_USAGE);
			}
			strand = atoi(argv[i+1]);
		}
		else if (str.compare("--remove") == 0 || str.compare("-r") == 0) {
			if (argv[i+1] == NULL || argv[i+1][0] == '-') continue;
			if (!isdigit(argv[i+1][0])) {
				printf("wrong augument after option --remove.\n");
				exit(EX_USAGE);
			}
			write_remove = atoi(argv[i+1]);
		}
		else if (str.compare("--gz") == 0 || str.compare("-z") == 0) {
			if (argv[i+1] == NULL || argv[i+1][0] == '-') continue;
			if (!isdigit(argv[i+1][0])) {
				printf("wrong augument after option --gz.\n");
				exit(EX_USAGE);
			}
			gz = atoi(argv[i+1]);
		}
		else {
			fprintf(stderr,"Bad option %s\n", argv[i]);
			exit(EX_USAGE);
		}
	}
	
	if (qtf) {
		if (gz > 0) { // determine compression level
			compression_level = "wb" + to_string(gz);
			file_name_suffix = ".uniq.fastq.gz";
			removed_name_suffix = ".removed.fastq.gz";
		}
		if (fin1.empty()) { // SE reads
			if (fout.empty() || fremoved.empty()) { // means we need to figure out the file name by ourselves
				string fin_copy(fin); // copy fin
				auto search_path = fin_copy.rfind("/");
				if (search_path != string::npos) { // if find '/'
					fin_copy.erase(0, search_path + 1);
				}
				auto search_suffix = fin_copy.rfind(".fastq");
				if (search_suffix == string::npos) { // if fail to find "fastq", go to find "fq"
					search_suffix = fin_copy.rfind(".fq");
				}
				if (search_suffix != string::npos) { // if can find "fastq" or "fq"
					fin_copy.erase(search_suffix);
				}

				char buff[FILENAME_MAX];
				getcwd(buff, FILENAME_MAX);
				string current_working_dir(buff); // copy buff
				if (fout.empty()) { // which is default, else use what suggested by user
					fout = current_working_dir + "/" + fin_copy + file_name_suffix;
				}
				if (write_remove) {
					if (fremoved.empty()) { // which is default, else use what suggested by user
						fremoved = current_working_dir + "/" + fin_copy + removed_name_suffix;
					}
					clog << "Output removed duplicated reads to " << fremoved << endl;
				}
			}
			clog << "Output unique reads to " << fout << endl;

			if (strand) {
				if (write_remove) {
					quantify_reads_rc(fin, fout, fremoved, compression_level);
				} else {
					quantify_reads_rc(fin, fout, compression_level);
				}
			} else {
				if (write_remove) {
					quantify_reads(fin, fout, fremoved, compression_level);
				} else {
					quantify_reads(fin, fout, compression_level);
				}
			}
		} else { // PE reads
			if (fout1.empty() || fremoved1.empty()) { // means we need to figure out the file name by ourselves
				string fin_copy1(fin1); // copy fin1
				auto search_path1 = fin_copy1.rfind("/");
				if (search_path1 != string::npos) { // if find '/'
					fin_copy1.erase(0, search_path1 + 1);
				}
				auto search_suffix1 = fin_copy1.rfind(".fastq");
				if (search_suffix1 == string::npos) { // if fail to find "fastq", go to find "fq"
					search_suffix1 = fin_copy1.rfind(".fq");
				}
				if (search_suffix1 != string::npos) { // if can find "fastq" or "fq"
					fin_copy1.erase(search_suffix1);
				}
				string fin_copy2(fin2); // copy fin2
				auto search_path2 = fin_copy2.rfind("/");
				if (search_path2 != string::npos) { // if find '/'
					fin_copy2.erase(0, search_path2 + 1);
				}
				auto search_suffix2 = fin_copy2.rfind(".fastq");
				if (search_suffix2 == string::npos) { // if fail to find "fastq", go to find "fq"
					search_suffix2 = fin_copy2.rfind(".fq");
				}
				if (search_suffix2 != string::npos) { // if can find "fastq" or "fq"
					fin_copy2.erase(search_suffix2);
				}
				char buff[FILENAME_MAX];
				getcwd(buff, FILENAME_MAX);
				string current_working_dir(buff); // copy buff
				if (fout1.empty()) { // which is default, else use what suggested by user
					fout1 = current_working_dir + "/" + fin_copy1 + file_name_suffix;
					fout2 = current_working_dir + "/" + fin_copy2 + file_name_suffix;
				}
				if (write_remove) {
					if (fremoved1.empty()) { // which is default, else use what suggested by user
						fremoved1 = current_working_dir + "/" + fin_copy1 + removed_name_suffix;
						fremoved2 = current_working_dir + "/" + fin_copy2 + removed_name_suffix;
					}
					clog << "Output removed duplicated read pairs read 1 to " << fremoved1 << "\nOutput removed duplicated read pairs read 2 to " << fremoved2 << endl;
				}
			}
			clog << "Output unique read pairs read 1 to " << fout1 << "\nOutput unique read pairs read 2 to " << fout2 << endl;
			
			if (strand) {
				if (write_remove) {
					quantify_reads_rc(fin1, fin2, fout1, fout2, fremoved1, fremoved2, compression_level);
				} else {
					quantify_reads_rc(fin1, fin2, fout1, fout2, compression_level); 
				}
			} else {
				if (write_remove) {
					quantify_reads(fin1, fin2, fout1, fout2, fremoved1, fremoved2, compression_level);
				} else {
					quantify_reads(fin1, fin2, fout1, fout2, compression_level); 
				}
			}
		}
		return 0; 
	} else {
		cerr << "No input file(s). Use ./nubeam-dedup -h for more options.\n"; 
		exit(EX_USAGE);
		return EX_USAGE;
	}
	return 0; 
}
