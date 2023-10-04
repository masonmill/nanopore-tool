/*
 *      Test input:
 *     -------------
 *          -o Nanopal-out
 *          -b /home/masonmil/MEI_AD/data/MEI_AD_flongle_opt.workspace.LINE/Nanopore.sorted.bam
 *          -t 10
 *          -m /home/masonmil/MEI_AD/lib/L1.3
 */

#include <cstring>
#include <filesystem> // as of C++17
#include <fstream>
#include <getopt.h>
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {
    string FC;
    string BAM;
    int nprocs;
    string MEI;

    // process command line
    opterr = false;
    int opt, index = 0;
    static struct option long_options[] = {
        { "FC",     required_argument, nullptr, 'o'  },
        { "BAM",    required_argument, nullptr, 'b'  },
        { "nprocs", required_argument, nullptr, 't'  },
        { "MEI",    required_argument, nullptr, 'm'  },
        { nullptr,  0,                 nullptr, '\0' }
    };
    while((opt = getopt_long(argc, argv, "o:b:t:m:", long_options, &index)) != -1) {
        switch (opt) {
            case 'o':
                FC = optarg;
                break;
            case 'b':
                BAM = optarg;
                break;
            case 't':
                nprocs = atoi(optarg);
                break;
            case 'm':
                MEI = optarg;
                break;
            default:
                cerr << "Invalid option: -" << static_cast<char>(opt) << endl;
                exit(1);
        }
    }

    cout << "Nano-Pal Input\n" 
         << "BAM: " << BAM << "\n" 
         << "MEI: " << MEI << "\n"
         << "Directory: " << FC << "\n";

    // Other input params
    string lib = "/home/masonmil/software/NanoPal/NanoPal-and-Cas9-targeted-enrichment-pipelines-1.0/lib";
    string REF = "/data/genomes/hg38/seq/hg38.fa";
    string MINIMAP = "/home/masonmil/software/minimap2-2.17_x64-linux/minimap2";

    // MEI reference sequence
    string L1Hs = lib + "/L1.3";
    string S_refLINE = lib + "/hg38.RM.L1.ref";

    // PALMER call for NA12878
    string S_origLINE = lib + "/PALMER.NA12878.L1.txt";

    // P&P call for NA12878
    string S_PP_LINE = lib + "/union.1214/L1.inter.fi";

    // constuct workspace
    string dir = FC + ".workspace";
    /*
     *  Note: workspace directory name modified from ${FC}.workspace.${MEI} to ${FC}.workspace
     *      due to errors with directories and file paths 
     */
    string mkdir = "mkdir " + dir;
    system(mkdir.c_str());
    filesystem::current_path(dir);

    // RAW fastq data & generate fasta
    string fastq = "samtools fastq " + BAM + " > batch.fastq";
    system(fastq.c_str());
    system("cat batch.fastq | awk '{if(NR%4==1) {printf(\">%s\\n\",substr($0,2));} else if(NR%4==2) print;}' > batch.fasta");

    // alignment
    string align = MINIMAP + " -ax map-ont --secondary=no --eqx -Y -t " + to_string(nprocs) + " " + REF + " batch.fastq " + 
                   " | samtools sort -o Nanopore.sorted.bam";
    system(align.c_str());
    system("samtools index Nanopore.sorted.bam");

    system("mkdir palmer");

    // PALMER
    string chr_list = "/home/masonmil/software/NanoPal/NanoPal-and-Cas9-targeted-enrichment-pipelines-1.0/lib/chr.list";
    string line, do_palmer;
    ifstream list_in;
    list_in.open(chr_list);
    if (list_in.fail()) {
        cerr << "could not open ${lib}/chr.list" << endl;
        exit(1);
    }

    while (getline(list_in, line)) {
        do_palmer = "/home/masonmil/software/PALMER/PALMER --input Nanopore.sorted.bam --mode raw --workdir palmer/ --output " 
                    + line + " --ref_ver GRCh38 --ref_fa '/data/genomes/hg38/seq/hg38.fa' --type LINE --chr " + line;
        system(do_palmer.c_str());
    }
    list_in.close();

    system("rm blastn_refine.all.txt");
    system("cat palmer/chr*_*_*/blastn_refine.txt >> blastn_refine.all.txt");

    // don't need this dir again
    system("rm -r palmer");

    // Tool: Cigar information parsing
    string cp = "cp " + lib + "/cigar.parser.1106.o .";
    system(cp.c_str());

    system("samtools view Nanopore.sorted.bam -q 10 -F 0x100 -F 0x200 -F 0x800 -F 0x400 | awk '{print $1,$3,$4,$6}' > mapped.info.txt");

    ifstream mappedinfoin("mapped.info.txt");
    string a, b, c, d;
    while (mappedinfoin >> a >> b >> c >> d) {
        ofstream outfile("input_cigar");
        outfile << d;
        outfile.close();
        system("./cigar.parser.1106.o");
        ifstream cigar_results("cigar_results");
        ofstream cigar_results_all("cigar_results.all.txt", ios::app);
        cigar_results_all << cigar_results.rdbuf();
        cigar_results.close();
        cigar_results_all.close();
    }

    // Get the information of start and end for each read with L1Hs signal
    system("awk '{print $4+$6+$10}' cigar_results.all.txt > cigar.ref.txt");
    system("paste mapped.info.txt cigar.ref.txt | awk '{print $1,$2,$3,$3+$5}' >  mapped.info.final.txt");
    
    // Given BAM file, convert to FASTA
    // Resulting FASTA file will be saved in current directory
    // FASTA file: all_reads.fa
    filesystem::current_path("/home/masonmil/MEI_AD/cblastn");
    string bamtofasta = "samtools fasta " + BAM + " > all_reads.fasta"; // input must have file ext .fasta for cblastn
    system(bamtofasta.c_str());

    // Use resulting file to create a BLAST database
    // Resulting BLAST database will be saved in current directory
    // Database files: all_reads.nhr, all_reads.nin, all_reads.nsq
    system("makeblastdb -in all_reads.fasta -title all_reads -dbtype nucl -out all_reads");

    // Create input to CBlastn::do_blastn
    // The all_reads database files and the MEI file must be accessible by cblastn (same directory or full path)
    // CBlastn output file: read.all.txt
    string do_blastn = "/home/masonmil/MEI_AD/cblastn/cblastn -db all_reads -in " + MEI + 
                       " -out /home/masonmil/MEI_AD/" + dir + "/read.all.txt -evalue 0.001";
    system(do_blastn.c_str()); 

    /*
     *      Note: cblastn output will require reformatting before continuing with the translation of 
     *          the Nano-Pal_CM.LINE.parallel script to C/C++.
     * 
     */

    return 0;
}
