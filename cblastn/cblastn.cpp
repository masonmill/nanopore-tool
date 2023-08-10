/*
    Input:
   --------
        -db         Name of BLAST database from makeblastdb (no file ext)
        -in         Query FASTA file (formatted as <read_name> \n <sequence> \n ...)
        -out        Output file
        -evalue     E-value threshold for saving hits

    Output:
   ---------
        Summary of CBlastn results

*/

#include <cstring>
#include <fstream>

#include <ncbi_pch.hpp>
#include <corelib/ncbiapp.hpp>
#include <corelib/ncbienv.hpp>
#include <corelib/ncbiargs.hpp>

#include <objmgr/object_manager.hpp>

#include <objects/seqalign/Seq_align_set.hpp>

#include <algo/blast/api/sseqloc.hpp>
#include <algo/blast/api/local_blast.hpp>
#include <algo/blast/api/uniform_search.hpp>
#include <algo/blast/api/blast_types.hpp>
#include <algo/blast/api/blast_aux.hpp>
#include <algo/blast/api/objmgr_query_data.hpp>
#include <algo/blast/api/blast_options_handle.hpp>
#include <algo/blast/api/blast_nucl_options.hpp>

#include <algo/blast/blastinput/blast_input.hpp>
#include <algo/blast/blastinput/blast_fasta_input.hpp>

USING_NCBI_SCOPE;
USING_SCOPE(blast);

/*
 * Boilerplate command-line argument processing code automatically generated
 * by new_project.sh script in NCBI C++ Toolkit
 *
 *     Authors:  Denis Vakatov, Vladimir Ivanov
 */

/////////////////////////////////////////////////////////////////////////////
//  CBlastn::

class CBlastn : public CNcbiApplication
{
private:
    virtual void Init(void);
    virtual int  Run(void);
    virtual void Exit(void);
    CSearchResultSet full_blast(TSeqLocVector query_loc, CRef<CBlastOptionsHandle> opts, Uint8 *full_db_size);
    void ProcessCommandLineArgs(CRef<CBlastOptionsHandle> opts_handle);
};


/////////////////////////////////////////////////////////////////////////////
//  Init test for all different types of arguments

void CBlastn::Init(void)
{
    // Create command-line argument descriptions class
    unique_ptr<CArgDescriptions> arg_desc(new CArgDescriptions);

    // Specify USAGE context
    arg_desc->SetUsageContext(GetArguments().GetProgramBasename(), "C++ blastn program");

    arg_desc->AddKey
        ("db", "Database",
         "This is the name of the database",
         CArgDescriptions::eString);

    arg_desc->AddKey("in", "QueryFile",
                        "A file with the query", CArgDescriptions::eInputFile);

    arg_desc->AddKey("out", "OutputFile",
                        "The output file", CArgDescriptions::eOutputFile);

    arg_desc->AddKey("evalue", "evalue",
                        "E-value threshold for final hits", CArgDescriptions::eDouble);

    // Setup arg.descriptions for this application
    SetupArgDescriptions(arg_desc.release());
}


/// Modify BLAST options from defaults based upon command-line args.
///
/// @param opts_handle already created CBlastOptionsHandle to modify [in]

void CBlastn::ProcessCommandLineArgs(CRef<CBlastOptionsHandle> opts_handle)
{
	const CArgs& args = GetArgs();

    // Expect value is a supported option for all flavors of BLAST.
    if(args["evalue"].AsDouble())
        opts_handle->SetEvalueThreshold(args["evalue"].AsDouble());

    return;
}

CSearchResultSet CBlastn::full_blast(TSeqLocVector query_loc, CRef<CBlastOptionsHandle> opts, Uint8 *db_size) 
{
    const CArgs& args = GetArgs();

    const CSearchDatabase target_db(args["db"].AsString(), CSearchDatabase::eBlastDbIsNucleotide);

    *db_size = target_db.GetSeqDb()->GetTotalLength();
    CNcbiOstream& out = args["out"].AsOutputFile();
    out << "Database size: " << *db_size << "\n\n";

    CRef<IQueryFactory> query_factory(new CObjMgr_QueryFactory(query_loc));

    CLocalBlast blaster(query_factory, opts, target_db);

    CSearchResultSet results = *blaster.Run();

    // print hits
    ifstream file;
    string filename = args["db"].AsString() + ".fasta";
    file.open(filename);
    int err = 0;
    if (file.fail()) {
        cerr << "Could not open file " << filename << endl;
        err = 1;
    }
    else {
        file.close();
    }

    if (!err) {
        int seq_id = 0;
        int curr_line = 0;
        string line;

        out << "Read name: <read_name>\nHit at pos <start-end>\nEvalue: <evalue> \n\n";
        
        for (unsigned int i = 0; i < results.GetNumResults(); i++) {
            CConstRef<CSeq_align_set> sas = results[i].GetSeqAlign();
            // cout << MSerial_AsnText << *sas;

            const list <CRef<CSeq_align>> &seqAlignList = sas->Get();
            for (list <CRef<CSeq_align>>::const_iterator seqAlign_it = ((ncbi::s_ITERATE_ConstRef(seqAlignList)).begin()); 
                seqAlign_it != ncbi::s_ITERATE_ConstRef(seqAlignList).end(); ++seqAlign_it ) {
                file.open(filename);
                seq_id = 0;
                curr_line = 0;
                sscanf((*seqAlign_it)->GetSeq_id(1).GetSeqIdString().c_str(), "%*[^:]:%d", &seq_id);
                while (!file.eof()) {
                    ++curr_line;
                    getline(file, line);
                    if (curr_line == (seq_id * 2) - 1) break;
                }
                file.close();
                out << "Read name: " << line << "\n"
                    << "Hit at pos " << (*seqAlign_it)->GetSeqStart(1) << "-" << (*seqAlign_it)->GetSeqStop(1) << "\n";
                double evalue;
                (*seqAlign_it)->GetNamedScore(CSeq_align::eScore_EValue, evalue);
                out << "E-value: " << evalue << "\n\n";
            }
        }
    }

    // Make a vector the size of the number of reads with each read name in ordered location
    // Apply p-value filter

    return results;
}

/////////////////////////////////////////////////////////////////////////////
//  Run test (printout arguments obtained from command-line)

int CBlastn::Run(void)
{
    // Get arguments
    const CArgs& args = GetArgs();

    EProgram program = ProgramNameToEnum("blastn");

    CRef<CBlastOptionsHandle> opts(CBlastOptionsFactory::Create(program));

    ProcessCommandLineArgs(opts);

    opts->Validate();  // Can throw CBlastException::eInvalidOptions for invalid option.

    // This will dump the options to stderr.
    // opts->GetOptions().DebugDumpText(cerr, "opts", 1);

    CRef<CObjectManager> objmgr = CObjectManager::GetInstance();
    if (!objmgr) {
         throw std::runtime_error("Could not initialize object manager");
    }

    SDataLoaderConfig dlconfig(false); // not working with proteins
    CBlastInputSourceConfig iconfig(dlconfig);
    CBlastFastaInputSource fasta_input(args["in"].AsInputFile(), iconfig);
    CScope scope(*objmgr);

    CBlastInput blast_input(&fasta_input);

    TSeqLocVector query_loc = blast_input.GetAllSeqLocs(scope);

    Uint8 db_size;
    CSearchResultSet results = full_blast(query_loc, opts, &db_size);

    // output number of BLAST hits
    int num_results = 0;
    for (unsigned int i = 0; i < results.GetNumResults(); i++) {
        const list <CRef<CSeq_align>> &seqAlignList = results[i].GetSeqAlign()->Get();
        num_results += seqAlignList.size();
    }
    CNcbiOstream& out = args["out"].AsOutputFile();
    out << "num_results: " << num_results << "\n";

    return 0;
}


/////////////////////////////////////////////////////////////////////////////
//  Cleanup

void CBlastn::Exit(void)
{
    SetDiagStream(0);
}


/////////////////////////////////////////////////////////////////////////////
//  MAIN

#ifndef SKIP_DOXYGEN_PROCESSING
int NcbiSys_main(int argc, char* argv[])
{
    // Execute main application function
    return CBlastn().AppMain(argc, argv, 0, eDS_Default, 0);
}
#endif /* SKIP_DOXYGEN_PROCESSING */
