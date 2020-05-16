#include <iostream>

//local libraries
#include "Contig.h"
#include "MPCollection.h"

//Graph clases
#include "GraphS.h"
#include "TReduction.h"
#include "GreedyS.h"
#include "MatchingS.h"
#include "G2Seq.h"
#include "ValidateBackbone.h"
#include "BFilling.h"
#include "GPolisher.h"
#include "SPolisher.h"

//options for liger
#include "lgopts.h"

int main (int argc, char* argv[]){

    gengetopt_args_info ai;
    if (cmdline_parser (argc, argv, &ai) != 0) {
        //cmdline_parser_print_help ();
        std::cout << "Run "<<argv[0]<<" -h to see the list of options. "<<std::endl;
        return EXIT_FAILURE;
    }


    string prefix=ai.prefix_arg;//prefix for output
    int cores=ai.cpu_arg;
    //load contigs from fai file
    auto contigs=new Contig(ai.contigs_arg);
    //contigs->print_contigs_file(prefix);
    //creates and reads the libraries
    auto libs=new MPCollection(ai.samlist_arg);
    //we read the libraries
    libs->read_libs(contigs,cores);
    //we check that the links were loaded correctly
    //libs->print_link_libs();
    //we use the contig coverage from short-reads and we use it to mark the repetitive contigs;
    if(ai.ccoverage_given){
            contigs->mark_repeats_astat(ai.ccoverage_arg,ai.rcn_arg);
    }else{
        //we compute the contig coverage and we mark the ones that are repeats from long reads
        contigs->mark_repeats();
    }

    //we are ready to create the scaffolding graph using the libraries and the contigs
    auto gori=new GraphS(contigs,libs,ai.mcs_arg);
    //This object perform the transitive reduction of the graph it compute the biconnected components and the search for paths
    TReduction tr;
    auto reduced=tr.transitive_reduction(ai.lme_arg,ai.mit_arg,gori);

    //Matching Cover of the graph
    MatchingS mc;
    G2Seq g2s;

    //we compute various parameters
    auto avg_ctg=contigs->getAvg_ctg_cov();
    //default run default
    auto mcg=mc.matching_cover(reduced,ai.rct_arg*avg_ctg,ai.mcr_arg);
    //we write an AGP file for each cover
    g2s.graph2agp(mcg,prefix+".MCover.raw.agp");
    //we obtained the chains
    //we create the Validation object
    ValidateBackbone eval;
    //matching cover evaluation
    auto vmcg=eval.Check_Backbone(ai.mlp_arg,ai.nlm_arg,mcg, tr.get_reduced_edges(), mc.get_circular_nodes());
    //we write an AGP file of the validated lines
    g2s.graph2agp(vmcg,prefix+".MCover.val.agp");
    //class that fill the gaps in the backbone it use SPOA and FASTER-SMITh-WATERMAN, coupled of with alignment-free and graph methods.
    //We fill the GAPs using the long-reads consensus and the reduced Graph to place repetitive contigs
    auto gapfiller=new BFilling(vmcg,ai.longreads_arg,prefix);
    //we create the edge consensus unsing the number of cores specified
    gapfiller->create_edge_consensus(cores);
    //then links are not available anymore
    libs->clear_links();
    //Graph polisher
    auto polishing = new GPolisher(reduced,vmcg,cores);
    polishing->polish();
     //Contigs using in the GPolisher step
    auto ctg_used_in_gpolished=polishing->get_ctg_used_in_polishing();
    //SPolisher
    auto sp=new SPolisher(ai.minimizer_size_arg,ai.minimizer_window_arg,ai.minimizer_freq_arg,vmcg);//k,w.mf
    //we build the kmer-index
    sp->build_edge_index();
    //we map the short-contig on the edge and select the best layout for each mate-edge
    sp->polish(ctg_used_in_gpolished);
    //we write the fasta file of the SPolished sequence and include the singleton contigs > 5kb
    g2s.graph2seq(vmcg,prefix+".SPolished",true);
    //hybrid assembly done
    return 0;
}
