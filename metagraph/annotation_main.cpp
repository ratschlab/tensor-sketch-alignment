#include <fstream>
#include <ctime>
#include <zlib.h>

#include "dbg_hash.hpp"
#include "config.hpp"
#include "helpers.hpp"
#include "utils.hpp"
#include "vcf_parser.hpp"
#include "annotate.hpp"
#include "unix_tools.hpp"


KSEQ_INIT(gzFile, gzread);


const std::vector<std::string> annots = {
  "AC_AFR", "AC_EAS", "AC_AMR", "AC_ASJ",
  "AC_FIN", "AC_NFE", "AC_SAS", "AC_OTH"
};


/*
 * Distribute the annotation of all k-mers in a sequence over
 * a number of parallel bins.
 */
// void* parallel_annotate_wrapper(void *) {

//     uint64_t curr_idx, start, end;

//     while (true) {
//         pthread_mutex_lock (&mutex_bin_idx);
//         if (anno_data->idx == anno_data->total_bins) {
//             pthread_mutex_unlock (&mutex_bin_idx);
//             break;
//         } else {
//             curr_idx = anno_data->idx;
//             anno_data->idx++;
//             pthread_mutex_unlock (&mutex_bin_idx);

//             start = curr_idx * anno_data->binsize;
//             //end = std::min(((curr_idx + 1) * anno_data->binsize) + anno_data->graph->k - 1, anno_data->seq->l);
//             end = std::min((curr_idx + 1) * anno_data->binsize,
//                            static_cast<uint64_t>(anno_data->seq->l));
//             //std::cerr << "start " << start << " end " << end << std::endl;
//             annotate::annotate_seq(anno_data->graph,
//                                    anno_data->config,
//                                    *(anno_data->seq),
//                                    *(anno_data->label),
//                                    start, end, &mutex_annotate);
//         }
//     }
//     pthread_exit((void*) 0);
// }


DBG_succ* load_critical_graph_from_file(const std::string &filename) {
    std::string filetype = ".dbg";
    std::string filebase =
        std::equal(filetype.rbegin(), filetype.rend(), filename.rbegin())
        ? filename.substr(0, filename.size() - filetype.size())
        : filename;

    DBG_succ *graph = new DBG_succ();
    if (!graph->load(filebase)) {
        std::cerr << "ERROR: input file "
                  << filename << " corrupted" << std::endl;
        delete graph;
        exit(1);
    }
    return graph;
}


int main(int argc, const char *argv[]) {

    Timer timer;

    // parse command line arguments and options
    std::unique_ptr<Config> config(new Config(argc, argv));

    if (config->verbose) {
        std::cout << "#############################\n"
                  << "### Welcome to AnnoGraph! ###\n"
                  << "#############################\n" << std::endl;
    }

    const auto &files = config->fname;

    if (config->identity != Config::BUILD) {
        std::cerr << "Error: Only BUILD mode is supported" << std::endl;
        exit(1);
    }

    annotate::BloomAnnotator *annotator = NULL;
    annotate::PreciseAnnotator *precise_annotator = NULL;

    DBGHash hashing_graph(config->k);

    if (config->bloom_fpp > -0.5
            || config->bloom_bits_per_edge > -0.5
            || config->bloom_num_hash_functions > 0) {
        if (config->bloom_fpp > -0.5) {
            // Expected FPP is set, optimize other parameters automatically
            annotator = new annotate::BloomAnnotator(hashing_graph,
                                                     config->bloom_fpp,
                                                     config->verbose);
        } else {
            assert(config->bloom_bits_per_edge >= 0);
            // Experiment mode, estimate FPP given other parameters,
            // optimize the number of hash functions if it's set to zero
            annotator = new annotate::BloomAnnotator(
                hashing_graph,
                config->bloom_bits_per_edge,
                config->bloom_num_hash_functions,
                config->verbose
            );
        }
        if (config->bloom_test_num_kmers > 0) {
            precise_annotator = new annotate::PreciseAnnotator(
                hashing_graph
            );
        }
    }
    if (annotator) {
        std::cout << "Bloom filter settings" << std::endl;
        std::cout << "\tBits per edge:\t" << annotator->size_factor() << std::endl;
        std::cout << "\tNum hash functions:\t" << annotator->num_hash_functions() << std::endl;
        std::cout << "\tApprox false pos prob:\t" << annotator->approx_false_positive_rate() << std::endl;
    }

    if (config->verbose)
        std::cerr << "k is " << config->k << std::endl;

    if (config->fast) {
        //enumerate all suffices
        assert(DBG_succ::alph_size > 1);
        size_t suffix_len = 0;
        auto suffices = utils::generate_strings(
            DBG_succ::alphabet.substr(0, DBG_succ::alph_size),
            suffix_len
        );

        //one pass per suffix
        double file_read_time = 0;
        double graph_const_time = 0;
        double bloom_const_time = 0;
        double precise_const_time = 0;
        std::map<std::string, size_t> annot_map;
        for (const std::string &suffix : suffices) {
            if (suffix.size())
                std::cout << "Suffix: " << suffix << std::endl;

            std::cout << "Start reading data and extracting k-mers..." << std::endl;

            //add sink nodes
            // graph->add_sink(config->parallel, suffix);

            // iterate over input files
            for (unsigned int f = 0; f < files.size(); ++f) {
                if (config->verbose) {
                    std::cout << std::endl << "Parsing " << files[f] << std::endl;
                }

                if (utils::get_filetype(files[f]) == "VCF") {
                    if (suffix.find('$') != std::string::npos)
                        continue;

                    //READ FROM VCF
                    Timer data_reading_timer;

                    vcf_parser vcf;
                    if (!vcf.init(config->refpath, files[f], config->k)) {
                        std::cerr << "ERROR reading VCF " << files[f] << std::endl;
                        exit(1);
                    }
                    std::cerr << "Loading VCF with " << config->parallel
                                                     << " threads per line\n";
                    std::string sequence;
                    std::vector<std::string> annotation;
                    data_reading_timer.reset();
                    for (size_t i = 1; vcf.get_seq(annots, &sequence, annotation); ++i) {
                        file_read_time += data_reading_timer.elapsed();
                        if (config->verbose && i % 10'000 == 0) {
                            std::cout << "." << std::flush;
                            if (i % 1'000'000 == 0) {
                                std::cout << std::endl;
                            }
                        }
                        for (size_t j = 0; j < 2; ++j) {
                            data_reading_timer.reset();
                            hashing_graph.add_sequence(sequence, true);
                            graph_const_time += data_reading_timer.elapsed();

                            for (auto it = annotation.begin(); it != annotation.end(); ++it) {
                                auto map_ins = annot_map.insert(
                                        std::pair<std::string,size_t>(
                                            *it,
                                            annot_map.size()));
                                data_reading_timer.reset();
                                if (map_ins.second) {
                                    annotator->add_column(
                                            sequence,
                                            (sequence.length() - hashing_graph.get_k())
                                                * (static_cast<size_t>(config->reverse) + 1)
                                    );
                                } else {
                                    annotator->add_sequence(
                                            sequence,
                                            map_ins.first->second,
                                            (sequence.length() - hashing_graph.get_k())
                                                * (static_cast<size_t>(config->reverse) + 1)
                                    );
                                }
                                bloom_const_time += data_reading_timer.elapsed();
                                if (precise_annotator) {
                                    data_reading_timer.reset();
                                    if (map_ins.second) {
                                        precise_annotator->add_column(sequence);
                                    } else {
                                        precise_annotator->add_sequence(sequence, map_ins.first->second);
                                    }
                                    precise_const_time += data_reading_timer.elapsed();
                                }
                            }
                            annotation.clear();
                            if (config->reverse) {
                                reverse_complement(sequence.begin(), sequence.end());
                            }
                            else
                                break;
                        }
                    }
                } else if (utils::get_filetype(files[f]) == "FASTA"
                            || utils::get_filetype(files[f]) == "FASTQ") {
                    // open stream
                    gzFile input_p = gzopen(files[f].c_str(), "r");
                    if (input_p == Z_NULL) {
                        std::cerr << "ERROR no such file " << files[f] << std::endl;
                        exit(1);
                    }

                    //TODO: handle read_stream->qual
                    kseq_t *read_stream = kseq_init(input_p);
                    if (read_stream == NULL) {
                        std::cerr << "ERROR while opening input file "
                                  << files[f] << std::endl;
                        exit(1);
                    }
                    Timer result_timer;
                    bool fastq = utils::get_filetype(files[f]) == "FASTQ";
                    std::pair<std::map<std::string, size_t>::iterator, bool> map_ins;
                    if (fastq) {
                        map_ins = annot_map.insert(
                                std::pair<std::string,size_t>(
                                    files[f],
                                    annot_map.size()));
                    }

                    while (kseq_read(read_stream) >= 0) {
                        {
                            file_read_time += result_timer.elapsed();

                            for (size_t j = 0; j < 2; ++j) {

                                result_timer.reset();
                                hashing_graph.add_sequence(read_stream->seq.s);
                                graph_const_time += result_timer.elapsed();

                                //constructor->add_read(read_stream->seq.s);
                                if (fastq) {
                                    //assume each FASTQ file is one column

                                    result_timer.reset();
                                    //TODO: for reads this leads to a much smaller filter than needed
                                    //TODO: tell the annotator to allocate a much larger space for reads
                                    annotator->add_sequence(
                                            read_stream->seq.s,
                                            map_ins.first->second);
                                    bloom_const_time += result_timer.elapsed();

                                    if (precise_annotator) {
                                        result_timer.reset();
                                        precise_annotator->add_sequence(read_stream->seq.s, map_ins.first->second);
                                        precise_const_time += result_timer.elapsed();
                                    }
                                } else {
                                    //assume each sequence in each FASTA file is a column
                                    result_timer.reset();
                                    map_ins = annot_map.insert(
                                        std::pair<std::string,size_t>(
                                            std::string(read_stream->name.s),
                                            annot_map.size()));
                                    if (map_ins.second)
                                        annotator->add_column(
                                                read_stream->seq.s,
                                                (read_stream->seq.l - hashing_graph.get_k())
                                                    * (static_cast<size_t>(config->reverse) + 1)
                                        );
                                    else
                                        annotator->add_sequence(
                                                read_stream->seq.s,
                                                map_ins.first->second,
                                                (read_stream->seq.l - hashing_graph.get_k())
                                                    * (static_cast<size_t>(config->reverse) + 1)
                                        );
                                    bloom_const_time += result_timer.elapsed();

                                    if (precise_annotator) {
                                        result_timer.reset();
                                        if (map_ins.second)
                                            precise_annotator->add_column(read_stream->seq.s);
                                        else
                                            precise_annotator->add_sequence(read_stream->seq.s, map_ins.first->second);
                                        precise_const_time += result_timer.elapsed();
                                    }
                                }
                                if (config->reverse)
                                    reverse_complement(read_stream->seq);
                                else
                                    break;
                            }
                        }
                        result_timer.reset();
                    }
                    kseq_destroy(read_stream);
                    gzclose(input_p);
                } else {
                    std::cerr << "ERROR: Filetype unknown for file "
                              << files[f] << std::endl;
                    exit(1);
                }
                //fprintf(stdout, "current mem usage: %lu MB\n", get_curr_mem() / (1<<20));
            }
            get_RAM();
            std::cout << "Reading data finished\t" << timer.elapsed() << "sec" << std::endl;
        }

        //Runtime stats
        std::cout << "Runtime statistics" << std::endl;
        std::cout << "File reading\t" << file_read_time << std::endl;
        std::cout << "Graph construction\t" << graph_const_time << std::endl;
        std::cout << "Bloom filter\t" << bloom_const_time << std::endl;
        std::cout << "Precise filter\t" << precise_const_time << std::endl;

        if (precise_annotator) {
            //Check FPP
            std::cout << "Approximating FPP...\t" << std::flush;
            timer.reset();
            //TODO: set this value
            annotator->test_fp_all(*precise_annotator, config->bloom_test_num_kmers);
            std::cout << timer.elapsed() << "sec" << std::endl;
        }

        if (config->state == Config::DYN) {
            std::cerr << "Converting static graph to dynamic...\t" << std::flush;
            timer.reset();
            std::cout << timer.elapsed() << "sec" << std::endl;
        }
    } else {
        //slower method
        //TODO: merge in optimizations from seqmerge branch
        //TODO: add VCF parsing (needs previous merge)
        for (unsigned int f = 0; f < files.size(); ++f) {
            if (config->verbose) {
                std::cout << std::endl << "Parsing " << files[f] << std::endl;
            }
            // open stream
            if (utils::get_filetype(files[f]) == "VCF") {
                std::cerr << "ERROR: this method of reading VCFs not yet implemented" << std::endl;
                exit(1);
            } else if (utils::get_filetype(files[f]) == "FASTA"
                        || utils::get_filetype(files[f]) == "FASTQ") {
                gzFile input_p = gzopen(files[f].c_str(), "r");
                if (input_p == Z_NULL) {
                    std::cerr << "ERROR no such file " << files[f] << std::endl;
                    exit(1);
                }
                kseq_t *read_stream = kseq_init(input_p);
                if (read_stream == NULL) {
                    std::cerr << "ERROR while opening input file "
                              << files[f] << std::endl;
                    exit(1);
                }
                for (size_t i = 1; kseq_read(read_stream) >= 0; ++i) {

                    if (config->reverse) {
                        reverse_complement(read_stream->seq);
                    }
                }
                kseq_destroy(read_stream);
                gzclose(input_p);
            } else {
                std::cerr << "ERROR: Filetype unknown for file "
                          << files[f] << std::endl;
                exit(1);
            }
        }
    }

    std::cout << "Graph size(edges):\t" << hashing_graph.get_num_edges() << std::endl;

    config->infbase = config->outfbase;
    // graph->annotationToFile(config->infbase + ".anno.dbg");
    //graph->print_state();

    // output and cleanup

    // graph output
    if (!config->outfbase.empty() && annotator)
        annotator->serialize(config->outfbase);

    if (annotator)
        delete annotator;
    if (precise_annotator)
        delete precise_annotator;

    return 0;
}