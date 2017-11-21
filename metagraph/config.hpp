#ifndef __CONFIG_HPP__
#define __CONFIG_HPP__

#include <cstring>
#include <string>
#include <vector>

class Config {
    public:
        Config();
        Config(int argc, const char *argv[]);

        ~Config();

        void print_usage(std::string prog_name, int identity = -1);
        //void print_call(string prog_name);

        bool verbose;
        bool quiet;
        bool integrate;
        bool print_graph;
        bool print_graph_succ;
        bool reverse;
        bool fast;
        bool add_anno;

        unsigned int k;
        unsigned int distance;
        unsigned int parallel;
        unsigned int bins_per_thread;
        unsigned int parts_total;
        unsigned int part_idx;
        unsigned int collect;
        unsigned int frequency;
        unsigned int nsplits;

        std::vector<std::string> fname;
        std::string outfbase;
        std::string infbase;
        std::string sqlfbase;
        std::string dbpath;
        std::string refpath;
        std::string suffix;

        enum identities {noidentity = -1,
                         build = 1,
                         merge,
                         compare,
                         align,
                         stats,
                         annotate,
                         dump,
                         classify};
        int identity;

        enum state_type {stat = 1,
                         dyn,
                         cstr};
        state_type state;

    private:
        void init();
};

#endif
