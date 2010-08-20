/*
  This getopt_long() is compatible with GNU's, however, added original
  extention (short 1 byte option).


  Copyright (c) 2004 Koji Arai

  Permission is hereby granted, free of charge, to any person
  obtaining a copy of this software and associated documentation files
  (the "Software"), to deal in the Software without restriction,
  including without limitation the rights to use, copy, modify, merge,
  publish, distribute, sublicense, and/or sell copies of the Software,
  and to permit persons to whom the Software is furnished to do so,
  subject to the following conditions:

  The above copyright notice and this permission notice shall be
  included in all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
  BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
  ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
  CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.


  Compilation for Test:

      GNU:
      cc -DUSE_GNU -DDEBUG getopt_long.c -o test_getopt_long_gnu

      not GNU:
      cc -I. -DDEBUG getopt_long.c -o test_getopt_long

      ./test_getopt_long
      ./test_getopt_long_gnu

  BUGS:
    * not implemented any features for getopt() and getopt_long().
*/

#include <assert.h>
#include <stdio.h>
#include <string.h>

#if DEBUG
static int
puts_argv(char **argv)
{
    int i;

    for (i = 0; argv[i]; i++) {
        if (i) printf(" ");

        printf("%s", argv[i]);
    }
    printf("\n");

    return 0;
}
#endif

#ifndef USE_GNU
#include <stdio.h>
#include "getopt_long.h"

char *optarg;
int optind;

int opterr;
int optopt;

/*
  return value 0: no option (include '-')
               1: short option like '-x'
               2: long option like '--xxx' and just '--'
*/
static int
is_option(char *arg)
{
    if (arg[0] == '-') {
        switch (arg[1]) {
        case 0:                 /* just "-" */
            return 0;
        case '-':               /* long option (include just "--")*/
            return 2;
        default:                /* short option */
            return 1;
        }
    }
    return 0;
}

static int
insert_argv(char **argv, int src, int dest)
{
    int i;
    char *tmp = argv[src];

    if (src > dest) {
        for (i = src; i > dest; i--)
            argv[i] = argv[i-1];
    }
    if (src < dest) {
        for (i = src; i < dest; i++)
            argv[i] = argv[i+1];
    }

    argv[dest] = tmp;

    return 0;
}

static int
search_longopt(char *arg, struct option *longopts)
{
    int i, found = -1;
    int len;
    for (len = 0; arg[len] && arg[len] != '='; len++)
        ;

    for (i = 0; longopts[i].name; i++) {
        if (strncmp(arg, longopts[i].name, len) == 0) {
            if (found != -1)
                return -1;      /* found some candidate */
            found = i;
        }
    }
    return found;
}

/*
 * implemented my extention feature.
 * optional 1 byte argument with [...]
 *   e.g.) shortopts = "a[0123]b"
 *          accepts "-a0 -a1b" (same as "-a0 -a1 -b")
 */
static int
has_argument_short(char *arg, const char *shortopts)
{
    int i;
    int open_bracket = 0;
    for (i = 0; shortopts[i]; i++) {
        switch (shortopts[i]) {
        case '[':
            open_bracket++;
            continue;
        case ']':
            if (open_bracket <= 0) {
                fprintf(stderr, "getopt_long() -- unbalanced bracket in short options");
                return -1;
            }
            open_bracket--;
            continue;
        }
        if (open_bracket) continue;
        if (*arg != shortopts[i]) continue;

        switch (shortopts[i+1]) {
        case ':':
            if (shortopts[i+2] != ':') {
                if (arg[1])
                    return 1; /* following string is argument */
                else
                    return 2; /* next argv is argument */
            }
            else {
                /* '::' means optional argument (GNU extention) */
                if (arg[1])
                    return 1;
                else
                    return 0; /* no argument */
            }
        case '[':
            if (arg[1] == '\0')
                return 0;   /* no argument */
            /* my extention */
            for (i++; shortopts[i] && shortopts[i] != ']'; i++) {
                if (arg[1] == shortopts[i])
                    return 3; /* has 1 byte argument */
            }
            if (!shortopts[i]) {
                fprintf(stderr, "getopt_long() -- unbalanced bracket in short options");
                return -1;
            }
            break;
        default:
            return 0;   /* no argument */
        }
    }
    /* Invalid option */
    return -1;
}

static int
has_argument_long(char *arg, struct option *longopts)
{
    int i;

    i = search_longopt(arg, longopts);
    if (i == -1) {
        /* Invalid option */
        return -1;
    }
    else {
        int len = strlen(arg);
        char *p = strchr(arg, '=');
        if (p) {
            len = p - arg;
        }

        switch (longopts[i].has_arg) {
        case no_argument:
            return 0;
        case required_argument:
            if (arg[len] == '=')
                return 1;
            else
                return 2;
        case optional_argument:
            if (arg[len] == '=')
                return 1;
            else
                return 0;
        default:
            assert(0);
        }
    }
}

/*
  -1: no option
   0: no argument
   1: has argument in this argv
   2: has argument in next argv
   3: has 1 byte argument in this argv
*/
static int
has_argument(char *arg,
             const char *shortopts,
             struct option *longopts)
{
    int i, n;

    switch (is_option(arg)) {
    case 0:                     /* no option */
        return -1;
    case 1:
        /* short option */
        n = -1;
        for (i = 1; arg[i]; i++) {
            n = has_argument_short(arg+i, shortopts);
            if (n == 0 && arg[i+1]) continue;
            if (n == 3 && arg[i+2]) { i++; continue; }
            break;
        }
        return n;
    case 2:
        /* long option */
        return has_argument_long(arg+2, longopts);
        break;
    default:
        assert(0);
    }
}

int
getopt_long(int argc, char **argv,
            const char *shortopts,
            struct option *longopts,
            int *indexptr)
{
    char *opt;
    int i;
    static int shortoptind;
    static int no_optind = 0;
    
    if (optind == 0) {            /* skip first argument (command name) */
        optind++;
        no_optind = 0;
        shortoptind = 0;
    }

    optarg = 0;

    if (no_optind && !shortoptind) {
        while (!is_option(argv[no_optind]))
            insert_argv(argv, no_optind, optind-1);

        if (has_argument(argv[no_optind], shortopts, longopts) == 2)
            no_optind += 2;
        else
            no_optind++;

        if (argv[optind] && strcmp(argv[optind], "--") == 0) {
            while (!is_option(argv[no_optind]))
                insert_argv(argv, no_optind, optind);
            optind = no_optind;
            no_optind = 0;
        }
    }

    if (optind >= argc)
        goto end_of_option;

 retry:
    /*
    puts_argv(&argv[optind]);
    */
    opt = argv[optind];

    if (shortoptind == 0 && is_option(opt) == 1) {
        shortoptind++;
    }

    if (shortoptind) {
        /* short option */
        char *p = &opt[shortoptind];

        if (*p == '\0')
            assert(0);

        switch (has_argument_short(p, shortopts)) {
        case 0:
            /* no argument */
            optarg = 0;

            shortoptind++;
            if (opt[shortoptind] == '\0')
                optind++, shortoptind = 0;
            return *p;
        case 1:
            /* following character is argument */
            optind++, shortoptind = 0;
            optarg = &p[1];
            return *p;
        case 2:
            /* next argv is argument */
            optind++, shortoptind = 0;
            optarg = argv[optind++];
            return *p;
        case 3:
            /* has 1 byte argument */
            optarg = &p[1];
            if (p[2] == 0)
                optind++, shortoptind = 0;
            else
                shortoptind += 2;
            return *p;
        default:
            /* Invalid option */
            if (opterr)
                fprintf(stderr,
                        "%s: invalid option -- %c\n",
                        argv[0],
                        *p);

            optind++, shortoptind = 0;
            optopt = *p;
            return '?';
        }
    }
    else if (opt[0] == '-' && opt[1] == '-') {
        /* long option */

        if (opt[2] == '\0') {
            /* end of command line switch */
            optind++;
            return -1;
        }

        opt += 2;

        i = search_longopt(opt, longopts);
        if (i == -1) {
            optind++;
            optopt = 0;
            return '?';
        }
        else {
            int len = strlen(opt);
            char *p = strchr(opt, '=');
            if (p) {
                len = p - opt;
            }

            switch (longopts[i].has_arg) {
            case no_argument:
                break;
            case required_argument:
                if (opt[len] == '=')
                    optarg = opt + len + 1;
                else {
                    optind++;
                    optarg = argv[optind];
                    if (optarg == 0) {
                        if (opterr)
                            fprintf(stderr,
                                    "%s: option `--%s' requires an argument\n",
                                    argv[0],
                                    opt);

                        optopt = 0;
                        return '?'; /* no argument */
                    }
                }
                break;
            case optional_argument:
                if (opt[len] == '=')
                    optarg = opt + len + 1;
                else {
                    optarg = 0;
                }
                break;
            default:
                break;
            }

            *indexptr = i;
            optind++;
            if (longopts[i].flag) {
                *longopts[i].flag = longopts[i].val;
                return 0;
            }
            else {
                return longopts[i].val;
            }
        }

        optind++;
        optopt = 0;
        return '?';
    }

    /* not option */
    if (no_optind == 0)
        no_optind = optind;

    for (i = optind; argv[i]; i++) {
        if (is_option(argv[i])) {
            optind = i;
            goto retry;
        }
    }

 end_of_option:
    if (no_optind) {
        optind = no_optind;
        no_optind = 0;
    }

    return -1;
}
#endif /* USE_GNU */

#if DEBUG

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#if USE_GNU
#include <getopt.h>  /* use GNU getopt_long() */
#endif

static int verbose_flag;
static int option_index;
int argc;
char *argv[50];
char **p;
int c;
static struct option long_options[] = {
    {"verbose", no_argument, &verbose_flag, 1},
    {"brief", no_argument, &verbose_flag, 0},
    {"add", required_argument, 0, 'a'},
    {"append", no_argument, 0, 0},
    {"delete", required_argument, 0, 0},
    {"create", optional_argument, 0, 0},
    {"change", optional_argument, 0, 0},
    {0, 0, 0, 0}
};

int
call_getopt_long(int argc, char **argv,
                 const char *shortopts,
                 struct option *longopts,
                 int *indexptr)
{
    int c;
    c = getopt_long(argc, argv, shortopts, longopts, indexptr);
    puts_argv(argv);
    printf("ret=%d(%c) option_index=%d ", c, c, option_index);
    printf("optind=%d optarg=[%s] opterr=%d optopt=%d(%c)\n",
           optind, optarg, opterr, optopt, optopt);
    if (c == 0) {
        struct option *opt;
        opt = &longopts[*indexptr];
        printf("long option: --%s has_arg=%d\n", opt->name, opt->has_arg);
        if (opt->flag)
            printf("           flag=[%8p] val=%d\n", opt->flag, *opt->flag);
    }

    return c;
}

void
test1()
{
    optind = 0;
    argc = 0;
    p = argv;

    argc++; *p++ = "command_name";
    argc++; *p++ = "-a";
    argc++; *p++ = "-bcd";
    argc++; *p++ = "-d";
    argc++; *p++ = "-e";
    argc++; *p++ = "-f";
    argc++; *p++ = "-g";
    *p = 0;

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    assert(c == 'a');
    assert(option_index == 0);
    assert(optind == 2);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    assert(c == 'b');
    assert(option_index == 0);
    assert(optind == 2);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    assert(c == 'c');
    assert(option_index == 0);
    assert(optind == 3);
    assert(optarg == &argv[2][3]);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    assert(c == 'd');
    assert(option_index == 0);
    assert(optind == 5);
    assert(optarg == argv[4]);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    assert(c == '?');
    assert(option_index == 0);
    assert(optind == 6);
    assert(optarg == 0);
    assert(optopt == 'f');

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    assert(c == '?');
    assert(option_index == 0);
    assert(optind == 7);
    assert(optarg == 0);
    assert(optopt == 'g');

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    assert(c == -1);
    assert(option_index == 0);
    assert(optind == 7);
    assert(optarg == 0);
    assert(optopt == 'g');      /* no changed */
}

void
test2()
{
    optind = 0;
    argc = 0;
    p = argv;

    argc++; *p++ = "command_name";
    argc++; *p++ = "--verbose";
    argc++; *p++ = "--brief";
    argc++; *p++ = "--add";
    argc++; *p++ = "add_argument";
    argc++; *p++ = "--add=add_argument";
    argc++; *p++ = "--append";
    argc++; *p++ = "--delete=del_argument";
    argc++; *p++ = "--create=cre_argument";
    argc++; *p++ = "--create";
    argc++; *p++ = "files...";
    *p = 0;

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    assert(c == 0);
    assert(option_index == 0);
    assert(optind == 2);
    assert(optarg == 0);
    assert(optopt == 'g');      /* no changed */
    assert(strcmp(long_options[option_index].name, "verbose") == 0);
    assert(*long_options[option_index].flag == 1);

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    assert(c == 0);
    assert(option_index == 1);
    assert(optind == 3);
    assert(optarg == 0);
    assert(optopt == 'g');      /* no changed */
    assert(strcmp(long_options[option_index].name, "brief") == 0);
    assert(*long_options[option_index].flag == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    assert(c == 'a');
    assert(option_index == 2);
    assert(optind == 5);
    assert(optarg == argv[4]);
    assert(optopt == 'g');      /* no changed */
    assert(strcmp(long_options[option_index].name, "add") == 0);
    assert(long_options[option_index].flag == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    assert(c == 'a');
    assert(option_index == 2);
    assert(optind == 6);
    assert(optarg == argv[5]+6);
    assert(optopt == 'g');      /* no changed */
    assert(strcmp(long_options[option_index].name, "add") == 0);
    assert(long_options[option_index].flag == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    assert(c == 0);
    assert(option_index == 3);
    assert(optind == 7);
    assert(optarg == 0);
    assert(optopt == 'g');      /* no changed */
    assert(strcmp(long_options[option_index].name, "append") == 0);
    assert(long_options[option_index].flag == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    assert(c == 0);
    assert(option_index == 4);
    assert(optind == 8);
    assert(optarg == argv[7]+9);
    assert(optopt == 'g');      /* no changed */
    assert(strcmp(long_options[option_index].name, "delete") == 0);
    assert(long_options[option_index].flag == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    assert(c == 0);
    assert(option_index == 5);
    assert(optind == 9);
    assert(optarg == argv[8]+9);
    assert(optopt == 'g');      /* no changed */
    assert(strcmp(long_options[option_index].name, "create") == 0);
    assert(long_options[option_index].flag == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    assert(c == 0);
    assert(option_index == 5);
    assert(optind == 10);
    assert(optarg == 0);
    assert(optopt == 'g');      /* no changed */
    assert(strcmp(long_options[option_index].name, "create") == 0);
    assert(long_options[option_index].flag == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    assert(c == -1);
    assert(option_index == 5);  /* no changed */
    assert(optind == 10);
    assert(optarg == 0);
    assert(optopt == 'g');      /* no changed */
    assert(strcmp(argv[optind], "files...") == 0);

}

void
test3()
{
    optind = 0;
    argc = 0;
    p = argv;

    argc++; *p++ = "command_name";
    argc++; *p++ = "--delete";  /* required argument has no argument */
    *p = 0;

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    assert(c == '?');
    assert(option_index == 5);  /* no changed */
    assert(optind == 2);        /* changed */
    assert(optarg == 0);
    assert(optopt == 0);       /* changed */
    assert(argv[optind] == 0);

    /* */
    optind = 0;
    argc = 0;
    p = argv;

    argc++; *p++ = "command_name";
    argc++; *p++ = "--file";  /* not option */
    *p = 0;

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    assert(c == '?');
    assert(option_index == 5);  /* no changed */
    assert(optind == 2);
    assert(optarg == 0);
    assert(optopt == 0);
    assert(argv[optind] == 0);
}

void
test4()
{
    optind = 0;
    argc = 0;
    p = argv;

    argc++; *p++ = "command_name";
    argc++; *p++ = "-a";
    argc++; *p++ = "a1";
    argc++; *p++ = "a2";
    argc++; *p++ = "-b";
    argc++; *p++ = "b";
    argc++; *p++ = "-efg";      /* some options in a argument */
    argc++; *p++ = "g";
    argc++; *p++ = "-c";
    argc++; *p++ = "c";
    argc++; *p++ = "d";
    *p = 0;

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:efg:", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "a1") == 0);
    assert(strcmp(*p++, "a2") == 0);
    assert(strcmp(*p++, "-b") == 0);
    assert(strcmp(*p++, "b") == 0);
    assert(strcmp(*p++, "-efg") == 0);
    assert(strcmp(*p++, "g") == 0);
    assert(strcmp(*p++, "-c") == 0);
    assert(strcmp(*p++, "c") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(*p == 0);

    assert(c == 'a');
    assert(option_index == 5);  /* no changed */
    assert(optind == 2);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:efg:", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "a1") == 0);
    assert(strcmp(*p++, "a2") == 0);
    assert(strcmp(*p++, "-b") == 0);
    assert(strcmp(*p++, "b") == 0);
    assert(strcmp(*p++, "-efg") == 0);
    assert(strcmp(*p++, "g") == 0);
    assert(strcmp(*p++, "-c") == 0);
    assert(strcmp(*p++, "c") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(*p == 0);

    assert(c == 'b');
    assert(option_index == 5);  /* no changed */
    assert(optind == 5);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:efg:", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-b") == 0);
    assert(strcmp(*p++, "a1") == 0);
    assert(strcmp(*p++, "a2") == 0);
    assert(strcmp(*p++, "b") == 0);
    assert(strcmp(*p++, "-efg") == 0);
    assert(strcmp(*p++, "g") == 0);
    assert(strcmp(*p++, "-c") == 0);
    assert(strcmp(*p++, "c") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(*p == 0);

    assert(c == 'e');
    assert(option_index == 5);  /* no changed */
    assert(optind == 6);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:efg:", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-b") == 0);
    assert(strcmp(*p++, "a1") == 0);
    assert(strcmp(*p++, "a2") == 0);
    assert(strcmp(*p++, "b") == 0);
    assert(strcmp(*p++, "-efg") == 0);
    assert(strcmp(*p++, "g") == 0);
    assert(strcmp(*p++, "-c") == 0);
    assert(strcmp(*p++, "c") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(*p == 0);

    assert(c == 'f');
    assert(option_index == 5);  /* no changed */
    assert(optind == 6);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:efg:", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-b") == 0);
    assert(strcmp(*p++, "a1") == 0);
    assert(strcmp(*p++, "a2") == 0);
    assert(strcmp(*p++, "b") == 0);
    assert(strcmp(*p++, "-efg") == 0);
    assert(strcmp(*p++, "g") == 0);
    assert(strcmp(*p++, "-c") == 0);
    assert(strcmp(*p++, "c") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(*p == 0);

    assert(c == 'g');
    assert(option_index == 5);  /* no changed */
    assert(optind == 8);
    assert(optarg == argv[7]);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:efg:", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-b") == 0);
    assert(strcmp(*p++, "-efg") == 0);
    assert(strcmp(*p++, "g") == 0);
    assert(strcmp(*p++, "a1") == 0);
    assert(strcmp(*p++, "a2") == 0);
    assert(strcmp(*p++, "b") == 0);
    assert(strcmp(*p++, "-c") == 0);
    assert(strcmp(*p++, "c") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(*p == 0);

    assert(c == 'c');
    assert(option_index == 5);  /* no changed */
    assert(optind == 10);
    assert(optarg == argv[9]);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:efg:", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-b") == 0);
    assert(strcmp(*p++, "-efg") == 0);
    assert(strcmp(*p++, "g") == 0);
    assert(strcmp(*p++, "-c") == 0);
    assert(strcmp(*p++, "c") == 0);
    assert(strcmp(*p++, "a1") == 0);
    assert(strcmp(*p++, "a2") == 0);
    assert(strcmp(*p++, "b") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(*p == 0);

    assert(c == -1);
    assert(option_index == 5);  /* no changed */
    assert(optind == 7);
    assert(optarg == 0);
    assert(optopt == 0);

}

void
test5()
{
    optind = 0;
    argc = 0;
    p = argv;

    argc++; *p++ = "command_name";
    argc++; *p++ = "-a";
    argc++; *p++ = "-";
    argc++; *p++ = "-b";
    *p = 0;

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-") == 0);
    assert(strcmp(*p++, "-b") == 0);
    assert(*p == 0);

    assert(c == 'a');
    assert(option_index == 5);  /* no changed */
    assert(optind == 2);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-") == 0);
    assert(strcmp(*p++, "-b") == 0);
    assert(*p == 0);

    assert(c == 'b');
    assert(option_index == 5);  /* no changed */
    assert(optind == 4);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-b") == 0);
    assert(strcmp(*p++, "-") == 0);
    assert(*p == 0);

    assert(c == -1);
    assert(option_index == 5);  /* no changed */
    assert(optind == 3);
    assert(optarg == 0);
    assert(optopt == 0);
}

void
test6()
{
    optind = 0;
    argc = 0;
    p = argv;

    argc++; *p++ = "command_name";
    argc++; *p++ = "-a";
    argc++; *p++ = "-";
    argc++; *p++ = "-";
    argc++; *p++ = "-b";
    *p = 0;

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-") == 0);
    assert(strcmp(*p++, "-") == 0);
    assert(strcmp(*p++, "-b") == 0);
    assert(*p == 0);

    assert(c == 'a');
    assert(option_index == 5);  /* no changed */
    assert(optind == 2);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-") == 0);
    assert(strcmp(*p++, "-") == 0);
    assert(strcmp(*p++, "-b") == 0);
    assert(*p == 0);

    assert(c == 'b');
    assert(option_index == 5);  /* no changed */
    assert(optind == 5);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-b") == 0);
    assert(strcmp(*p++, "-") == 0);
    assert(strcmp(*p++, "-") == 0);
    assert(*p == 0);

    assert(c == -1);
    assert(option_index == 5);  /* no changed */
    assert(optind == 3);
    assert(optarg == 0);
    assert(optopt == 0);
}

void
test7()
{
    optind = 0;
    argc = 0;
    p = argv;

    argc++; *p++ = "command_name";
    argc++; *p++ = "-a";
    argc++; *p++ = "-";
    argc++; *p++ = "-";
    argc++; *p++ = "-c";
    argc++; *p++ = "c";
    *p = 0;

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-") == 0);
    assert(strcmp(*p++, "-") == 0);
    assert(strcmp(*p++, "-c") == 0);
    assert(strcmp(*p++, "c") == 0);
    assert(*p == 0);

    assert(c == 'a');
    assert(option_index == 5);  /* no changed */
    assert(optind == 2);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-") == 0);
    assert(strcmp(*p++, "-") == 0);
    assert(strcmp(*p++, "-c") == 0);
    assert(strcmp(*p++, "c") == 0);
    assert(*p == 0);

    assert(c == 'c');
    assert(option_index == 5);  /* no changed */
    assert(optind == 6);
    assert(optarg == argv[5]);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-c") == 0);
    assert(strcmp(*p++, "c") == 0);
    assert(strcmp(*p++, "-") == 0);
    assert(strcmp(*p++, "-") == 0);
    assert(*p == 0);

    assert(c == -1);
    assert(option_index == 5);  /* no changed */
    assert(optind == 4);
    assert(optarg == 0);
    assert(optopt == 0);
}

void
test8()
{
    optind = 0;
    argc = 0;
    p = argv;

    argc++; *p++ = "command_name";
    argc++; *p++ = "-a";
    argc++; *p++ = "-c";
    argc++; *p++ = "c";
    argc++; *p++ = "--";
    argc++; *p++ = "-d";
    argc++; *p++ = "d";
    *p = 0;

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-c") == 0);
    assert(strcmp(*p++, "c") == 0);
    assert(strcmp(*p++, "--") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(*p == 0);

    assert(c == 'a');
    assert(option_index == 5);  /* no changed */
    assert(optind == 2);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-c") == 0);
    assert(strcmp(*p++, "c") == 0);
    assert(strcmp(*p++, "--") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(*p == 0);

    assert(c == 'c');
    assert(option_index == 5);  /* no changed */
    assert(optind == 4);
    assert(optarg == argv[3]);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-c") == 0);
    assert(strcmp(*p++, "c") == 0);
    assert(strcmp(*p++, "--") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(*p == 0);

    assert(c == -1);
    assert(option_index == 5);  /* no changed */
    assert(optind == 5);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-c") == 0);
    assert(strcmp(*p++, "c") == 0);
    assert(strcmp(*p++, "--") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(*p == 0);

    assert(c == 'd');
    assert(option_index == 5);  /* no changed */
    assert(optind == 7);
    assert(optarg == argv[6]);
    assert(optopt == 0);
}

void
test9()
{
    optind = 0;
    argc = 0;
    p = argv;

    argc++; *p++ = "command_name";
    argc++; *p++ = "-a";
    argc++; *p++ = "-";
    argc++; *p++ = "-";
    argc++; *p++ = "-c";
    argc++; *p++ = "c";
    argc++; *p++ = "--";
    argc++; *p++ = "-d";
    argc++; *p++ = "d";
    *p = 0;

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-") == 0);
    assert(strcmp(*p++, "-") == 0);
    assert(strcmp(*p++, "-c") == 0);
    assert(strcmp(*p++, "c") == 0);
    assert(strcmp(*p++, "--") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(*p == 0);

    assert(c == 'a');
    assert(option_index == 5);  /* no changed */
    assert(optind == 2);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-") == 0);
    assert(strcmp(*p++, "-") == 0);
    assert(strcmp(*p++, "-c") == 0);
    assert(strcmp(*p++, "c") == 0);
    assert(strcmp(*p++, "--") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(*p == 0);

    assert(c == 'c');
    assert(option_index == 5);  /* no changed */
    assert(optind == 6);
    assert(optarg == argv[5]);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-c") == 0);
    assert(strcmp(*p++, "c") == 0);
    assert(strcmp(*p++, "--") == 0);
    assert(strcmp(*p++, "-") == 0);
    assert(strcmp(*p++, "-") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(*p == 0);

    assert(c == -1);
    assert(option_index == 5);  /* no changed */
    assert(optind == 5);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc:d:", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-c") == 0);
    assert(strcmp(*p++, "c") == 0);
    assert(strcmp(*p++, "--") == 0);
    assert(strcmp(*p++, "-") == 0);
    assert(strcmp(*p++, "-") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(*p == 0);

    assert(c == 'd');
    assert(option_index == 5);  /* no changed */
    assert(optind == 9);
    assert(optarg == argv[8]);
    assert(optopt == 0);
}

void
test10()
{
    optind = 0;
    argc = 0;
    p = argv;

    argc++; *p++ = "command_name";
    argc++; *p++ = "-a";
    argc++; *p++ = "-cc";
    argc++; *p++ = "-d";
    argc++; *p++ = "d";
    argc++; *p++ = "-c";        /* no argument */
    argc++; *p++ = "-d";        /* at last */
    *p = 0;

    /*************************/
    c = call_getopt_long(argc, argv, "abc::d::", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-cc") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(strcmp(*p++, "-c") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(*p == 0);

    assert(c == 'a');
    assert(option_index == 5);  /* no changed */
    assert(optind == 2);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc::d::", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-cc") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(strcmp(*p++, "-c") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(*p == 0);

    assert(c == 'c');
    assert(option_index == 5);  /* no changed */
    assert(optind == 3);
    assert(optarg == argv[2]+2);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc::d::", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-cc") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(strcmp(*p++, "-c") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(*p == 0);

    assert(c == 'd');
    assert(option_index == 5);  /* no changed */
    assert(optind == 4);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc::d::", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-cc") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(strcmp(*p++, "-c") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(*p == 0);

    assert(c == 'c');
    assert(option_index == 5);  /* no changed */
    assert(optind == 6);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc::d::", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-cc") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(strcmp(*p++, "-c") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(*p == 0);

    assert(c == 'd');
    assert(option_index == 5);  /* no changed */
    assert(optind == 7);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc::d::", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-cc") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(strcmp(*p++, "-c") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(*p == 0);

    assert(c == -1);
    assert(option_index == 5);  /* no changed */
    assert(optind == 6);
    assert(optarg == 0);
    assert(optopt == 0);
}

void
test11()
{
    optind = 0;
    argc = 0;
    p = argv;

    argc++; *p++ = "command_name";
    argc++; *p++ = "--verbose";
    argc++; *p++ = "--create=c";
    argc++; *p++ = "--change";
    argc++; *p++ = "d";
    argc++; *p++ = "--create";  /* no argument */
    argc++; *p++ = "--change";  /* at last */
    *p = 0;

    /*************************/
    c = call_getopt_long(argc, argv, "abc::d::", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "--verbose") == 0);
    assert(strcmp(*p++, "--create=c") == 0);
    assert(strcmp(*p++, "--change") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(strcmp(*p++, "--create") == 0);
    assert(strcmp(*p++, "--change") == 0);
    assert(*p == 0);

    assert(c == 0);
    assert(option_index == 0);
    assert(optind == 2);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc::d::", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "--verbose") == 0);
    assert(strcmp(*p++, "--create=c") == 0);
    assert(strcmp(*p++, "--change") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(strcmp(*p++, "--create") == 0);
    assert(strcmp(*p++, "--change") == 0);
    assert(*p == 0);

    assert(c == 0);
    assert(option_index == 5);
    assert(optind == 3);
    assert(optarg == argv[2]+9);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc::d::", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "--verbose") == 0);
    assert(strcmp(*p++, "--create=c") == 0);
    assert(strcmp(*p++, "--change") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(strcmp(*p++, "--create") == 0);
    assert(strcmp(*p++, "--change") == 0);
    assert(*p == 0);

    assert(c == 0);
    assert(option_index == 6);
    assert(optind == 4);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc::d::", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "--verbose") == 0);
    assert(strcmp(*p++, "--create=c") == 0);
    assert(strcmp(*p++, "--change") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(strcmp(*p++, "--create") == 0);
    assert(strcmp(*p++, "--change") == 0);
    assert(*p == 0);

    assert(c == 0);
    assert(option_index == 5);
    assert(optind == 6);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc::d::", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "--verbose") == 0);
    assert(strcmp(*p++, "--create=c") == 0);
    assert(strcmp(*p++, "--change") == 0);
    assert(strcmp(*p++, "--create") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(strcmp(*p++, "--change") == 0);
    assert(*p == 0);

    assert(c == 0);
    assert(option_index == 6);
    assert(optind == 7);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc::d::", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "--verbose") == 0);
    assert(strcmp(*p++, "--create=c") == 0);
    assert(strcmp(*p++, "--change") == 0);
    assert(strcmp(*p++, "--create") == 0);
    assert(strcmp(*p++, "--change") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(*p == 0);

    assert(c == -1);
    assert(option_index == 6);
    assert(optind == 6);
    assert(optarg == 0);
    assert(optopt == 0);
}

void
test12()
{
    optind = 0;
    argc = 0;
    p = argv;

    argc++; *p++ = "command_name";
    argc++; *p++ = "--verbose";
    argc++; *p++ = "--create=c";
    argc++; *p++ = "files...";
    argc++; *p++ = "--delete";  /* required argument */
    argc++; *p++ = "d";
    argc++; *p++ = "--create";  /* no argument */
    argc++; *p++ = "--change";  /* at last */
    *p = 0;

    /*************************/
    c = call_getopt_long(argc, argv, "abc::d::", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "--verbose") == 0);
    assert(strcmp(*p++, "--create=c") == 0);
    assert(strcmp(*p++, "files...") == 0);
    assert(strcmp(*p++, "--delete") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(strcmp(*p++, "--create") == 0);
    assert(strcmp(*p++, "--change") == 0);
    assert(*p == 0);

    assert(c == 0);
    assert(option_index == 0);
    assert(optind == 2);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc::d::", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "--verbose") == 0);
    assert(strcmp(*p++, "--create=c") == 0);
    assert(strcmp(*p++, "files...") == 0);
    assert(strcmp(*p++, "--delete") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(strcmp(*p++, "--create") == 0);
    assert(strcmp(*p++, "--change") == 0);
    assert(*p == 0);

    assert(c == 0);
    assert(option_index == 5);
    assert(optind == 3);
    assert(optarg == argv[2]+9);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc::d::", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "--verbose") == 0);
    assert(strcmp(*p++, "--create=c") == 0);
    assert(strcmp(*p++, "files...") == 0);
    assert(strcmp(*p++, "--delete") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(strcmp(*p++, "--create") == 0);
    assert(strcmp(*p++, "--change") == 0);
    assert(*p == 0);

    assert(c == 0);
    assert(option_index == 4);
    assert(optind == 6);
    assert(optarg == argv[5]);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc::d::", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "--verbose") == 0);
    assert(strcmp(*p++, "--create=c") == 0);
    assert(strcmp(*p++, "--delete") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(strcmp(*p++, "files...") == 0);
    assert(strcmp(*p++, "--create") == 0);
    assert(strcmp(*p++, "--change") == 0);
    assert(*p == 0);

    assert(c == 0);
    assert(option_index == 5);
    assert(optind == 7);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc::d::", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "--verbose") == 0);
    assert(strcmp(*p++, "--create=c") == 0);
    assert(strcmp(*p++, "--delete") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(strcmp(*p++, "--create") == 0);
    assert(strcmp(*p++, "files...") == 0);
    assert(strcmp(*p++, "--change") == 0);
    assert(*p == 0);

    assert(c == 0);
    assert(option_index == 6);
    assert(optind == 8);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc::d::", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "--verbose") == 0);
    assert(strcmp(*p++, "--create=c") == 0);
    assert(strcmp(*p++, "--delete") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(strcmp(*p++, "--create") == 0);
    assert(strcmp(*p++, "--change") == 0);
    assert(strcmp(*p++, "files...") == 0);
    assert(*p == 0);

    assert(c == -1);
    assert(option_index == 6);
    assert(optind == 7);
    assert(optarg == 0);
    assert(optopt == 0);
}

void
test13()
{
    optind = 0;
    argc = 0;
    p = argv;

    argc++; *p++ = "command_name";
    argc++; *p++ = "--verbose";
    argc++; *p++ = "--create=c";
    argc++; *p++ = "files...";
    argc++; *p++ = "--delete";
    argc++; *p++ = "d";
    argc++; *p++ = "--";        /* option terminator */
    argc++; *p++ = "--change";
    *p = 0;

    /*************************/
    c = call_getopt_long(argc, argv, "abc::d::", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "--verbose") == 0);
    assert(strcmp(*p++, "--create=c") == 0);
    assert(strcmp(*p++, "files...") == 0);
    assert(strcmp(*p++, "--delete") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(strcmp(*p++, "--") == 0);
    assert(strcmp(*p++, "--change") == 0);
    assert(*p == 0);

    assert(c == 0);
    assert(option_index == 0);
    assert(optind == 2);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc::d::", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "--verbose") == 0);
    assert(strcmp(*p++, "--create=c") == 0);
    assert(strcmp(*p++, "files...") == 0);
    assert(strcmp(*p++, "--delete") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(strcmp(*p++, "--") == 0);
    assert(strcmp(*p++, "--change") == 0);
    assert(*p == 0);

    assert(c == 0);
    assert(option_index == 5);
    assert(optind == 3);
    assert(optarg == argv[2]+9);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc::d::", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "--verbose") == 0);
    assert(strcmp(*p++, "--create=c") == 0);
    assert(strcmp(*p++, "files...") == 0);
    assert(strcmp(*p++, "--delete") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(strcmp(*p++, "--") == 0);
    assert(strcmp(*p++, "--change") == 0);
    assert(*p == 0);

    assert(c == 0);
    assert(option_index == 4);
    assert(optind == 6);
    assert(optarg == argv[5]);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc::d::", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "--verbose") == 0);
    assert(strcmp(*p++, "--create=c") == 0);
    assert(strcmp(*p++, "--delete") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(strcmp(*p++, "--") == 0);
    assert(strcmp(*p++, "files...") == 0);
    assert(strcmp(*p++, "--change") == 0);
    assert(*p == 0);

    assert(c == -1);
    assert(option_index == 4);
    assert(optind == 6);
    assert(optarg == 0);
    assert(optopt == 0);

}

void
test14()
{
    optind = 0;
    argc = 0;
    p = argv;

    argc++; *p++ = "command_name";
    argc++; *p++ = "-o5";
    argc++; *p++ = "files...";
    *p = 0;

    /*************************/
    c = call_getopt_long(argc, argv, "o[567]", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-o5") == 0);
    assert(strcmp(*p++, "files...") == 0);
    assert(*p == 0);

    assert(c == 'o');
    assert(option_index == 4);  /* no changed */
    assert(optind == 2);
    assert(optarg == argv[1]+2);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc[cde]d[fgh]", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-o5") == 0);
    assert(strcmp(*p++, "files...") == 0);
    assert(*p == 0);

    assert(c == -1);
    assert(option_index == 4);  /* no changed */
    assert(optind == 2);
    assert(optarg == 0);
    assert(optopt == 0);

}

void
test15()
{
    optind = 0;
    argc = 0;
    p = argv;

    argc++; *p++ = "command_name";
    argc++; *p++ = "-a";
    argc++; *p++ = "-ccd";
    argc++; *p++ = "-ce";
    argc++; *p++ = "-d";
    argc++; *p++ = "d";
    argc++; *p++ = "-cdd";
    argc++; *p++ = "-d";
    *p = 0;

    /*************************/
    c = call_getopt_long(argc, argv, "abc[cde]d[fgh]", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-ccd") == 0);
    assert(strcmp(*p++, "-ce") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(strcmp(*p++, "-cdd") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(*p == 0);

    assert(c == 'a');
    assert(option_index == 4);  /* no changed */
    assert(optind == 2);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc[cde]d[fgh]", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-ccd") == 0);
    assert(strcmp(*p++, "-ce") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(strcmp(*p++, "-cdd") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(*p == 0);

    assert(c == 'c');
    assert(option_index == 4);  /* no changed */
    assert(optind == 2);
    assert(optarg == argv[2]+2);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc[cde]d[fgh]", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-ccd") == 0);
    assert(strcmp(*p++, "-ce") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(strcmp(*p++, "-cdd") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(*p == 0);

    assert(c == 'd');
    assert(option_index == 4);  /* no changed */
    assert(optind == 3);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc[cde]d[fgh]", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-ccd") == 0);
    assert(strcmp(*p++, "-ce") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(strcmp(*p++, "-cdd") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(*p == 0);

    assert(c == 'c');
    assert(option_index == 4);  /* no changed */
    assert(optind == 4);
    assert(optarg == argv[3]+2);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc[cde]d[fgh]", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-ccd") == 0);
    assert(strcmp(*p++, "-ce") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(strcmp(*p++, "-cdd") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(*p == 0);

    assert(c == 'd');
    assert(option_index == 4);  /* no changed */
    assert(optind == 5);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc[cde]d[fgh]", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-ccd") == 0);
    assert(strcmp(*p++, "-ce") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(strcmp(*p++, "-cdd") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(*p == 0);

    assert(c == 'c');
    assert(option_index == 4);  /* no changed */
    assert(optind == 6);
    assert(optarg == argv[6]+2);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc[cde]d[fgh]", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-ccd") == 0);
    assert(strcmp(*p++, "-ce") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(strcmp(*p++, "-cdd") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(*p == 0);

    assert(c == 'd');
    assert(option_index == 4);  /* no changed */
    assert(optind == 7);
    assert(optarg == 0);
    assert(optopt == 0);

    /*************************/
    c = call_getopt_long(argc, argv, "abc[cde]d[fgh]", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-ccd") == 0);
    assert(strcmp(*p++, "-ce") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(strcmp(*p++, "-cdd") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(*p == 0);

    assert(c == 'd');
    assert(option_index == 4);  /* no changed */
    assert(optind == 8);
    assert(optarg == 0);
    assert(optopt == 0);


    /*************************/
    c = call_getopt_long(argc, argv, "abc[cde]d[fgh]", long_options, &option_index);

    p = argv;
    assert(strcmp(*p++, "command_name") == 0);
    assert(strcmp(*p++, "-a") == 0);
    assert(strcmp(*p++, "-ccd") == 0);
    assert(strcmp(*p++, "-ce") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(strcmp(*p++, "-cdd") == 0);
    assert(strcmp(*p++, "-d") == 0);
    assert(strcmp(*p++, "d") == 0);
    assert(*p == 0);

    assert(c == -1);
    assert(option_index == 4);  /* no changed */
    assert(optind == 7);
    assert(optarg == 0);
    assert(optopt == 0);


}

int
main()
{
    opterr = 0;
    optopt = 0;

    test1();
    test2();
    test3();
    test4();
    test5();
    test6();
    test7();
    test8();
    test9();
    test10();
    test11();
    test12();
    test13();
#ifndef USE_GNU
    test14();
    test15();
#endif

    return 0;
}
#endif
