/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#ifndef SPK_INPUT_H
#define SPK_INPUT_H

#include "stdio.h"
#include "pointers.h"

namespace SPPARKS_NS {

class Input : protected Pointers {
 public:
  int narg;                    // # of command args
  char **arg;                  // parsed args for command
  class Variable *variable;    // defined variables

  Input(class SPPARKS *, int, char **);
  ~Input();
  void file();                   // process all input
  void file(const char *);       // process an input script
  char *one(const char *);       // process a single command
  void substitute(char *, int);  // substitute for variables in a string

 private:
  int me;                      // proc ID
  char *command;               // ptr to current command
  int maxarg;                  // max # of args in arg
  char *line,*copy,*work;      // input line & copy of it
  int echo_screen;             // 0 = no, 1 = yes
  int echo_log;                // 0 = no, 1 = yes
  int nfile,maxfile;           // current # and max # of open input files
  int label_active;            // 0 = no label, 1 = looking for label
  char *labelstr;              // label string being looked for
  int jump_skip;               // 1 if skipping next jump, 0 otherwise

  FILE **infiles;              // list of open input files

  void parse();                // parse an input text line
  int execute_command();       // execute a single command

  void clear();                // individual commands
  void echo();
  void ifthenelse();
  void include();
  void jump();
  void label();
  void log();
  void next_command();
  void print();
  void variable_command();

  void app_style();
  void boundary();
  void diag_style();
  void dimension();
  void dump();
  void dump_modify();
  void dump_one();
  void lattice();
  void pair_coeff();
  void pair_style();
  void processors();
  void region();
  void reset_time();
  void run();
  void seed();
  void solve_style();
  void stats();
  void undump();
};

}

#endif

/* ERROR/WARNING messages:

E: Label wasn't found in input script

Self-explanatory.

E: Input line too long: %s

This is a hard (very large) limit defined in the input.cpp file.

E: Unknown command: %s

The command is not known to SPPARKS.  Check the input script.

E: Another input script is already being processed

Cannot attempt to open a 2nd input script, when the original file is
still being processed.

E: Cannot open input script %s

Self-explanatory.

E: Unbalanced quotes in input line

No matching end double quote was found following a leading double
quote.

E: Invalid variable name

Variable name used in an input script line is invalid.

E: Substitution for illegal variable

Self-explanatory.

E: Input line too long after variable substitution

This is a hard (very large) limit defined in the input.cpp file.

E: App_style specific command before app_style set

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Cannot open logfile %s

Self-explanatory.

E: App_style command after simulation box is defined

Self-explanatory.

E: Boundary command after simulation box is defined

UNDOCUMENTED

E: Diag_style command before app_style set

Self-explanatory.

E: Dimension command after simulation box is defined

Self-explanatory.

E: Dimension command after lattice is defined

Self-explanatory.

E: Dump command before app_style set

Self-explanatory.

E: Dump_modify command before app_style set

Self-explanatory.

E: Dump_one command before app_style set

Self-explanatory.

E: Lattice command before app_style set

Self-explanatory.

E: Pair_coeff command before app_style set

Self-explanatory.

E: Pair_coeff command before pair_style is defined

Self-explanatory.

E: Pair_style command before app_style set

Self-explanatory.

E: Processors command after simulation box is defined

Self-explanatory.

E: Region command before app_style set

Self-explanatory.

E: Reset_time command before app_style set

Self-explanatory.

E: Run command before app_style set

Self-explanatory.

E: Solve_style command before app_style set

Self-explanatory.

E: Stats command before app_style set

Self-explanatory.

E: Undump command before app_style set

Self-explanatory.

*/
