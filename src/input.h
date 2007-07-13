/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef INPUT_H
#define INPUT_H

#include "stdio.h"
#include "sysptr.h"

namespace SPPARKS {

class Input : protected SysPtr {
 public:
  int narg;                    // # of command args
  char **arg;                  // parsed args for command
  class Variable *variable;    // defined variables

  Input(class SPK *, int, char **);
  ~Input();
  void file();                   // process all input
  void file(char *);             // process an input script
  char *one(char *);             // process a single command
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
  void include();
  void jump();
  void label();
  void log();
  void next_command();
  void print();
  void variable_command();

  void app_style();
  void run();
  void solve_style();
  void sweep_style();
};

}

#endif
