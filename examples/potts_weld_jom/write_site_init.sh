# You must run 'potts_init.in' before running this script.

write_site_init(){
  # If only 1 time step is written to the dump file by 'spparks', 
  # then the following single command will skip the header portion 
  # of the dump file and then take the site id and spin from each line 
  # and then append 0.0 to the line and print.  Redirect the output 
  # to a file of your choosing.  In example, the spparks generated 
  # dump file (with only 1 time step) is command argument $1 and the 
  # the output is redirected to potts_init.init

  input_file=$1
  output_file=$2
  printf "# This line is ignored\n" >$output_file
  printf "Values\n\n" >> $output_file
  tail -n +10 $input_file | cut -d' ' -f1,2 | sed -e 's/$/ 0.0/' >> $output_file
}

SITE_INIT_FILE=site.init

# Translate dump file from to site initialization file in preparation for 
#   running weld model; see function 'write_site_init' above in this file;
write_site_init potts_init.dump $SITE_INIT_FILE

