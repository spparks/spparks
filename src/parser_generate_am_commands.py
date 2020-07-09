import sys
def count_args(params):
    nargs=0
    keys=list(params.keys())
    for k in keys:
        cmd=params[k]
        nargs+=(1+cmd['len'])
    return nargs

def write_path_parser():
    line={"am_path":{"type":"int","len":1},\
            "params":{"start":{"type":"double","len":2},\
            "end":{"type":"double","len":2},\
            "speed":{"type":"double","len":1}\
            } }
    cmd=list(line.keys())[0]
    params=line['params']
    # Initial '1' is for command 'id'
    nargs=1+count_args(params)
    s= "if (strcmp(command,\"am_path\") == 0) {\n"
    s+="   if (narg < "+str(nargs)+") error->all(FLERR,\"Illegal 'am_path' command; wrong num args.\");\n"
    s+="   int id=std::atoi(arg[0]);\n"
    s+="   double x0,y0,x1,y1;\n"
    s+="   double speed;\n"
    s+="   if(strcmp(arg[1],\"start\")==0){\n"
    s+="      x0=atof(arg[2]);\n"
    s+="      y0=atof(arg[3]);\n"
    s+="   }\n"
    s+="   else {error->all(FLERR,\"Illegal path command. Expected keyword 'start'\");}\n"
    s+="   if(strcmp(arg[4],\"end\")==0){\n"
    s+="      x1=atof(arg[5]);\n"
    s+="      y1=atof(arg[6]);\n"
    s+="   }\n"
    s+="   else {error->all(FLERR,\"Illegal path command. Expected keyword 'end'\");}\n"
    s+="   if(strcmp(arg[7],\"speed\")==0){\n"
    s+="      speed=atof(arg[8]);\n"
    s+="   }\n"
    s+="   else {error->all(FLERR,\"Illegal path command. Expected keyword 'speed'\");}\n"
    s+="   paths[id]=Path(Point(x0,y0),Point(x1,y1),speed);\n"
    s+="}\n"
    return s

def write_path_layer_parser():
    line={"am_path_layer":{"type":"int","len":1},\
            "params":{"num_paths":{"type":"int","len":1},\
                      "path_ids":{"type":"int","len":0},\
                      "thickness":{"type":"double","len":1}\
                     } \
        }
    cmd=list(line.keys())[0]
    params=line['params']
    s= "else if (strcmp(command,\"am_path_layer\") == 0) {\n"
    s+="   if (narg < 2) error->all(FLERR,\"Illegal 'am_path_layer' command; wrong num args.\");\n"
    s+="   int id=std::atoi(arg[0]);\n"
    s+="   int p;\n"
    s+="   if(strcmp(arg[1],\"num_paths\")==0){\n"
    s+="      int num_paths=std::atoi(arg[2]);\n"
    s+="      std::vector<int> path_ids(num_paths);\n"
    s+="      p=3;\n"
    s+="      for(int i(0);i<num_paths;i++){\n"
    s+="         path_ids[i]=std::atoi(arg[p]);\n"
    s+="         p+=1;\n"
    s+="      }\n"
    s+="   }\n"
    s+="   else {error->all(FLERR,\"Illegal path_layer command. Expected keyword 'num_paths'\");}\n"
    s+="   if(strcmp(arg[p],\"thickness\")==0){\n"
    s+="      int num_repeats=std::atoi(arg[p+1]);\n"
    s+="   }\n"
    s+="   else {error->all(FLERR,\"Illegal path_layer command. Expected keyword 'thickness'.\");}\n"
    s+="}\n"
    return s

def write_pattern_parser():
    line={"am_pattern":{"type":"int","len":1},\
            "params":{"num_layers":{"type":"int","len":1},\
                      "layer_ids":{"type":"int","len":0},\
                      "repeat":{"type":"int","len":1}\
                     } \
        }
    cmd=list(line.keys())[0]
    params=line['params']
    s= "else if (strcmp(command,\"am_pattern\") == 0) {\n"
    s+="   if (narg < 2) error->all(FLERR,\"Illegal 'am_pattern' command; wrong num args.\");\n"
    s+="   int id=std::atoi(arg[0]);\n"
    s+="   int p;\n"
    s+="   if(strcmp(arg[1],\"num_layers\")==0){\n"
    s+="      int num_layers=std::atoi(arg[2]);\n"
    s+="      if (narg != 2+num_layers+2) error->all(FLERR,\"Illegal 'am_pattern' command; wrong num args.\");\n"
    s+="      std::vector<int> layer_ids(num_layers);\n"
    s+="      p=3;\n"
    s+="      for(int i(0);i<num_layers;i++){\n"
    s+="         layer_ids[i]=std::atoi(arg[p]);\n"
    s+="         p+=1;\n"
    s+="      }\n"
    s+="   }\n"
    s+="   else {error->all(FLERR,\"Illegal pattern command. Expected keyword 'num_layers'\");}\n"
    s+="   if(strcmp(arg[p],\"repeat\")==0){\n"
    s+="      int num_repeats=std::atoi(arg[p+1]);\n"
    s+="   }\n"
    s+="   else {error->all(FLERR,\"Illegal pattern command. Expected keyword 'repeat'.\");}\n"
    s+="}\n"
    return s

def write_build_parser():
    line={"am_build":"", "params":{"z_start":{"type":"double","len":1} } }
    cmd=list(line.keys())[0]
    params=line['params']
    s= "else if (strcmp(command,\"am_build\") == 0) {\n"
    s+="   if (narg != 2) error->all(FLERR,\"Illegal 'am_build' command; wrong num args; should have 2 arguments.\");\n"
    s+="   double z_start;\n"
    s+="   if(strcmp(arg[0],\"z_start\")==0){\n"
    s+="      z_start=std::atoi(arg[1]);\n"
    s+="   } else {error->all(FLERR,\"Illegal am_build command. Expected keyword 'z_start'\");}\n"
    s+="   build_layer_z=z_start;\n"
    s+="}\n"
    return s

def write_parser():
    s=write_path_parser();
    s+=write_path_layer_parser()
    #s+=write_pattern_parser()
    s+=write_build_parser()
    print(s)

if __name__ == '__main__':
    sys.exit(write_parser())
