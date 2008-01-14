/* ----------------------------------------------------------------------

------------------------------------------------------------------------- */
#include <iostream>
#include "random_park.h"
#include "node.h"
#include "tree.h"
#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include "plus_node.h"
#include "minus_node.h"
#include "times_node.h"
#include "divide_node.h"
#include "pow_node.h"
#include "const_node.h"
#include "var_node.h"
#include "dbl_var_node.h"
#include "int_var_node.h"
/* ---------------------------------------------------------------------- */
using namespace std;

using namespace SPPARKS;
Tree::Tree()
{
}
/* ---------------------------------------------------------------------- */
Tree::~Tree()
{
  delete [] parents;
  delete [] children;

  for (int c = 0; c < 50; c++) delete [] tokens[c];
  delete [] tokens;
  for (int c = 0; c < 50; c++) delete [] postfix[c];
  delete [] postfix;
  for (int c = 0; c < 50; c++) delete [] op_stack[c];
  delete [] op_stack;

  delete [] input_line;
}
/* ---------------------------------------------------------------------- */
void Tree::init(double const_lo_in, double const_hi_in, 
 int num_var, int max_depth_in)
{

  const_hi = const_hi_in;
  const_lo = const_lo_in;
  range = const_hi - const_lo;
  num_variables = num_var;
  max_depth = max_depth_in;

  nch = static_cast<int>
    (pow(2.0,static_cast<double>(max_depth)));
  npr = static_cast<int>
    (pow(2.0,static_cast<double>(max_depth-1)));

  children = new Node *[nch];
  parents = new Node *[npr];

  tokens = new char*[50];
  postfix = new char*[50];
  op_stack = new char*[50];
  for (int c = 0; c < 50; c++) {
    tokens[c] = new char[16];
    postfix[c] = new char[16];
    op_stack[c] = new char[16];
  }

  clear_tokens(tokens);
  clear_tokens(postfix);
  clear_tokens(op_stack);

  input_line = new char[250];
  state = NULL;
}
/* ---------------------------------------------------------------------- */
Node *Tree::random_node(RandomPark* random, int ran_type)
{
  Node *rn;
  int op;
  double c;
  int v;

  //type 0 select all, type 1 select leaf
  if(ran_type == 0) op = static_cast<int>((VARIABLE+1)*random->uniform());
  else if (ran_type == 1) 
    op = CONSTANT + static_cast<int>(2*random->uniform());
  if(op == CONSTANT)   {
                           rn = new ConstNode();
    c = const_lo + range * random->uniform();
    static_cast<ConstNode*> (rn)->set_value(c);
  }
  else if(op == VARIABLE)   {
                           rn = new VarNode();
    v = static_cast<int>(num_variables * random->uniform());
    static_cast<VarNode*> (rn)->set_value(v);
  }
  else if(op == PLUS)      rn = new PlusNode();
  else if(op == MINUS)     rn = new MinusNode();
  else if(op == TIMES)     rn = new TimesNode();
  else if(op == DIVIDE)    rn = new DivideNode();
  else if(op == POW)       rn = new PowNode();
    
  return rn;
}
/* ---------------------------------------------------------------------- */
Node *Tree::build_tree(RandomPark* random, int max_depth)
{
  int depth = 0;
  bool terminated = false;
  int t;
  int finish = 0;
  int ccnt, pcnt;

  for(t=0;t<npr;t++) parents[t] = NULL;
  for(t=0;t<nch;t++) children[t] = NULL;

  root = random_node(random,0);
  parents[0] = root;
  pcnt = 1;

  while (depth < max_depth & !terminated){
    terminated = true;
    if (depth >= max_depth - 2) finish = 1;
    //find available terminals and make children
    ccnt = 0;
    for(t=0;t<pcnt;t++)
      if(parents[t]->type < CONSTANT){
	terminated = false;
	children[ccnt] = random_node(random,finish);
	children[ccnt+1] = random_node(random,finish);
	parents[t]->set_left_child(children[ccnt]);
	parents[t]->set_right_child(children[ccnt+1]);
	ccnt += 2;
      }
    //display parents
//     cout << "depth "<<depth<<" : ";
//     for(t=0;t<pcnt;t++) parents[t]->write_stack(screen);
//     cout << endl;

    //promote children
    for(t=0;t<ccnt;t++) {
      parents[t] = children[t];
      children[t] = NULL;
    }
    pcnt = ccnt;
    //increment depth
    depth++;
  }
  return root;
}
/* ---------------------------------------------------------------------- */
void Tree::mutate_branch(RandomPark* random, Node *&root)
{
  int dir;
  int depth;
  int ran_depth;
  Node *current_node, *current_parent;

  current_node = root;
  current_parent = NULL;

  //descend to random depth
  ran_depth = static_cast<int>(max_depth*random->uniform());
  depth = 0;

  while (depth < ran_depth){
    if(current_node->type < CONSTANT){ 
      current_parent = current_node;
      dir = static_cast<int>(random->uniform()*2);
      if(dir == 0)
	current_node = static_cast<PlusNode*>
	  (current_node)->get_left_child();
      else if(dir == 1)
	current_node = static_cast<PlusNode*>
	  (current_node)->get_right_child();
      depth ++;
    }
    else break;
  }

  //create tree of remaining depth

  depth = max_depth - depth;
  current_node->clear();
  current_node = build_tree(random,depth);
  //  current_node->write();

  if(current_parent){
    if(dir == 0) current_parent->set_left_child(current_node);
    else if(dir == 1) current_parent->set_right_child(current_node);
  }
  else root = current_node;
}
/* ---------------------------------------------------------------------- */

void Tree::mutate_constant(RandomPark* random, Node *&root, 
			   double temperature)
{
  Node *leaves[nch];
  bool done = false;
  Node *current_node;
  int npr, nch, nlv;
  int depth, pos;

  if(root->type == VARIABLE) return;
  double new_value = random->uniform();
  if(root->type == CONSTANT){
    new_value = static_cast<ConstNode*>(root)->get_value() - const_lo + 
      temperature * new_value;
    pos = static_cast<int>(new_value / range);
    new_value -= pos*range;
    new_value += const_lo;
    static_cast<ConstNode*>(root)->set_value(new_value); 
    return;
  }
  current_node = root;
  parents[0] = root;
  npr = 1;
  nch = 0;
  nlv = 0;
  depth = 0;

  //find leaves

  while (!done & depth < max_depth){
    done = true;
    nch = 0;
    //find children and record constants
    for (int n = 0; n < npr; n++){
      current_node = static_cast<PlusNode*>
	(parents[n])->get_left_child();
      if(current_node->type < POW+1)
	{children[nch] = current_node; nch++; done = false;}
      if (current_node->type == CONSTANT) {leaves[nlv]=current_node;nlv++;}

      current_node = static_cast<PlusNode*>
	(parents[n])->get_right_child();
      if(current_node->type<4)
	{children[nch] = current_node; nch++; done = false;}
      if (current_node->type == CONSTANT) {leaves[nlv]=current_node;nlv++;}
    }
    //promote children
    for (int n = 0; n < nch; n++) parents[n] = children[n];
    npr = nch;
    depth ++;
  }

  if (nlv > 0){
    int leaf = static_cast<int>(nlv * random->uniform());
    new_value = static_cast<ConstNode*>(leaves[leaf])->get_value() 
      - const_lo + temperature * new_value;
    pos = static_cast<int>(new_value / range);
    new_value -= pos*range;
    new_value += const_lo;
    //     cout << "Changed constant "<< leaf << " from "; 
    //     leaves[leaf]->write_stack();
    //     cout  << " to "<< new_value<<endl;
    static_cast<ConstNode*>(leaves[leaf])->set_value(new_value);
  }
}
/* ---------------------------------------------------------------------- */
void Tree::mutate_variable(RandomPark* random, Node *&root)
{
  Node *leaves[nch];
  bool done = false;
  Node *current_node;
  int npr, nch, nlv;
  int depth;

  if(root->type == CONSTANT || num_variables == 1) return;
  int v = static_cast<int>(num_variables * random->uniform());
  if(root->type == VARIABLE) 
    {static_cast<VarNode*>(root)->set_value(v);; return;}

  current_node = root;
  parents[0] = root;
  npr = 1;
  nch = 0;
  nlv = 0;
  depth = 0;

  //find leaves

  while (!done & depth < max_depth){
    done = true;
    nch = 0;
    //find children and record constants
    for (int n = 0; n < npr; n++){
      current_node = static_cast<PlusNode*>
	(parents[n])->get_left_child();
      if(current_node->type<POW+1)
	{children[nch] = current_node; nch++; done = false;}
      if (current_node->type == VARIABLE) {leaves[nlv]=current_node;nlv++;}

      current_node = static_cast<PlusNode*>
	(parents[n])->get_right_child();
      if(current_node->type<POW+1)
	{children[nch] = current_node; nch++; done = false;}
      if (current_node->type == VARIABLE) {leaves[nlv]=current_node;nlv++;}
    }
    //promote children
    for (int n = 0; n < nch; n++) parents[n] = children[n];
    npr = nch;
    depth ++;
  }

  if(nlv > 0){
    
    int leaf = static_cast<int>(nlv * random->uniform());

//     cout << "Changed variable "<< leaf << " from "; 
//     leaves[leaf]->write_stack();
//     cout  << " to ";
    
    static_cast<VarNode*>(leaves[leaf])->set_value(v);
    
//     leaves[leaf]->write_stack();
//     cout  <<endl;
  }
}
/* ---------------------------------------------------------------------- */
void Tree::mutate_operator(RandomPark* random, Node *&root)
{
  Node *leaves[nch];
  Node *daddy[nch];
  int dir[nch];
  bool done = false;
  Node *current_node;
  int npr, nch, nlv;
  int depth;

  if(root->type == CONSTANT || root->type == VARIABLE) return;

  current_node = root;
  parents[0] = root;
  npr = 1;
  nch = 0;
  nlv = 0;
  depth = 0;

  //find leaves

  leaves[nlv]=current_node; 
  daddy[nlv] = NULL;
  dir[nlv] = -1; 
  nlv++;

  while (!done & depth < max_depth){
    done = true;
    nch = 0;
    //find children and record constants
    for (int n = 0; n < npr; n++){
      current_node = static_cast<PlusNode*>
	(parents[n])->get_left_child();
      if(current_node->type<POW+1){
	children[nch] = current_node; nch++; 
	done = false;
	leaves[nlv] = current_node;
	daddy[nlv] = parents[n];
	dir[nlv] = 0; 
	nlv++;
      }
      current_node = static_cast<PlusNode*>
	(parents[n])->get_right_child();
      if(current_node->type<POW+1){
	children[nch] = current_node; nch++; 
	done = false;
	leaves[nlv]=current_node; 
	daddy[nlv] = parents[n];
	dir[nlv] = 1; 
	nlv++;
      }
    }
    //promote children
    for (int n = 0; n < nch; n++) parents[n] = children[n];
    npr = nch;
    depth ++;
  }

  if(nlv > 0){
    Node *rn;
    
    int leaf = static_cast<int>(nlv * random->uniform());
    int op = static_cast<int>((POW+1) * random->uniform());

    if     (op == PLUS)      rn = new PlusNode();
    else if(op == MINUS)     rn = new MinusNode();
    else if(op == TIMES)     rn = new TimesNode();
    else if(op == DIVIDE)    rn = new DivideNode();
    else if(op == POW)       rn = new PowNode();

    Node *left = static_cast<PlusNode*>(leaves[leaf])->get_left_child();
    rn->set_left_child(left);
    Node *right = static_cast<PlusNode*>(leaves[leaf])->get_right_child();
    rn->set_right_child(right);

    if(dir[leaf] == 0)      daddy[leaf]->set_left_child(rn);
    else if(dir[leaf] == 1) daddy[leaf]->set_right_child(rn);
    else if(dir[leaf] == -1) root = rn;

    delete leaves[leaf];
    
//     cout << "Changed operator " << leaf << " from "; 
//     leaves[leaf]->write_stack(screen);
//     cout  << " to ";
    
//     rn->write_stack(screen);
//     cout  <<endl;
  }
}
/* ---------------------------------------------------------------------- */
void Tree::crossover(RandomPark* random, Node *&root1, Node *&root2, int ran_depth)
{
  int dir1, dir2;
  int depth;
  Node *current_node, *current_parent;
  Node * node1, * node2;
  Node * daddy1, * daddy2;
  
  //descend to given depth of tree 1
  current_node = root1;
  current_parent = NULL;
  depth = 0;

  while (depth < ran_depth){
    if(current_node->type<POW+1){ 
      current_parent = current_node;
      dir1 = static_cast<int>(random->uniform()*2);
      if(dir1 == 0)
	current_node = static_cast<PlusNode*>
	  (current_node)->get_left_child();
      else
	current_node = static_cast<PlusNode*>
	  (current_node)->get_right_child();
      depth ++;
    }
    else break;
  }

  node1 = current_node;
  daddy1 = current_parent;

  //descend to same depth of tree 2
  current_node = root2;
  current_parent = NULL;
  depth = 0;

  while (depth < ran_depth){
    if(current_node->type<POW+1){ 
      current_parent = current_node;
      dir2 = static_cast<int>(random->uniform()*2);
      if(dir2 == 0)
	current_node = static_cast<PlusNode*>
	  (current_node)->get_left_child();
      else
	current_node = static_cast<PlusNode*>
	  (current_node)->get_right_child();
      depth ++;
    }
    else break;
  }

  node2 = current_node;
  daddy2 = current_parent;

  //swap the trees
  if(daddy1){
    if(dir1 == 0) daddy1->set_left_child(node2);
    else daddy1->set_right_child(node2);
  }
  else root1 = node2;
  if(daddy2){
    if(dir2 == 0) daddy2->set_left_child(node1);
    else daddy2->set_right_child(node1);
  }
  else root2 = node1;

}
/* ---------------------------------------------------------------------- */
void Tree::swap(Node *&root1, Node *&root2)
{
  Node *temp;

  temp = root1;
  root1 = root2;
  root2 = temp;
}
/* ---------------------------------------------------------------------- */
Node * Tree::copy(Node *root)
{
  int op = root->type;
  Node *rn;
  Node *nt;
  int type;
  char name[16];
  void *var;
  
  if     (op == PLUS)      rn = new PlusNode();
  else if(op == MINUS)     rn = new MinusNode();
  else if(op == TIMES)     rn = new TimesNode();
  else if(op == DIVIDE)    rn = new DivideNode();
  else if(op == POW)       rn = new PowNode();
  else if(op == CONSTANT)  rn = new ConstNode();
  else if(op == VARIABLE)  rn = new VarNode();
  else if(op == VAR_DBL)  rn = new DblVarNode();
  else if(op == VAR_INT)  rn = new IntVarNode();

  if (op < POW+1) {
    nt = copy(static_cast<PlusNode*>(root)->get_left_child());
    rn->set_left_child(nt);
    nt = copy(static_cast<PlusNode*>(root)->get_right_child());
    rn->set_right_child(nt);
  }
  else if (op == CONSTANT) static_cast<ConstNode*>(rn)->
      set_value(static_cast<ConstNode*>(root)->get_value());
  else if (op == VARIABLE) static_cast<VarNode*>(rn)->
      set_value(static_cast<VarNode*>(root)->get_value());

  else if(op == VAR_DBL || op == VAR_INT)  {
    strcpy(name,static_cast<VarNode*>(root)->get_name());
    static_cast<VarNode*>(rn)->set_name(name);
    var = state->get_attribute(name, type);
    if(op == VAR_DBL)
      static_cast<DblVarNode*>(rn)->set_data_pointer(var);
    if(op == VAR_INT)
      static_cast<IntVarNode*>(rn)->set_data_pointer(var);
    static_cast<VarNode*>(rn)->set_data_type(type);
  }
  return rn;
}
/* ---------------------------------------------------------------------- */
Node *Tree::from_buffer(char **buf, int &n)
{
  int op = 0;
  double c = 0.0;
  int index = 0;
  Node *rn;
  Node *nt;
  int type;

  op = atoi(buf[n]);

  if     (op == PLUS)      rn = new PlusNode();
  else if(op == MINUS)     rn = new MinusNode();
  else if(op == TIMES)     rn = new TimesNode();
  else if(op == DIVIDE)    rn = new DivideNode();
  else if(op == POW)       rn = new PowNode();
  else if(op == CONSTANT){ 
    rn = new ConstNode(); 
    sscanf(buf[n],"%d %lg ",&op, &c);
  }
  else if(op == VARIABLE){  
    rn = new VarNode();
    sscanf(buf[n],"%d %d ",&op, &index);
  }
  else if(op == VAR_DBL)  {
    rn = new DblVarNode();
  }
  else if(op == VAR_INT)  {
    rn = new IntVarNode();
  }

 if (op < POW+1) {
    n++; nt = from_buffer(buf, n);
    rn->set_left_child(nt);
    n++; nt = from_buffer(buf, n);
    rn->set_right_child(nt);
  }
  else if (op == CONSTANT) {static_cast<ConstNode*>(rn)->set_value(c);}
  else if (op == VARIABLE) {static_cast<VarNode*>(rn)->set_value(index);}
  else if(op == VAR_DBL)  {

  }
  else if(op == VAR_INT)  {

  }
  return rn;
}
/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
Node *Tree::from_string(char *str)
{
  Node *rn = NULL;

  //  clean_string(str);
  strcpy(input_line,str);
  //  strcat(input_line,"\0");
  //convert input formula string to tokens
  string2tokens();
  //convert input tokens to postfix
  infix2postfix();
  //convert postfix to a tree
  rn = postfix2tree();

  return rn;
}
/* ---------------------------------------------------------------------- */
void Tree::clean_string(char *str)
{
  int i = 0;
  int j = 0;

  cout << "length of string incoming to clean = "<<strlen(str)<<endl;
  cout <<"the characters are: "<<endl;
  while(i != strlen(str)){
    cout <<i<<"  "<<str[i]<<endl;
      i++;
  }
  while (i != strlen(str))
    {
      if (str[i] == 32)
        {
	  i++;
        }
      else
        {
	  input_line[j] = str[i];
	  i++;
	  j++;
        }
    }
  strcat(input_line,"\0");
  input_line[j] = 0;

}
/* ---------------------------------------------------------------------- */
void Tree::infix2postfix()
{
  int k = 0;
  int n = 0;
  int m = 0;
  int t = -1;
  bool flag = true;

  char current[250];
  memset(current,0,250);

  clear_tokens(op_stack);
  clear_tokens(postfix);


  t = op_priority(tokens[0]);
 
  while(t != 99){

    t = op_priority(tokens[k]);

    if(t == 0){
      flag = push(tokens[k],op_stack);}

    else if(t == 1){
      while (!stack_empty(op_stack) && t < stack_peek(op_stack))
	{
	  flag = pop(current, op_stack);
	  flag = push(current,postfix);
	}
      flag = pop(current,op_stack);
      if(strcmp(current,"(")>0) 
	cout <<"paren pair failed to anihilate!"<<endl;
    }

    else if (t == 2 || t == 3 || t == 4){
      while (!stack_empty(op_stack) & stack_peek(op_stack)>=t)
	{
	  flag = pop(current, op_stack);
	  flag = push(current,postfix);
	}
      flag = push(tokens[k],op_stack);
    }

    else if (t == 8 || t == 10){
      flag = push(tokens[k],postfix);}


    else flag = false;

//     cout << "token "<<k<<"  "<<tokens[k]<<endl;
//     cout <<"op_stack: "<<endl;
//     if(!stack_empty(op_stack)) write_tokens(op_stack);
//     else cout <<"empty."<<endl;
//     //    cout <<"peek priority = "<<stack_peek(op_stack)<<endl;
//     cout <<"postfix: "<<endl;
//     if(!stack_empty(postfix)) write_tokens(postfix);
//     else cout <<"empty."<<endl;
//     //    cout <<"peek priority = "<<stack_peek(postfix)<<endl;

    k++;
  }
  while (!stack_empty(op_stack))
    {
      flag = pop(current, op_stack);
      flag = push(current,postfix);
    }
//   cout <<"op_stack: "<<endl;
//   if(!stack_empty(op_stack)) write_tokens(op_stack);
//   else cout <<"empty."<<endl;
//   cout <<"peek priority = "<<stack_peek(op_stack)<<endl;
//   cout <<"postfix: "<<endl;
//   if(!stack_empty(postfix)) write_tokens(postfix);
//   else cout <<"empty."<<endl;
//   cout <<"peek priority = "<<stack_peek(postfix)<<endl;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
void Tree::string2tokens()
{
  int m = 0;
  char *temp;

  // cout  <<"input_line "<<input_line<<" of length "<<strlen(input_line)<<endl;
  strcat(input_line,"#\0");

  char test[250];
  memset(test,0,250);

  int pos1 = 0;
  int pos2 = 0;

  while(pos1 < strlen(input_line)){
    if(strncmp(&input_line[pos1],"#",1)==0 ||
       strncmp(&input_line[pos1],"+",1)==0 ||
       strncmp(&input_line[pos1],"-",1)==0 ||
       strncmp(&input_line[pos1],"*",1)==0 ||
       strncmp(&input_line[pos1],"/",1)==0 ||
       strncmp(&input_line[pos1],"^",1)==0 ||
       strncmp(&input_line[pos1],"(",1)==0 ||
       strncmp(&input_line[pos1],")",1)==0
       ) 
      {
	strncpy(&test[pos2]," ",1); pos2++;
	test[pos2] = input_line[pos1]; pos2++;
	strncpy(&test[pos2]," ",1); pos2++;
	pos1++;
      }
    else{
      test[pos2] = input_line[pos1]; 
      pos1++;pos2++;
    }
  }
  while(test[pos2]==32) pos2--; //remove trailing zeroes
  test[pos2+1] = 0; //terminate string

  clear_tokens(tokens);

  char* tkn;

  tkn = strtok(test, " ");
  while(tkn != NULL){
    //    printf ("%s\n",tkn);
    strcpy(tokens[m],tkn);
    strcat(tokens[m],"\0");
    //    printf ("token %d %s\n",m, tokens[m]);
    tkn = strtok(NULL, " ");
    m++;
  }
  strcpy(tokens[m],"#");

  cout <<"tokens: "<<endl;
  if(!stack_empty(tokens)) write_tokens(tokens);
  else cout <<"empty."<<endl;
  //  cout <<"peek priority = "<<stack_peek(postfix)<<endl;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
Node *Tree::postfix2tree()
{
  vector<Node *> node_stack;
  int k = 0;
  int op = 0;
  string current_string;
  string var_name;
  char var_name_char[16];
  Node *current_nd;
  Node *tnd;

  bool flag = true;
  double tmp = 0.0;
  void *var_p = NULL;
  int type = -1;

  while(op_priority(postfix[k])<99){
    current_string.clear();
    current_string.assign(postfix[k]);

    op = op_priority(postfix[k]);

    if(current_string.compare(0,1,"+")==0)
      current_nd = new PlusNode();
    else if (current_string.compare(0,1,"-")==0)
      current_nd = new MinusNode();
    else if (current_string.compare(0,1,"*")==0)
      current_nd = new TimesNode();
    else if (current_string.compare(0,1,"/")==0)
      current_nd = new DivideNode();
    else if (current_string.compare(0,1,"^")==0)
      current_nd = new PowNode();
    else if (current_string.compare(0,1,"@")==0)
      {
	if(state==NULL){
	  cout <<"state is not initialized in tree."<<endl;
	  break;
	}
	var_name.clear();
	var_name.assign(postfix[k]);
	var_name.erase(0,1);
	strcpy(var_name_char,var_name.c_str());
	
	void *data = NULL;
	data = state->get_attribute(var_name_char, type);
	if(data==NULL) break;
	
	//	if (current_string.compare(1,3,"dbl")==0){
	if (type==2){
	  current_nd = new DblVarNode();
	  static_cast<DblVarNode*>(current_nd)->set_data_pointer(data);
	  //	  if(type!=2)cout <<"type mismatch in name desc."<<endl;
	}
	//	else if (current_string.compare(1,3,"int")==0){
	else if (type==1){
	  current_nd = new IntVarNode();
	  static_cast<IntVarNode*>(current_nd)->set_data_pointer(data);
	  //	  if(type!=1)cout <<"type mismatch in name desc."<<endl;
	}
	
	static_cast<VarNode*>(current_nd)->set_data_type(type);
	static_cast<VarNode*>(current_nd)->set_name(var_name_char);
      }
    
    else if (current_string.find_first_of("0123456789")<1){
      current_nd = new ConstNode();
      sscanf(postfix[k],"%lg",&tmp);
      static_cast<ConstNode*>(current_nd)->set_value(tmp);
    }
    else
      {
	cout << "Invalid token in formula "<<postfix[k]<<endl;
	break;
      }
    if(op<8){
      if(!node_stack.empty()){
	tnd = copy(node_stack.back());
	current_nd->set_right_child(tnd);
	node_stack.pop_back();
      }
      else cout <<"postfix error: no operands"<<endl;
      if(!node_stack.empty()){
	tnd = copy(node_stack.back());
	current_nd->set_left_child(tnd);
	node_stack.pop_back();
      }
      else cout <<"postfix error: no operands"<<endl;
    }
    node_stack.push_back(current_nd);
    k++;
  }
  if(!node_stack.empty()) return node_stack.back();
  else {cout << "postfix error. No root remains."<<endl; return NULL;}


}
/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
int Tree::op_priority(char *op)
{

  switch ( op[0] ) {
    
  case '(' :  
    return 0;
  case ')' :
    return 1;
  case '+' :  
    return 2;
  case '-' :
    return 2;
  case '*' :  
    return 3;
  case '/' :
    return 3;
  case '^' :
    return 4;
  case '@' :
    return 8;
  case '#' :
    return 99;
  default:
    return 10;
  }

}
/* ---------------------------------------------------------------------- */
void Tree::clear_tokens(char **tk)
{
  memset(tk[0],'#',1);
}
/* ---------------------------------------------------------------------- */
bool Tree::pop(char last[250], char **tk)
{
  int t = 0;

  while (op_priority(tk[t])<99) t++;

  if(t>0){
    strcpy(last, tk[t-1]);
    strcat(last,"\0");
    memset(tk[t-1],'#',1);
    return true;
  }
  return false;
}
/* ---------------------------------------------------------------------- */
bool Tree::push(char last[250], char **tk)
{
  int t = 0;
  
  while (op_priority(tk[t])<99) t++;

  if(t<50){
    strcpy(tk[t],last);
    strcat(tk[t],"\0");
    memset(tk[t+1],'#',1);
    return true;
  }
  return false;
}
/* ---------------------------------------------------------------------- */
void Tree::write_tokens(char **tk)
{
  int t = 0;
  
  while (op_priority(tk[t])<99) t++;

  if (t<50 & t> 0)
    for(int i = 0; i < t; i++)
      cout << "word "<<i<<"  "<<tk[i]<<endl;
  else cout <<"exceeded stack size with no terminating character"<<endl;
}
/* ---------------------------------------------------------------------- */
bool Tree::stack_empty(char **tk)
{
  int t = 0;
  while (op_priority(tk[t])<99) t++;

  if(t==0){
    return true;
  }
  return false;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int Tree::stack_peek(char **tk)
{
  int t = 0;
  while (op_priority(tk[t])<99) t++;

  if(t>0){
    return op_priority(tk[t-1]);
  }
  return 99;
}
/* ---------------------------------------------------------------------- */
