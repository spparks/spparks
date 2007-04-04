/* ----------------------------------------------------------------------

------------------------------------------------------------------------- */
#include <iostream>
#include "random_park.h"
#include "node.h"
#include "tree.h"
#include <cmath>
#include "string.h"
#include "plus_node.h"
#include "minus_node.h"
#include "times_node.h"
#include "divide_node.h"
#include "pow_node.h"
#include "const_node.h"
#include "var_node.h"

/* ---------------------------------------------------------------------- */
using namespace std;
using namespace SPPARKS;

Tree::Tree() : SysPtr(spk)
{
}
/* ---------------------------------------------------------------------- */
Tree::~Tree()
{
  delete [] parents;
  delete [] children;
  delete random;
}
/* ---------------------------------------------------------------------- */
void Tree::init(int seed, double const_lo_in, double const_hi_in, 
 int num_var, int max_depth_in)
{

  if (random != NULL)delete random;
  random = new RandomPark(seed);

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

}
/* ---------------------------------------------------------------------- */
Node *Tree::random_node(int ran_type)
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
Node *Tree::build_tree(int max_depth)
{
  int depth = 0;
  bool terminated = false;
  int t;
  int finish = 0;
  int ccnt, pcnt;

  for(t=0;t<npr;t++) parents[t] = NULL;
  for(t=0;t<nch;t++) children[t] = NULL;

  root = random_node(0);
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
	children[ccnt] = random_node(finish);
	children[ccnt+1] = random_node(finish);
	parents[t]->set_left_child(children[ccnt]);
	parents[t]->set_right_child(children[ccnt+1]);
	ccnt += 2;
      }
    //display parents
//     cout << "depth "<<depth<<" : ";
//     for(t=0;t<pcnt;t++) parents[t]->write_stack();
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
void Tree::mutate_branch(Node *&root)
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
  current_node = build_tree(depth);
  //  current_node->write();

  if(current_parent){
    if(dir == 0) current_parent->set_left_child(current_node);
    else if(dir == 1) current_parent->set_right_child(current_node);
  }
  else root = current_node;
}
/* ---------------------------------------------------------------------- */

void Tree::mutate_constant(Node *&root, double temperature)
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
    new_value = static_cast<ConstNode*>(leaves[leaf])->get_value() - const_lo + 
      temperature * new_value;
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
void Tree::mutate_variable(Node *&root)
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
void Tree::mutate_operator(Node *&root)
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
void Tree::crossover(Node *&root1, Node *&root2, int ran_depth)
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
  
  if     (op == PLUS)      rn = new PlusNode();
  else if(op == MINUS)     rn = new MinusNode();
  else if(op == TIMES)     rn = new TimesNode();
  else if(op == DIVIDE)    rn = new DivideNode();
  else if(op == POW)       rn = new PowNode();
  else if(op == CONSTANT)  rn = new ConstNode();
  else if(op == VARIABLE)  rn = new VarNode();

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

  if (op < POW+1) {
    n++; nt = from_buffer(buf, n);
    rn->set_left_child(nt);
    n++; nt = from_buffer(buf, n);
    rn->set_right_child(nt);
  }
  else if (op == CONSTANT) {static_cast<ConstNode*>(rn)->set_value(c);}
  else if (op == VARIABLE) {static_cast<VarNode*>(rn)->set_value(index);}

  return rn;
}
/* ---------------------------------------------------------------------- */
