#include <iostream>
#include <algorithm>

using namespace std;

#define RAND_MULT  214013u
#define RAND_ADD  2531011u

static unsigned SEED;

void init_random (unsigned int seed)
{ // set random seed, to generate deterministic sequence of pseudo-random numbers
  SEED = seed;
}

// generate random number between 1 and 32 (included)
unsigned myRandom()
{
  unsigned res = 1 + ((SEED>>13) & 31);
  SEED = RAND_MULT*SEED + RAND_ADD;
  return res;
}

unsigned Euclides( unsigned a, unsigned b )
{ // Compute Maximum Common Divisor
  unsigned c;
  if (a < b) // interchange values
  {
    c= b;  b= a;  a = c;
  }

  while (b)
  {
    c= b;  b= a % b;  a = c;
  }
  return a;
}


void oldChallenge ( unsigned Y, unsigned X, unsigned s )
{
  unsigned * board = new unsigned[Y*X];  // board matrix of Y*X cells
  unsigned * cost  = new unsigned[Y*X];  // cost matrix

  unsigned up, middle, down, min;
  int      y, x, arg_min;

  init_random(s);  // set initial random seed

  for ( x=0; x < X; x++ ) // Initialize first column of board with random values
    board[x*Y] = myRandom();

  for ( x=0; x < X; x++ ) // Initialize first column of cost matrix
    cost[x*Y] = board[x*Y];

  for (y=0; y<Y-1; y++) // for each column from 0 to Y-1
    for ( x=0; x<X; x++ ) // for each element x in column y
    {
      // read three elements in current column (y)
      if (x==0)   // boundary case
        up = cost[x*Y+y];
      else
        up = cost[(x-1)*Y+y];

      middle = cost[x*Y+y];

      if (x==X-1) // boundary case
        down = cost[x*Y+y];
      else
        down = cost[(x+1)*Y+y];
      // generate next data element in board (column y+1)
      if (y % 3 == 0)
        board[x*Y+y+1] = Euclides( myRandom(), 2*3*5*7*11*13);
      else if (y % 3 == 1)
      {
        unsigned v1, v2;
        v1 = myRandom();
        v2 = myRandom();
        board[x*Y+y+1] = Euclides( v1, v2 );
      }
      else
        board[x*Y+y+1] = myRandom();

      // calculate cost for this element and store it in column y+1
      min           = up < down? up+1: down+1;
      cost[x*Y+y+1] = board[x*Y+y+1] + (middle < min? middle: min);
    }

  // ******** Find minimum cost in last column *******************

  int init_pos = SEED % X; // generate random value from 0 to X-1
  arg_min = 0;

  arg_min = init_pos;
  min     = cost[init_pos*Y+Y-1];


  for ( x=init_pos+1; x<X; x++ ) // Look from init+1 to X-1
  {
    unsigned  v = cost[x*Y+Y-1];
    arg_min     = min > v? x: arg_min;
    min         = min > v? v: min;
  }

  for ( x=0; x< init_pos; x++ )  // look from 0 to init-1
  {
    unsigned  v = cost[x*Y+Y-1];
    arg_min     = min > v? x: arg_min;
    min         = min > v? v: min;
  }

  // ***** traceback the path leading to minimum cost *******
  //
  int      pos = arg_min;
  unsigned val, prev_cost = cost[pos*Y+Y-1];

  cout << "Path of cost= " << prev_cost << " ending at cell " << pos << " is ";

  for (y=Y-1; y>0; y--)
  {
    val = board[pos*Y+y];    // board value in path
    cout << val;

    prev_cost = prev_cost - val;

    if (pos==0)   // boundary case
      up = cost[  pos  *Y + y-1];
    else
      up = cost[(pos-1)*Y + y-1];

    middle = cost[  pos  *Y + y-1];
    if (pos!=X-1) // boundary case
      down = cost[(pos+1)*Y + y-1];
    else
      down = cost[  pos  *Y + y-1];

    if (prev_cost == middle)      // comes from middle cell
    {
      cout << "-";
    }
    else if (prev_cost-1 == up)   // comes from upper cell
    {
      pos--;
      cout << "-U-";
      prev_cost--;
    }
    else                          // comes from lower cell
    {
      pos++;
      cout << "-D-";
      prev_cost--;
    }
  }

  val = board[pos*Y];    // first board value in path

  cout << val << "\n";

  delete []board;
  delete []cost;
}

unsigned computeCost ( unsigned *V, unsigned Vsize, unsigned *W, unsigned Wsize )
{
  unsigned * cost        = new unsigned[Vsize*Wsize];  // cost matrix of Vsize*Wsize elements
  unsigned * addV        = new unsigned[Vsize];        // array of Vsize elements
  unsigned * addW        = new unsigned[Wsize];        // array of Vsize elements
  unsigned * addDiag     = new unsigned[Vsize+Wsize];  // array of Vsize+Wsize elements
  unsigned * addAntiDiag = new unsigned[Wsize+Vsize];  // array of Wsize+Vsize elements

  for ( int i=0; i < Vsize; i++ )
    addV[i] = 0;

  for ( int i=0; i < Wsize; i++ )
    addW[i] = 0;

  for ( int i=0; i < Vsize+Wsize; i++ ) {
    addDiag[i] = 0;
    addAntiDiag[i] = 0;
  }
  
    for ( int i=0; i < Vsize; i++ ) {
  for ( int j=0; j < Wsize; j++ ) {
      cost[i*Wsize+j] = V[i]*W[j];
    }
  }

    for ( int i=0; i < Vsize; i++ ) {
  for ( int j=0; j < Wsize; j++ ) {
      addV[i] += cost[i*Wsize+j];
      addW[j] += cost[i*Wsize+j];
    }
  }

    for ( int i=0; i < Vsize; i++ ) {
  for ( int j=0; j < Wsize; j++ ) {
      addDiag[Vsize-(i+1)+j] += cost[i*Wsize+j];
      addAntiDiag[i+j] += cost[i*Wsize+j];
    }
  }

  for ( int i=0; i < Vsize; i++ )
    addV[i] = std::max({addV[i], addDiag[i], addAntiDiag[i]});

  for ( int j=0; j < Wsize; j++ )
    addW[j] = std::max({addW[j], addDiag[Vsize+j-1], addAntiDiag[Vsize+j-1]});

   unsigned C = 0;
   for ( int i=0; i < Vsize; i++ )
    C += addV[i];

  for ( int j=0; j < Wsize; j++ )
    C += addW[j];

  delete []cost;
  delete []addV;
  delete []addW;
  delete []addDiag;
  delete []addAntiDiag;
  return C;
}


int main (int argc, char **argv)
{
  unsigned X=50000, Y= 5000, s=0, REP=1;

  // obtain arguments provided at Linux shell at run time
  if (argc>1) { REP = atoi(argv[1]); }
  if (argc>2) { X   = atoi(argv[2]); }
  if (argc>3) { Y   = atoi(argv[3]); }
  if (argc>4) { s   = atoi(argv[4]); }

  cout << "********** CHALLENGE # 1 ****************\n";
  cout << "Input vector sizes are X=" << X << " and Y=" << Y << "; Repeat " << REP << " times\n";

  unsigned * P = new unsigned[X];
  unsigned * T = new unsigned[Y];

  for (int t=0; t<REP; t++)
  {
    init_random(s);  // set initial random seed
    for (int x=0; x < X; x++ ) // Initialize P with random values
      P[x] = myRandom();
    for (int y=0; y < Y; y++ ) // Initialize T with random values
      T[y] = myRandom();

    unsigned C = computeCost ( P, X, T, Y );

    cout << "t=" << t << " Cost = " << C << "\n";
    s = s + 1;
  }

  return (0);
}
