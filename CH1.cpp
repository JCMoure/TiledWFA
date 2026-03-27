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

unsigned computeCost ( unsigned *V, unsigned Vsize, unsigned *W, unsigned Wsize )
{
  unsigned * cost        = new unsigned[Vsize*Wsize];  // cost matrix of Vsize*Wsize elements
  unsigned * addV        = new unsigned[Vsize];        // array of Vsize elements
  unsigned * addW        = new unsigned[Wsize];        // array of Vsize elements
  unsigned * addDiag     = new unsigned[Vsize+Wsize];  // array of Vsize+Wsize elements
  unsigned * addAntiDiag = new unsigned[Wsize+Vsize];  // array of Wsize+Vsize elements
  unsigned l, m, r, o;

  for ( int i=0; i < Vsize; i++ )
    addV[i] = 0;
    
  for ( int i=0; i < Wsize; i++ )
    addW[i] = 0;

  for ( int i=0; i < Vsize+Wsize; i++ )
    addDiag[i] = 0;
  
  for ( int i=0; i < Wsize+Vsize; i++ )
    addAntiDiag[i] = 0;
  
  for ( int i=0; i < Vsize; i++ )
    for ( int j=0; j < Wsize; j++ )
      cost[i*Wsize+j] = 0;

  for ( int j=0; j < Wsize; j++ ) {
    for ( int i=0; i < Vsize; i++ ) {
      cost[i*Wsize+j] += V[i]*W[j];
    }
  }

  for ( int j=0; j < Wsize; j++ ) {
    for ( int i=0; i < Vsize; i++ ) {
      addV[i] += cost[i*Wsize+j];
    }
  }

  for ( int j=0; j < Wsize; j++ ) {
    for ( int i=0; i < Vsize; i++ ) {
      addW[j] ^= cost[i*Wsize+j];
    }
  }

  for ( int j=0; j < Wsize; j++ ) {
    for ( int i=0; i < Vsize; i++ ) {
      addDiag[Vsize-(i+1)+j] = (addDiag[Vsize-(i+1)+j] ^ cost[i*Wsize+j])+1;
    }
  }

  for ( int j=0; j < Wsize; j++ ) {
    for ( int i=0; i < Vsize; i++ ) {
      addAntiDiag[i+j] = 2*addAntiDiag[i+j] + cost[i*Wsize+j];
    }
  }

  o = 0;
  for ( int i=0; i < Vsize-1; i++ ) {
    l = addDiag[i];
    r = addAntiDiag[i];
    m = addV[i];
    o = addV[i+1];

    addV[i+1] = o ^ ( (l+r)*(l+m) + (l+r)*(m+r) );
  }

  l = addDiag[Vsize-1];
  r = addAntiDiag[Vsize-1];
  m = addV[Vsize-1];
  o = addW[0];

  for ( int j=0; j < Wsize-1; j++ ) {
    addW[j] = o ^ ( (l+r)*(l+m) + (l+r)*(m+r) );
    l = addDiag[Vsize+j];
    r = addAntiDiag[Vsize+j];
    m = addW[j];
    o = addW[j+1];
  }
  addW[Wsize-1] = o ^ ( (l+r)*(l+m) + (l+r)*(m+r) );

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
  return C & 0xFFFF;
}


int main (int argc, char **argv)
{
  unsigned X=2000, Y= 1000, s=0, REP=100;

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
