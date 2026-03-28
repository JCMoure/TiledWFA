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

void zeroArray( unsigned *V, int Size) {
  for ( int i=0; i < Size; i++ )
    V[i] = 0;
}

void joinVectors (unsigned V1[], unsigned V2[], unsigned V[], int V1sz, int V2sz) {
  for ( int i=0; i < V1sz; i++ )
    V[i] = V1[i];
  for ( int i=0; i < V2sz; i++ )
    V[V1sz+i] = V2[i];
}

void computeMult(unsigned C[], unsigned V[], unsigned W[], int Vsz, int Wsz) {
  for ( int j=0; j < Wsz; j++ ) {
    for ( int i=0; i < Vsz; i++ ) {
      C[i*Wsz+j] += V[i]*W[j];
    }
  }
}

void addRows(unsigned C[], unsigned V[], int Vsz, int Wsz) {
  for ( int j=0; j < Wsz; j++ ) {
    for ( int i=0; i < Vsz; i++ ) {
      V[i] += C[i*Wsz+j];
    }
  }
}

void xorCols(unsigned C[], unsigned V[], int Vsz, int Wsz) {
  for ( int j=0; j < Wsz; j++ ) {
    for ( int i=0; i < Vsz; i++ ) {
      V[j] ^= C[i*Wsz+j];
    }
  }
}

void xorOneDiags(unsigned C[], unsigned V[], int Vsz, int Wsz) {
  for ( int j=0; j < Wsz; j++ ) {
    for ( int i=0; i < Vsz; i++ ) {
      V[Vsz-(i+1)+j] = (V[Vsz-(i+1)+j] ^ C[i*Wsz+j])+1;
    }
  }
}

void doubleAddAntiDiags(unsigned C[], unsigned V[], int Vsz, int Wsz) {
  for ( int j=0; j < Wsz; j++ ) {
    for ( int i=0; i < Vsz; i++ ) {
      V[i+j] = 2*V[i+j] + C[i*Wsz+j];
    }
  }
}

void combine1(unsigned Left[], unsigned Right[], unsigned V[], int Sz) {
  unsigned l, m, r, o;
  o = V[Sz-1];
  for ( int i=0; i < Sz-1; i++ ) {
    l = Left[i];
    r = Right[i];
    m = V[i];
    V[i] = o ^ ( (l+r)*(l+m) + (l+r)*(m+r) );
    o = m;
  }
}

void combine2(unsigned Left[], unsigned Right[], unsigned V[], int Sz) {
  unsigned l, m, r, o;
  o = V[Sz-1];
  for ( int i=0; i < Sz-1; i++ ) {
    l = Left[i];
    r = Right[i];
    m = V[i];
    V[i] = o ^ ( (l+r)*(l+m) + (l+r)*(m+r) );
    o = V[i];
  }
}

unsigned addVect(unsigned V[], int Sz) {
  unsigned C=0;
  for ( int i=0; i < Sz; i++ )
    C += V[i];
  return C;
}

unsigned addBits(unsigned V[], int Sz) {
  unsigned C=0;
  for ( int i=0; i < Sz; i++ )
    for ( int j=0; j<32; j++ )
      C += (V[i] & 1<<j) != 0;
  return C;
}

unsigned computeCost ( unsigned *V, unsigned Vsize, unsigned *W, unsigned Wsize )
{
  unsigned * cost        = new unsigned[Vsize*Wsize];  // cost matrix of Vsize*Wsize elements
  unsigned * addV        = new unsigned[Vsize];        // array of Vsize elements
  unsigned * addW        = new unsigned[Wsize];        // array of Vsize elements
  unsigned * addDiag     = new unsigned[Vsize+Wsize];  // array of Vsize+Wsize elements
  unsigned * addAntiDiag = new unsigned[Wsize+Vsize];  // array of Wsize+Vsize elements
  unsigned * VandW       = new unsigned[Vsize+Wsize];  // array of Wsize+Vsize elements

  zeroArray(addV,        Vsize);
  zeroArray(addW,        Wsize);
  zeroArray(addDiag,     Vsize+Wsize);
  zeroArray(addAntiDiag, Wsize+Vsize);
  zeroArray(VandW,       Vsize+Wsize);
  zeroArray(cost,        Vsize*Wsize);

  computeMult       (cost, V, W, Vsize, Wsize);
  addRows           (cost, addV, Vsize, Wsize);
  xorCols           (cost, addW, Vsize, Wsize);
  xorOneDiags       (cost, addDiag, Vsize, Wsize);
  doubleAddAntiDiags(cost, addDiag, Vsize, Wsize);

  joinVectors (addV, addW, VandW, Vsize, Wsize);
  for (int k=0; k<5*Vsize/Wsize; k++) {
    combine1 (addDiag, addAntiDiag, VandW, Vsize+Wsize);
    combine1 (VandW, addDiag, addAntiDiag, Vsize+Wsize);
    combine1 (addAntiDiag, VandW, addDiag, Vsize+Wsize);
    combine2 (addDiag, addAntiDiag, VandW, Vsize+Wsize);
    combine2 (VandW, addDiag, addAntiDiag, Vsize+Wsize);
    combine2 (addAntiDiag, VandW, addDiag, Vsize+Wsize);
  }

  unsigned C = 0;
  C += addVect( addDiag, Vsize+Wsize );
  C += addVect( addAntiDiag, Vsize+Wsize );
  C += addVect( VandW, Vsize+Wsize );
  C += addBits( addDiag, Vsize+Wsize );
  C += addBits( addAntiDiag, Vsize+Wsize );
  C += addBits( VandW, Vsize+Wsize );

  delete []cost;
  delete []addV;
  delete []addW;
  delete []addDiag;
  delete []addAntiDiag;
  delete []VandW;

  return C & 0xFFFF;
}

int main (int argc, char **argv)
{
  unsigned X=4000, Y= 500, s=0, REP=1000;

  // obtain arguments provided at Linux shell at run time
  if (argc>1) { REP = atoi(argv[1]); }
  if (argc>2) { X   = atoi(argv[2]); }
  if (argc>3) { Y   = atoi(argv[3]); }
  if (argc>4) { s   = atoi(argv[4]); }

  cout << "********** CHALLENGE # 1 ****************\n";
  cout << "Input vector sizes are X=" << X << " and Y=" << Y << "; Repeat " << REP << " times\n";

  unsigned * P = new unsigned[X];
  unsigned * T = new unsigned[Y];

  unsigned C = 0;
  for (int t=0; t<REP; t++)
  {
    init_random(s);  // set initial random seed
    for (int x=0; x < X; x++ ) // Initialize P with random values
      P[x] = myRandom();
    for (int y=0; y < Y; y++ ) // Initialize T with random values
      T[y] = myRandom();

    C += computeCost ( P, X, T, Y );
    
    if (t % 100 == 99) {
      cout << "t=" << t+1 << " Cost = " << C << "\n";    
    }
    s = s + 1;
  }

  return (0);
}
