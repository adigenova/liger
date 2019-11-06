#ifndef ASSERTS
#define NDEBUG // disable asserts, they're computationnally intensive
#endif 

#include <stdio.h>
#include <assert.h>
#include <algorithm> // for min
#include "Kmer.h"
#include "lut.h"

using namespace std;

int sizeKmer;
uint64_t nsolids = 0;
kmer_type kmerMask;    
kmer_type kmerMaskm1;

int NT2int(char nt)
{
    int i;
    i = nt;
    i = (i>>1)&3; // that's quite clever, guillaume.
    return i;
}

int revcomp_int(int nt_int)
{   
    return (nt_int<2)?nt_int+2:nt_int-2;
}


unsigned char  code4NT(char *seq)
{
    int i;
    unsigned char x;
    
    x=0;
    for (i=0; i<4; ++i)
    {
        x = x*4 + NT2int(seq[i]);
    }
    return x;
}

unsigned char  code_n_NT(char *seq, int nb)
{
    int i;
    unsigned char x;
    
    x=0;
    for (i=0; i<nb; ++i)
    {
        x = x*4 + NT2int(seq[i]);
    }
    x = x << ((4-nb)*2) ;
    return x;
}


kmer_type  codeSeed(char *seq)
{
    return codeSeed(seq, sizeKmer, kmerMask);
}

kmer_type  codeSeed(char *seq, int sizeKmer, kmer_type kmerMask)
{
    int i;
    kmer_type x;
    
    x=0;
    for (i=0; i<sizeKmer; ++i)
    {
        x = x*4 + NT2int(seq[i]);
    }
    return x;
}
kmer_type  codeSeedRight(char *seq, kmer_type  val_seed, bool new_read)
{
    return codeSeedRight(seq, val_seed, new_read, sizeKmer, kmerMask);
}

kmer_type  codeSeedRight(char *seq, kmer_type  val_seed, bool new_read, int sizeKmer, kmer_type kmerMask)
{
    if (new_read) return codeSeed(seq, sizeKmer, kmerMask);
    // shortcut    
    return (val_seed * 4 + NT2int(seq[sizeKmer-1])) & kmerMask ;
}

kmer_type  codeSeedRight_revcomp(char *seq, kmer_type  val_seed, bool new_read)
{
    return codeSeedRight_revcomp(seq, val_seed, new_read, sizeKmer, kmerMask);    
}

kmer_type  codeSeedRight_revcomp(char *seq, kmer_type  val_seed, bool new_read, int sizeKmer, kmer_type kmerMask)
{
    
  if (new_read) return revcomp( codeSeed(seq, sizeKmer, kmerMask), sizeKmer );
  // shortcut    
  return ((val_seed >> 2) +  ( ((kmer_type) comp_NT[NT2int(seq[sizeKmer-1])]) <<  (2*(sizeKmer-1))  )  ) & kmerMask;
    
}

// warning: only call this function for sequential enumeration of kmers (no arbitrary position)
kmer_type extractKmerFromRead(char *readSeq, int position, kmer_type *graine, kmer_type *graine_revcomp)
{
    return extractKmerFromRead(readSeq, position, graine, graine_revcomp, true);
}

kmer_type extractKmerFromRead(char *readSeq, int position, kmer_type *graine, kmer_type *graine_revcomp, bool sequential)
{

    return extractKmerFromRead(readSeq, position, graine, graine_revcomp, sequential, sizeKmer, kmerMask);
}

kmer_type extractKmerFromRead(char *readSeq, int position, kmer_type *graine, kmer_type *graine_revcomp, bool sequential, int sizeKmer, kmer_type kmerMask)
{
    assert(graine != graine_revcomp); // make sure two different pointers
    bool new_read = (position == 0) || (!sequential); // faster computation for immediately overlapping kmers

        *graine = codeSeedRight(&readSeq[position], *graine, new_read, sizeKmer, kmerMask);
        *graine_revcomp = codeSeedRight_revcomp(&readSeq[position], *graine_revcomp, new_read, sizeKmer, kmerMask);
    //printf("%d sizeKmer %d %d\n",sizeKmer, graine,graine_revcomp);
    return  min(*graine,*graine_revcomp);
}


int first_nucleotide(kmer_type kmer)
{
    int result;
   result = kmer&3;
    return result;
}

int code2seq (kmer_type code, char *seq)
{
    return code2seq (code, seq, sizeKmer, kmerMask);
}

int code2seq (kmer_type code, char *seq, int sizeKmer, kmer_type kmerMask)
{
    int i;
    kmer_type temp = code;
    char bin2NT[4] = {'A','C','T','G'};

    for (i=sizeKmer-1; i>=0; i--)
    {
        seq[i]=bin2NT[first_nucleotide(temp&3)];
        temp = temp>>2;
    }
    //printf("sizeKmer = %d\n", sizeKmer);
    seq[sizeKmer]='\0';
    return sizeKmer;
}

// return the i-th nucleotide of the kmer_type kmer
int code2nucleotide( kmer_type code, int which_nucleotide)
{
    kmer_type temp = code;
    temp = temp >> (2*(sizeKmer-1-which_nucleotide));
    return first_nucleotide(temp&3);
}

uint64_t revcomp(uint64_t x, int size) {
  int i;

  uint64_t revcomp = x;
  //  printf("x %x    revcomp  %x \n",x,revcomp);
  unsigned char * kmerrev  = (unsigned char *) (&revcomp);
  unsigned char * kmer  = (unsigned char *) (&x);

  for (i=0; i<8; ++i)
    {
       kmerrev[7-i] =     revcomp_4NT[kmer[i]];
    }

    return (revcomp >> (2*( 4*sizeof(uint64_t) - size))  ) ;
}

uint64_t revcomp(uint64_t x) {
  return revcomp(x,sizeKmer);
}

// will be used by assemble()
void revcomp_sequence(char s[], int len)
{
#define CHAR_REVCOMP(a,b) {switch(a){\
	case 'A': b='T';break;case 'C': b='G';break;case 'G': b='C';break;case 'T': b='A';break;default: b=a;break;}}
		  int i;
		  unsigned char t;
		  for (i=0;i<len/2;i++)
		  {
			  t=s[i];
			  CHAR_REVCOMP(s[len-i-1],s[i]);
			  CHAR_REVCOMP(t,s[len-i-1]);
		  }
		  if (len%2==1)
			  CHAR_REVCOMP(s[len/2],s[len/2]);

}


kmer_type next_kmer(kmer_type graine, int added_nt, int *strand)
{
    assert(added_nt<4);
    assert(graine<=revcomp(graine));
    assert((strand == NULL) || (*strand<2));

    kmer_type new_graine;
    kmer_type temp_graine;

    if (strand != NULL && *strand == 1)// the kmer we're extending is actually a revcomp sequence in the bidirected debruijn graph node
        temp_graine = revcomp(graine);
    else
        temp_graine = graine;

    new_graine = (((temp_graine) * 4 )  + added_nt) & kmerMask;
    //new_graine = (((graine) >> 2 )  + ( ((kmer_type)added_nt) << ((sizeKmer-1)*2)) ) & kmerMask; // previous kmer
    kmer_type revcomp_new_graine = revcomp(new_graine);

    if (strand != NULL)
        *strand = (new_graine < revcomp_new_graine)?0:1;

    return min(new_graine,revcomp_new_graine);
}


//////////////////////////funcs for binary reads


kmer_type  codeSeed_bin(char *seq)
{
    int i;
    kmer_type x;
    
    x=0;
    for (i=0; i<sizeKmer; ++i)
    {
        x = x*4 + (seq[i]);
    }
    return x;
}

inline kmer_type  codeSeedRight_bin(char *seq, kmer_type  val_seed, bool new_read)
{
    
    if (new_read) return codeSeed_bin(seq);
    // shortcut
    return (val_seed * 4 + (seq[sizeKmer-1])) & kmerMask ;
    
}

inline kmer_type  codeSeedRight_revcomp_bin(char *seq, kmer_type  val_seed, bool new_read)
{
    
    if (new_read) return revcomp(codeSeed_bin(seq));
    // shortcut
  return ((val_seed >> 2) +  ( ((kmer_type) comp_NT[(int)(seq[sizeKmer-1])]) <<  (2*(sizeKmer-1))  )  ) & kmerMask;
    
}



kmer_type extractKmerFromRead_bin(char *readSeq, int position, kmer_type *graine, kmer_type *graine_revcomp, bool use_compressed)
{
    assert(graine != graine_revcomp); // make sure two different pointers

    bool new_read = (position == 0);
    if(!use_compressed)
    {
        *graine = codeSeedRight(&readSeq[position], *graine, new_read);
        *graine_revcomp = codeSeedRight_revcomp(&readSeq[position], *graine_revcomp, new_read);
    }
    else
    {
        *graine = codeSeedRight_bin(&readSeq[position], *graine, new_read);
        *graine_revcomp = codeSeedRight_revcomp_bin(&readSeq[position], *graine_revcomp, new_read);
    }
    return  min(*graine,*graine_revcomp);
}


// debug only: convert a kmer_type to char*
char debug_kmer_buffer[1024];
char* print_kmer(kmer_type kmer)
{
    return print_kmer(kmer,sizeKmer,kmerMask);
}

char* print_kmer(kmer_type kmer, int sizeKmer, kmer_type kmerMask)
{
    code2seq(kmer,debug_kmer_buffer, sizeKmer, kmerMask);
    return debug_kmer_buffer;
}
