#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "2bit.h"
#include "src/labeled_aiarray/labeled_augmented_array.h"
#include "interval_kmer.h"


int chrom_in(char *chrom, char **chrom_list, size_t n_chroms)
{   
    size_t x;
    for (x = 0; x < n_chroms; x++)
    {
         if (strcmp(chrom, chrom_list[x]) == 0)
         {
             return 1;
         }
    }
    return 0;
}


char* substr(const char *src, int m, int n)
{
    // get the length of the destination string
    int len = n - m;
 
    // allocate (len + 1) chars for destination (+1 for extra null character)
    char *dest = (char*)malloc(sizeof(char) * (len + 1));
 
    // extracts characters between m'th and n'th index from source string
    // and copy them into the destination string
    int i;
    for (i = m; i < n && (*(src + i) != '\0'); i++)
    {
        *dest = *(src + i);
        dest++;
    }
 
    // null-terminate the destination string
    *dest = '\0';
 
    // return the destination string
    return dest - len;
}


base_freq_t *base_freq_init(int n_bases)
{
    // Reserve memory
    base_freq_t *bf = malloc(1 * sizeof(base_freq_t));

    bf->A = calloc(n_bases, sizeof(float));
    bf->T = calloc(n_bases, sizeof(float));
    bf->G = calloc(n_bases, sizeof(float));
    bf->C = calloc(n_bases, sizeof(float));
    bf->n_intervals = 0;
    bf->n_bases = n_bases;

    bf->up = (int)(n_bases / 2);
    bf->down = n_bases - bf->up;

    return bf;
}


tribase_freq_t *tribase_freq_init(int n_bases)
{
    // Reserve memory
    tribase_freq_t *bf = malloc(1 * sizeof(tribase_freq_t));

    bf->AAA = calloc(n_bases, sizeof(float));
    bf->AAT = calloc(n_bases, sizeof(float));
    bf->AAG = calloc(n_bases, sizeof(float));
    bf->AAC = calloc(n_bases, sizeof(float));
    bf->ATA = calloc(n_bases, sizeof(float));
    bf->ATT = calloc(n_bases, sizeof(float));
    bf->ATG = calloc(n_bases, sizeof(float));
    bf->ATC = calloc(n_bases, sizeof(float));
    bf->AGA = calloc(n_bases, sizeof(float));
    bf->AGT = calloc(n_bases, sizeof(float));
    bf->AGG = calloc(n_bases, sizeof(float));
    bf->AGC = calloc(n_bases, sizeof(float));
    bf->ACA = calloc(n_bases, sizeof(float));
    bf->ACT = calloc(n_bases, sizeof(float));
    bf->ACG = calloc(n_bases, sizeof(float));
    bf->ACC = calloc(n_bases, sizeof(float));
    bf->TAA = calloc(n_bases, sizeof(float));
    bf->TAT = calloc(n_bases, sizeof(float));
    bf->TAG = calloc(n_bases, sizeof(float));
    bf->TAC = calloc(n_bases, sizeof(float));
    bf->TTA = calloc(n_bases, sizeof(float));
    bf->TTT = calloc(n_bases, sizeof(float));
    bf->TTG = calloc(n_bases, sizeof(float));
    bf->TTC = calloc(n_bases, sizeof(float));
    bf->TGA = calloc(n_bases, sizeof(float));
    bf->TGT = calloc(n_bases, sizeof(float));
    bf->TGG = calloc(n_bases, sizeof(float));
    bf->TGC = calloc(n_bases, sizeof(float));
    bf->TCA = calloc(n_bases, sizeof(float));
    bf->TCT = calloc(n_bases, sizeof(float));
    bf->TCG = calloc(n_bases, sizeof(float));
    bf->TCC = calloc(n_bases, sizeof(float));
    bf->GAA = calloc(n_bases, sizeof(float));
    bf->GAT = calloc(n_bases, sizeof(float));
    bf->GAG = calloc(n_bases, sizeof(float));
    bf->GAC = calloc(n_bases, sizeof(float));
    bf->GTA = calloc(n_bases, sizeof(float));
    bf->GTT = calloc(n_bases, sizeof(float));
    bf->GTG = calloc(n_bases, sizeof(float));
    bf->GTC = calloc(n_bases, sizeof(float));
    bf->GGA = calloc(n_bases, sizeof(float));
    bf->GGT = calloc(n_bases, sizeof(float));
    bf->GGG = calloc(n_bases, sizeof(float));
    bf->GGC = calloc(n_bases, sizeof(float));
    bf->GCA = calloc(n_bases, sizeof(float));
    bf->GCT = calloc(n_bases, sizeof(float));
    bf->GCG = calloc(n_bases, sizeof(float));
    bf->GCC = calloc(n_bases, sizeof(float));
    bf->CAA = calloc(n_bases, sizeof(float));
    bf->CAT = calloc(n_bases, sizeof(float));
    bf->CAG = calloc(n_bases, sizeof(float));
    bf->CAC = calloc(n_bases, sizeof(float));
    bf->CTA = calloc(n_bases, sizeof(float));
    bf->CTT = calloc(n_bases, sizeof(float));
    bf->CTG = calloc(n_bases, sizeof(float));
    bf->CTC = calloc(n_bases, sizeof(float));
    bf->CGA = calloc(n_bases, sizeof(float));
    bf->CGT = calloc(n_bases, sizeof(float));
    bf->CGG = calloc(n_bases, sizeof(float));
    bf->CGC = calloc(n_bases, sizeof(float));
    bf->CCA = calloc(n_bases, sizeof(float));
    bf->CCT = calloc(n_bases, sizeof(float));
    bf->CCG = calloc(n_bases, sizeof(float));
    bf->CCC = calloc(n_bases, sizeof(float));

    bf->n_intervals = 0;
    bf->n_bases = n_bases;
    bf->up = (int)(n_bases / 2);
    bf->down = n_bases - bf->up;

    return bf;
}


void base_freq_destroy(base_freq_t *bf)
{
    if (bf == 0) 
    {
        return;
    }

    free(bf->A);
    free(bf->T);
    free(bf->G);
    free(bf->C);

    free(bf);
}


void tribase_freq_destroy(tribase_freq_t *bf)
{
    if (bf == 0) 
    {
        return;
    }

    free(bf->AAA);
    free(bf->AAT);
    free(bf->AAG);
    free(bf->AAC);
    free(bf->ATA);
    free(bf->ATT);
    free(bf->ATG);
    free(bf->ATC);
    free(bf->AGA);
    free(bf->AGT);
    free(bf->AGG);
    free(bf->AGC);
    free(bf->ACA);
    free(bf->ACT);
    free(bf->ACG);
    free(bf->ACC);
    free(bf->TAA);
    free(bf->TAT);
    free(bf->TAG);
    free(bf->TAC);
    free(bf->TTA);
    free(bf->TTT);
    free(bf->TTG);
    free(bf->TTC);
    free(bf->TGA);
    free(bf->TGT);
    free(bf->TGG);
    free(bf->TGC);
    free(bf->TCA);
    free(bf->TCT);
    free(bf->TCG);
    free(bf->TCC);
    free(bf->GAA);
    free(bf->GAT);
    free(bf->GAG);
    free(bf->GAC);
    free(bf->GTA);
    free(bf->GTT);
    free(bf->GTG);
    free(bf->GTC);
    free(bf->GGA);
    free(bf->GGT);
    free(bf->GGG);
    free(bf->GGC);
    free(bf->GCA);
    free(bf->GCT);
    free(bf->GCG);
    free(bf->GCC);
    free(bf->CAA);
    free(bf->CAT);
    free(bf->CAG);
    free(bf->CAC);
    free(bf->CTA);
    free(bf->CTT);
    free(bf->CTG);
    free(bf->CTC);
    free(bf->CGA);
    free(bf->CGT);
    free(bf->CGG);
    free(bf->CGC);
    free(bf->CCA);
    free(bf->CCT);
    free(bf->CCG);
    free(bf->CCC);

    free(bf);

    return;
}


void base_freq_add(base_freq_t *bf, char *seq)
{
    int i;
    for (i = 0; i < bf->n_bases; i++)
    {
        if (*(seq + i) == 'A' || *(seq + i) == 'a')
        {
            bf->A[i] = bf->A[i] + 1;
        }
        else if (*(seq + i) == 'T' || *(seq + i) == 't')
        {
            bf->T[i] = bf->T[i] + 1;
        }
        else if (*(seq + i) == 'G' || *(seq + i) == 'g')
        {
            bf->G[i] = bf->G[i] + 1;
        }
        else if (*(seq + i) == 'C' || *(seq + i) == 'c')
        {
            bf->C[i] = bf->C[i] + 1;
        }
    }

    bf->n_intervals++;

    return;
}

void tribase_freq_add(tribase_freq_t *bf, char *seq)
{

    //printf("n_bases: %d\n", bf->n_bases);
    int i;
    for (i = 0; i < bf->n_bases; i++)
    {
        //printf("i: %d\n", i);
        //printf("seq: %c", *(seq + i));
        //printf("%c", *(seq + i + 1));
        //printf("%c\n", *(seq + i + 2));

        if (*(seq + i) == 'A' || *(seq + i) == 'a')
        {
            if (*(seq + i + 1) == 'A' || *(seq + i + 1) == 'a')
            {
                if (*(seq + i + 2) == 'A' || *(seq + i + 2) == 'a')
                {
                    bf->AAA[i] = bf->AAA[i] + 1;
                }
                else if (*(seq + i + 2) == 'T' || *(seq + i + 2) == 't')
                {
                    bf->AAT[i] = bf->AAT[i] + 1;
                }
                else if (*(seq + i + 2) == 'G' || *(seq + i + 2) == 'g')
                {
                    bf->AAG[i] = bf->AAG[i] + 1;
                }
                else if (*(seq + i + 2) == 'C' || *(seq + i + 2) == 'c')
                {
                    bf->AAC[i] = bf->AAC[i] + 1;
                }
            }
            else if (*(seq + i + 1) == 'T' || *(seq + i + 1) == 't')
            {
                if (*(seq + i + 2) == 'A' || *(seq + i + 2) == 'a')
                {
                    bf->ATA[i] = bf->ATA[i] + 1;
                }
                else if (*(seq + i + 2) == 'T' || *(seq + i + 2) == 't')
                {
                    bf->ATT[i] = bf->ATT[i] + 1;
                }
                else if (*(seq + i + 2) == 'G' || *(seq + i + 2) == 'g')
                {
                    bf->ATG[i] = bf->ATG[i] + 1;
                }
                else if (*(seq + i + 2) == 'C' || *(seq + i + 2) == 'c')
                {
                    bf->ATC[i] = bf->ATC[i] + 1;
                }
            }
            else if (*(seq + i + 1) == 'G' || *(seq + i + 1) == 'g')
            {
                if (*(seq + i + 2) == 'A' || *(seq + i + 2) == 'a')
                {
                    bf->AGA[i] = bf->AGA[i] + 1;
                }
                else if (*(seq + i + 2) == 'T' || *(seq + i + 2) == 't')
                {
                    bf->AGT[i] = bf->AGT[i] + 1;
                }
                else if (*(seq + i + 2) == 'G' || *(seq + i + 2) == 'g')
                {
                    bf->AGG[i] = bf->AGG[i] + 1;
                }
                else if (*(seq + i + 2) == 'C' || *(seq + i + 2) == 'c')
                {
                    bf->AGC[i] = bf->AGC[i] + 1;
                }
            }
            else if (*(seq + i + 1) == 'C' || *(seq + i + 1) == 'c')
            {
                if (*(seq + i + 2) == 'A' || *(seq + i + 2) == 'a')
                {
                    bf->ACA[i] = bf->ACA[i] + 1;
                }
                else if (*(seq + i + 2) == 'T' || *(seq + i + 2) == 't')
                {
                    bf->ACT[i] = bf->ACT[i] + 1;
                }
                else if (*(seq + i + 2) == 'G' || *(seq + i + 2) == 'g')
                {
                    bf->ACG[i] = bf->ACG[i] + 1;
                }
                else if (*(seq + i + 2) == 'C' || *(seq + i + 2) == 'c')
                {
                    bf->ACC[i] = bf->ACC[i] + 1;
                }
            }
        }
        else if (*(seq + i) == 'T' || *(seq + i) == 't')
        {
            if (*(seq + i + 1) == 'A' || *(seq + i + 1) == 'a')
            {
                if (*(seq + i + 2) == 'A' || *(seq + i + 2) == 'a')
                {
                    bf->TAA[i] = bf->TAA[i] + 1;
                }
                else if (*(seq + i + 2) == 'T' || *(seq + i + 2) == 't')
                {
                    bf->TAT[i] = bf->TAT[i] + 1;
                }
                else if (*(seq + i + 2) == 'G' || *(seq + i + 2) == 'g')
                {
                    bf->TAG[i] = bf->TAG[i] + 1;
                }
                else if (*(seq + i + 2) == 'C' || *(seq + i + 2) == 'c')
                {
                    bf->TAC[i] = bf->TAC[i] + 1;
                }
            }
            else if (*(seq + i + 1) == 'T' || *(seq + i + 1) == 't')
            {
                if (*(seq + i + 2) == 'A' || *(seq + i + 2) == 'a')
                {
                    bf->TTA[i] = bf->TTA[i] + 1;
                }
                else if (*(seq + i + 2) == 'T' || *(seq + i + 2) == 't')
                {
                    bf->TTT[i] = bf->TTT[i] + 1;
                }
                else if (*(seq + i + 2) == 'G' || *(seq + i + 2) == 'g')
                {
                    bf->TTG[i] = bf->TTG[i] + 1;
                }
                else if (*(seq + i + 2) == 'C' || *(seq + i + 2) == 'c')
                {
                    bf->TTC[i] = bf->TTC[i] + 1;
                }
            }
            else if (*(seq + i + 1) == 'G' || *(seq + i + 1) == 'g')
            {
                if (*(seq + i + 2) == 'A' || *(seq + i + 2) == 'a')
                {
                    bf->TGA[i] = bf->TGA[i] + 1;
                }
                else if (*(seq + i + 2) == 'T' || *(seq + i + 2) == 't')
                {
                    bf->TGT[i] = bf->TGT[i] + 1;
                }
                else if (*(seq + i + 2) == 'G' || *(seq + i + 2) == 'g')
                {
                    bf->TGG[i] = bf->TGG[i] + 1;
                }
                else if (*(seq + i + 2) == 'C' || *(seq + i + 2) == 'c')
                {
                    bf->TGC[i] = bf->TGC[i] + 1;
                }
            }
            else if (*(seq + i + 1) == 'C' || *(seq + i + 1) == 'c')
            {
                if (*(seq + i + 2) == 'A' || *(seq + i + 2) == 'a')
                {
                    bf->TCA[i] = bf->TCA[i] + 1;
                }
                else if (*(seq + i + 2) == 'T' || *(seq + i + 2) == 't')
                {
                    bf->TCT[i] = bf->TCT[i] + 1;
                }
                else if (*(seq + i + 2) == 'G' || *(seq + i + 2) == 'g')
                {
                    bf->TCG[i] = bf->TCG[i] + 1;
                }
                else if (*(seq + i + 2) == 'C' || *(seq + i + 2) == 'c')
                {
                    bf->TCC[i] = bf->TCC[i] + 1;
                }
            }
        }
        else if (*(seq + i) == 'C' || *(seq + i) == 'c')
        {
            if (*(seq + i + 1) == 'A' || *(seq + i + 1) == 'a')
            {
                if (*(seq + i + 2) == 'A' || *(seq + i + 2) == 'a')
                {
                    bf->CAA[i] = bf->CAA[i] + 1;
                }
                else if (*(seq + i + 2) == 'T' || *(seq + i + 2) == 't')
                {
                    bf->CAT[i] = bf->CAT[i] + 1;
                }
                else if (*(seq + i + 2) == 'G' || *(seq + i + 2) == 'g')
                {
                    bf->CAG[i] = bf->CAG[i] + 1;
                }
                else if (*(seq + i + 2) == 'C' || *(seq + i + 2) == 'c')
                {
                    bf->CAC[i] = bf->CAC[i] + 1;
                }
            }
            else if (*(seq + i + 1) == 'T' || *(seq + i + 1) == 't')
            {
                if (*(seq + i + 2) == 'A' || *(seq + i + 2) == 'a')
                {
                    bf->CTA[i] = bf->CTA[i] + 1;
                }
                else if (*(seq + i + 2) == 'T' || *(seq + i + 2) == 't')
                {
                    bf->CTT[i] = bf->CTT[i] + 1;
                }
                else if (*(seq + i + 2) == 'G' || *(seq + i + 2) == 'g')
                {
                    bf->CTG[i] = bf->CTG[i] + 1;
                }
                else if (*(seq + i + 2) == 'C' || *(seq + i + 2) == 'c')
                {
                    bf->CTC[i] = bf->CTC[i] + 1;
                }
            }
            else if (*(seq + i + 1) == 'G' || *(seq + i + 1) == 'g')
            {
                if (*(seq + i + 2) == 'A' || *(seq + i + 2) == 'a')
                {
                    bf->CGA[i] = bf->CGA[i] + 1;
                }
                else if (*(seq + i + 2) == 'T' || *(seq + i + 2) == 't')
                {
                    bf->CGT[i] = bf->CGT[i] + 1;
                }
                else if (*(seq + i + 2) == 'G' || *(seq + i + 2) == 'g')
                {
                    bf->CGG[i] = bf->CGG[i] + 1;
                }
                else if (*(seq + i + 2) == 'C' || *(seq + i + 2) == 'c')
                {
                    bf->CGC[i] = bf->CGC[i] + 1;
                }
            }
            else if (*(seq + i + 1) == 'C' || *(seq + i + 1) == 'c')
            {
                if (*(seq + i + 2) == 'A' || *(seq + i + 2) == 'a')
                {
                    bf->CCA[i] = bf->CCA[i] + 1;
                }
                else if (*(seq + i + 2) == 'T' || *(seq + i + 2) == 't')
                {
                    bf->CCT[i] = bf->CCT[i] + 1;
                }
                else if (*(seq + i + 2) == 'G' || *(seq + i + 2) == 'g')
                {
                    bf->CCG[i] = bf->CCG[i] + 1;
                }
                else if (*(seq + i + 2) == 'C' || *(seq + i + 2) == 'c')
                {
                    bf->CCC[i] = bf->CCC[i] + 1;
                }
            }
        }
        else if (*(seq + i) == 'G' || *(seq + i) == 'g')
        {
            if (*(seq + i + 1) == 'A' || *(seq + i + 1) == 'a')
            {
                if (*(seq + i + 2) == 'A' || *(seq + i + 2) == 'a')
                {
                    bf->GAA[i] = bf->GAA[i] + 1;
                }
                else if (*(seq + i + 2) == 'T' || *(seq + i + 2) == 't')
                {
                    bf->GAT[i] = bf->GAT[i] + 1;
                }
                else if (*(seq + i + 2) == 'G' || *(seq + i + 2) == 'g')
                {
                    bf->GAG[i] = bf->GAG[i] + 1;
                }
                else if (*(seq + i + 2) == 'C' || *(seq + i + 2) == 'c')
                {
                    bf->GAC[i] = bf->GAC[i] + 1;
                }
            }
            else if (*(seq + i + 1) == 'T' || *(seq + i + 1) == 't')
            {
                if (*(seq + i + 2) == 'A' || *(seq + i + 2) == 'a')
                {
                    bf->GTA[i] = bf->GTA[i] + 1;
                }
                else if (*(seq + i + 2) == 'T' || *(seq + i + 2) == 't')
                {
                    bf->GTT[i] = bf->GTT[i] + 1;
                }
                else if (*(seq + i + 2) == 'G' || *(seq + i + 2) == 'g')
                {
                    bf->GTG[i] = bf->GTG[i] + 1;
                }
                else if (*(seq + i + 2) == 'C' || *(seq + i + 2) == 'c')
                {
                    bf->GTC[i] = bf->GTC[i] + 1;
                }
            }
            else if (*(seq + i + 1) == 'G' || *(seq + i + 1) == 'g')
            {
                if (*(seq + i + 2) == 'A' || *(seq + i + 2) == 'a')
                {
                    bf->GGA[i] = bf->GGA[i] + 1;
                }
                else if (*(seq + i + 2) == 'T' || *(seq + i + 2) == 't')
                {
                    bf->GGT[i] = bf->GGT[i] + 1;
                }
                else if (*(seq + i + 2) == 'G' || *(seq + i + 2) == 'g')
                {
                    bf->GGG[i] = bf->GGG[i] + 1;
                }
                else if (*(seq + i + 2) == 'C' || *(seq + i + 2) == 'c')
                {
                    bf->GGC[i] = bf->GGC[i] + 1;
                }
            }
            else if (*(seq + i + 1) == 'C' || *(seq + i + 1) == 'c')
            {
                if (*(seq + i + 2) == 'A' || *(seq + i + 2) == 'a')
                {
                    bf->GCA[i] = bf->GCA[i] + 1;
                }
                else if (*(seq + i + 2) == 'T' || *(seq + i + 2) == 't')
                {
                    bf->GCT[i] = bf->GCT[i] + 1;
                }
                else if (*(seq + i + 2) == 'G' || *(seq + i + 2) == 'g')
                {
                    bf->GCG[i] = bf->GCG[i] + 1;
                }
                else if (*(seq + i + 2) == 'C' || *(seq + i + 2) == 'c')
                {
                    bf->GCC[i] = bf->GCC[i] + 1;
                }
            }
        }
    }

    bf->n_intervals++;

    return;
}


void base_freq_normalize(base_freq_t *bf)
{
    float sum = 0;
    int i;
    for (i = 0; i < bf->n_bases; i++)
    {
        sum = bf->A[i] + bf->T[i] + bf->G[i] + bf->C[i];
        if (sum == 0)
        {
            continue;
        }

        bf->A[i] = bf->A[i] / sum;
        bf->T[i] = bf->T[i] / sum;
        bf->G[i] = bf->G[i] / sum;
        bf->C[i] = bf->C[i] / sum;
    }

    return;
}


void tribase_freq_normalize(tribase_freq_t *bf)
{
    float sum = 0;
    int i;
    for (i = 0; i < bf->n_bases; i++)
    {
        sum = bf->AAA[i] + bf->AAT[i] + 
                bf->AAG[i] + bf->AAC[i] + 
                bf->ATA[i] + bf->ATT[i] + 
                bf->ATG[i] + bf->ATC[i] + 
                bf->AGA[i] + bf->AGT[i] + 
                bf->AGG[i] + bf->AGC[i] + 
                bf->ACA[i] + bf->ACT[i] + 
                bf->ACG[i] + bf->ACC[i] + 
                bf->TAA[i] + bf->TAT[i] + 
                bf->TAG[i] + bf->TAC[i] + 
                bf->TTA[i] + bf->TTT[i] + 
                bf->TTG[i] + bf->TTC[i] + 
                bf->TGA[i] + bf->TGT[i] + 
                bf->TGG[i] + bf->TGC[i] + 
                bf->TCA[i] + bf->TCT[i] + 
                bf->TCG[i] + bf->TCC[i] + 
                bf->GAA[i] + bf->GAT[i] + 
                bf->GAG[i] + bf->GAC[i] + 
                bf->GTA[i] + bf->GTT[i] + 
                bf->GTG[i] + bf->GTC[i] + 
                bf->GGA[i] + bf->GGT[i] + 
                bf->GGG[i] + bf->GGC[i] + 
                bf->GCA[i] + bf->GCT[i] + 
                bf->GCG[i] + bf->GCC[i] +
                bf->CAA[i] + bf->CAT[i] +
                bf->CAG[i] + bf->CAC[i] +
                bf->CTA[i] + bf->CTT[i] +
                bf->CTG[i] + bf->CTC[i] +
                bf->CGA[i] + bf->CGT[i] +
                bf->CGG[i] + bf->CGC[i] +
                bf->CCA[i] + bf->CCT[i] +
                bf->CCG[i] + bf->CCC[i];
        
        if (sum == 0)
        {
            continue;
        }

        bf->AAA[i] = bf->AAA[i] / sum;
        bf->AAT[i] = bf->AAT[i] / sum;
        bf->AAG[i] = bf->AAG[i] / sum;
        bf->AAC[i] = bf->AAC[i] / sum;
        bf->ATA[i] = bf->ATA[i] / sum;
        bf->ATT[i] = bf->ATT[i] / sum;
        bf->ATG[i] = bf->ATG[i] / sum;
        bf->ATC[i] = bf->ATC[i] / sum;
        bf->AGA[i] = bf->AGA[i] / sum;
        bf->AGT[i] = bf->AGT[i] / sum;
        bf->AGG[i] = bf->AGG[i] / sum;
        bf->AGC[i] = bf->AGC[i] / sum;
        bf->ACA[i] = bf->ACA[i] / sum;
        bf->ACT[i] = bf->ACT[i] / sum;
        bf->ACG[i] = bf->ACG[i] / sum;
        bf->ACC[i] = bf->ACC[i] / sum;
        bf->TAA[i] = bf->TAA[i] / sum;
        bf->TAT[i] = bf->TAT[i] / sum;
        bf->TAG[i] = bf->TAG[i] / sum;
        bf->TAC[i] = bf->TAC[i] / sum;
        bf->TTA[i] = bf->TTA[i] / sum;
        bf->TTT[i] = bf->TTT[i] / sum;
        bf->TTG[i] = bf->TTG[i] / sum;
        bf->TTC[i] = bf->TTC[i] / sum;
        bf->TGA[i] = bf->TGA[i] / sum;
        bf->TGT[i] = bf->TGT[i] / sum;
        bf->TGG[i] = bf->TGG[i] / sum;
        bf->TGC[i] = bf->TGC[i] / sum;
        bf->TCA[i] = bf->TCA[i] / sum;
        bf->TCT[i] = bf->TCT[i] / sum;
        bf->TCG[i] = bf->TCG[i] / sum;
        bf->TCC[i] = bf->TCC[i] / sum;
        bf->GAA[i] = bf->GAA[i] / sum;
        bf->GAT[i] = bf->GAT[i] / sum;
        bf->GAG[i] = bf->GAG[i] / sum;
        bf->GAC[i] = bf->GAC[i] / sum;
        bf->GTA[i] = bf->GTA[i] / sum;
        bf->GTT[i] = bf->GTT[i] / sum;
        bf->GTG[i] = bf->GTG[i] / sum;
        bf->GTC[i] = bf->GTC[i] / sum;
        bf->GGA[i] = bf->GGA[i] / sum;
        bf->GGT[i] = bf->GGT[i] / sum;
        bf->GGG[i] = bf->GGG[i] / sum;
        bf->GGC[i] = bf->GGC[i] / sum;
        bf->GCA[i] = bf->GCA[i] / sum;
        bf->GCT[i] = bf->GCT[i] / sum;
        bf->GCG[i] = bf->GCG[i] / sum;
        bf->GCC[i] = bf->GCC[i] / sum;
        bf->CAA[i] = bf->CAA[i] / sum;
        bf->CAT[i] = bf->CAT[i] / sum;
        bf->CAG[i] = bf->CAG[i] / sum;
        bf->CAC[i] = bf->CAC[i] / sum;
        bf->CTA[i] = bf->CTA[i] / sum;
        bf->CTT[i] = bf->CTT[i] / sum;
        bf->CTG[i] = bf->CTG[i] / sum;
        bf->CTC[i] = bf->CTC[i] / sum;
        bf->CGA[i] = bf->CGA[i] / sum;
        bf->CGT[i] = bf->CGT[i] / sum;
        bf->CGG[i] = bf->CGG[i] / sum;
        bf->CGC[i] = bf->CGC[i] / sum;
        bf->CCA[i] = bf->CCA[i] / sum;
        bf->CCT[i] = bf->CCT[i] / sum;
        bf->CCG[i] = bf->CCG[i] / sum;
        bf->CCC[i] = bf->CCC[i] / sum;
    }

    return;
}


interval_base_freq_t *interval_base_freq_init(int n_intervals, int n_bases)
{
    // Reserve memory
    interval_base_freq_t *ibf = malloc(1 * sizeof(interval_base_freq_t));

    ibf->start = base_freq_init(n_bases);
    ibf->end = base_freq_init(n_bases);

    return ibf;
}


interval_tribase_freq_t *interval_tribase_freq_init(int n_bases)
{
    // Reserve memory
    interval_tribase_freq_t *ibf = malloc(1 * sizeof(interval_tribase_freq_t));

    ibf->start = tribase_freq_init(n_bases);
    ibf->end = tribase_freq_init(n_bases);

    return ibf;
}


void interval_base_freq_destroy(interval_base_freq_t *ibf)
{
    if (ibf == 0) 
    {
        return;
    }

    base_freq_destroy(ibf->start);
    base_freq_destroy(ibf->end);

    free(ibf);
}


void interval_tribase_freq_destroy(interval_tribase_freq_t *ibf)
{
    if (ibf == 0) 
    {
        return;
    }

    tribase_freq_destroy(ibf->start);
    tribase_freq_destroy(ibf->end);

    free(ibf);
}


void interval_base_freq_add(interval_base_freq_t *ibf, TwoBit *tb, char *name, uint32_t start, uint32_t end)
{
    char *start_seq = twobitSequence(tb, name, start - ibf->start->up, start + ibf->start->down);
    char *end_seq = twobitSequence(tb, name, end - ibf->end->up, end + ibf->end->down);

    base_freq_add(ibf->start, start_seq);
    base_freq_add(ibf->end, end_seq);

    free(start_seq);
    free(end_seq);

    return;
}


void interval_tribase_freq_add(interval_tribase_freq_t *ibf, TwoBit *tb, char *name, uint32_t start, uint32_t end)
{
    char *start_seq = twobitSequence(tb, name, start - ibf->start->up, start + ibf->start->down + 3);
    char *end_seq = twobitSequence(tb, name, end - ibf->end->up, end + ibf->end->down + 3);

    tribase_freq_add(ibf->start, start_seq);
    tribase_freq_add(ibf->end, end_seq);

    free(start_seq);
    free(end_seq);

    return;
}


void interval_base_freq_normalize(interval_base_freq_t *ibf)
{
    base_freq_normalize(ibf->start);
    base_freq_normalize(ibf->end);

    return;
}


void interval_tribase_freq_normalize(interval_tribase_freq_t *ibf)
{
    tribase_freq_normalize(ibf->start);
    tribase_freq_normalize(ibf->end);

    return;
}


interval_base_freq_t *read_interval_base_freq(labeled_aiarray_t *laia, char *fname, int n_bases)
{
    // Open 2bit file
    TwoBit *tb = twobitOpen(fname, 0);
    interval_base_freq_t *ibf = interval_base_freq_init(laia->n_labels, n_bases);
    
    int t;
    for (t = 0; t < laia->n_labels; t++)
    {
        label_t *p = &laia->labels[t];
        char *name = p->name;
        uint32_t chrom_len = twobitChromLen(tb, name);

        // Check if chromosome is in the 2bit file
        if (chrom_in(name, tb->cl->chrom, tb->hdr->nChroms) == 0)
        {
            continue;
        }

        int i;
        for (i = 0; i < p->ail->nr; i++)
        {
            if (p->ail->interval_list[i].start > n_bases && p->ail->interval_list[i].end < (chrom_len - n_bases))
            {
                interval_base_freq_add(ibf, tb, name, p->ail->interval_list[i].start, p->ail->interval_list[i].end);
            }
        }

    }

    // Close 2bit file
    twobitClose(tb);

    interval_base_freq_normalize(ibf);

    return ibf;
}


interval_tribase_freq_t *read_interval_tribase_freq(labeled_aiarray_t *laia, char *fname, int n_bases)
{
    // Open 2bit file
    TwoBit *tb = twobitOpen(fname, 0);
    //printf("Initializing interval_tribase_freq\n");
    interval_tribase_freq_t *ibf = interval_tribase_freq_init(n_bases);
    //printf("   Done\n");

    //printf("Iterating over intervals\n");
    
    int t;
    for (t = 0; t < laia->n_labels; t++)
    {
        label_t *p = &laia->labels[t];
        char *name = p->name;
        uint32_t chrom_len = twobitChromLen(tb, name);

        // Check if chromosome is in the 2bit file
        if (chrom_in(name, tb->cl->chrom, tb->hdr->nChroms) == 0)
        {
            continue;
        }

        int i;
        for (i = 0; i < p->ail->nr; i++)
        {
            if (p->ail->interval_list[i].start > n_bases && p->ail->interval_list[i].end < (chrom_len - n_bases))
            {
                //printf("   Adding interval %d\n", i);
                interval_tribase_freq_add(ibf, tb, name, p->ail->interval_list[i].start, p->ail->interval_list[i].end);
                //printf("   Done\n");
            }
        }
    }

    //printf("   Done with intervals\n");

    // Close 2bit file
    twobitClose(tb);

    //printf("Normalizing\n");
    interval_tribase_freq_normalize(ibf);
    //printf("   Done\n");

    return ibf;
}


kmer_count_t *kmer_count_init(int kmer)
{
    // Reserve memory
    kmer_count_t *kc = malloc(1 * sizeof(kmer_count_t));

    kc->max_kmers = 64;
    kc->kmers = malloc(kc->max_kmers * sizeof(kmer_t));
    kc->n_kmers = 0;
    kc->kmer_lookup = kh_init(khStrInt);

    return kc;

}

void kmer_count_destroy(kmer_count_t *kc)
{
    int32_t i;
	if (kc == 0) 
    {
        return;
    }

	for (i = 0; i < kc->n_kmers; ++i)
    {
		free((char*)kc->kmers[i].name);
	}
	free(kc->kmers);
	kh_destroy(khStrInt, (strhash_t*)kc->kmer_lookup);

	free(kc);
}


void add_kmer(kmer_count_t *kc, char *kmer_name)
{
    khiter_t k;
	strhash_t *h = (strhash_t*)kc->kmer_lookup;

	k = kh_get(khStrInt, h, kmer_name);
	if (k == kh_end(h))
    {
        if (kc->n_kmers == kc->max_kmers)
        {
			EXPAND(kc->kmers, kc->max_kmers);
        }

        // Add label_name to label_map
        int ret_name = kh_name_set(khStrInt, h, kmer_name, kc->n_kmers);
        kc->n_kmers++;

        // Determine label code
		uint32_t t = kh_value(h, k);

		kmer_t *p = &kc->kmers[t];
		p->name = strdup(kmer_name);
		p->count = 0;
    }

    // Determine label code
	uint32_t t = kh_value(h, k);
    kmer_t *p = &kc->kmers[t];
    p->count++;

    return;
}


int32_t get_kmer(kmer_count_t *kc, char *kmer)
{   /* Return index for given label */
	
    khint_t k;
	strhash_t *h = (strhash_t*)kc->kmer_lookup;
	k = kh_get(khStrInt, h, kmer);
    
	return k == kh_end(h)? -1 : kh_val(h, k);
}


void append_kmers(kmer_count_t *kc, int kmer, char *seq)
{
    int seq_length = strlen(seq);
    int i;
    for (i = 0; i <= seq_length - kmer; i++)
    {
        char *kmer_seq = substr(seq, i, i+kmer);
        //printf("i: %d seq_length: %d kmer: %s\n", i, seq_length, kmer_seq);
        add_kmer(kc, kmer_seq);
        free(kmer_seq);
    }

    return;
}


int fetch_kmer(kmer_count_t *kc, char *seq)
{
    uint32_t t = get_kmer(kc, seq);

    int count = 0;
    if (t >= 0)
    {
        count = kc->kmers[t].count;
    }
    
    return count;
}


kmer_count_t *interval_kmer_count(labeled_aiarray_t *laia, char *fname, int kmer, int last_n)
{
    kmer_count_t *kc = kmer_count_init(kmer);
    TwoBit *tb = twobitOpen(fname, 0);
    
    int t;
    for (t = 0; t < laia->n_labels; t++)
    {
        label_t *p = &laia->labels[t];
        char *name = p->name;

        // Check if chromosome is in the 2bit file
        if (chrom_in(name, tb->cl->chrom, tb->hdr->nChroms) == 0)
        {
            continue;
        }

        int i;
        if (last_n == 0)
        {
            for (i = 0; i < p->ail->nr; i++)
            {
                char *seq = twobitSequence(tb, name, p->ail->interval_list[i].start, p->ail->interval_list[i].end);
                append_kmers(kc, kmer, seq);
                free(seq);
            }
        }
        else {
            for (i = 0; i < p->ail->nr; i++)
            {
                if ((p->ail->interval_list[i].end - p->ail->interval_list[i].start) < last_n)
                {
                    continue;
                }

                char *seq = twobitSequence(tb, name, p->ail->interval_list[i].end - last_n, p->ail->interval_list[i].end);
                append_kmers(kc, kmer, seq);
                free(seq);
            }
        }
    }

    // Close 2bit file
    twobitClose(tb);

    return kc;
}


char *fetch_sequence(char *fname, char *name, int start, int end)
{
    // Open 2bit file
    TwoBit *tb = twobitOpen(fname, 0);

    // Fetch sequence
    char *seq = twobitSequence(tb, name, start, end);

    // Close 2bit file
    twobitClose(tb);

    return seq;
}


int base2code(char base)
{
    int code = -1;
    if (base == 'A' || base == 'a')
    {
        code = 0;
    }
    else if (base == 'T' || base == 't')
    {
        code = 1;
    }
    else if (base == 'G' || base == 'g')
    {
        code = 2;
    }
    else if (base == 'C' || base == 'c')
    {
        code = 3;
    }

    return code;
}


void fetch_sequence_code(char *fname, char *name, int start, int end, int *seq_code)
{
    // Open 2bit file
    TwoBit *tb = twobitOpen(fname, 0);

    // Fetch sequence
    char *seq = twobitSequence(tb, name, start, end);

    // Convert sequence to code
    int seq_length = strlen(seq);
    int i;
    for (i = 0; i < seq_length; i++)
    {
        switch (*(seq + i))
        {
            case 'A':
            case 'a':
                seq_code[i] = 0;
                break;
            case 'T':
            case 't':
                seq_code[i] = 1;
                break;
            case 'G':
            case 'g':
                seq_code[i] = 2;
                break;
            case 'C':
            case 'c':
                seq_code[i] = 3;
                break;
            default:
                seq_code[i] = -1;
        }
    }

    // Close 2bit file
    twobitClose(tb);

    return;
}


void gc_content(labeled_aiarray_t *laia, char *fname, float gc[])
{
    TwoBit *tb = twobitOpen(fname, 0);

    labeled_aiarray_iter_t *iter = labeled_aiarray_iter_init(laia);
    while (labeled_aiarray_iter_next(iter) == 1)
    {
        // Check if chromosome is in the 2bit file
        char const *name = iter->intv->name;
        if (chrom_in(name, tb->cl->chrom, tb->hdr->nChroms) == 0)
        {
            continue;
        }

        // Check if interval is within chromosome length
        if (iter->intv->i->end > twobitChromLen(tb, name))
        {
            continue;
        }

        // Fetch sequence
        char *seq = twobitSequence(tb, name, iter->intv->i->start, iter->intv->i->end);
        int seq_length = strlen(seq);

        // Calculate GC content
        int i;
        for (i = 0; i <= seq_length - 1; i++)
        {
            if (*(seq + i) == 'G' || *(seq + i) == 'C' || *(seq + i) == 'g' || *(seq + i) == 'c')
            {
                gc[iter->n] = gc[iter->n] + 1;
            }
        }
        
        // Free sequence
        free(seq);

        // Normalize GC content
        gc[iter->n] = gc[iter->n] / (float)seq_length;
    }

    // Close 2bit file
    twobitClose(tb);
    labeled_aiarray_iter_destroy(iter);

    return;
}