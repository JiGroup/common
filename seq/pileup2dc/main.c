//
//  main.c
//  pileup2dc
//
//  Created by Patrick Flaherty on 3/27/12.
//  Copyright (c) 2012 Stanford University. All rights reserved.
//

#include <stdio.h>
#include <ctype.h>
#include <string.h>

#define MAXLINE 1000000 // maximum line length
enum base { A,C,G,T,N };

// data structure for a single position
typedef struct {
    int pos; // 1-base position in reference
    enum base refb; // reference base
    
    int depthF[5], depthR[5], depthT[5]; // forward/reverse depth in base enum order
} pile_t;

// get pileup file length
static int get_filelength ( FILE *fid ) {
    int numLines = 0;
    char tempString[MAXLINE];
    
    rewind(fid); // start counting from the file start
    while( !feof(fid) ) {
        fgets(tempString, MAXLINE, fid);
        numLines++;
    }
    rewind(fid);
    
    return numLines-1;
}

// compute total depth
static int sum_depth(const int *d) {
    return d[0]+d[1]+d[2]+d[3]+d[4];
}

// parse a single alignment string
static void parse_align_string(const char refb, const char *alignString, pile_t *pile) 
{
    enum base ref;
    char alignChar;
    int i;
    
    
    switch (refb) {
        case 'A':
            pile->refb=A;
            break;
        case 'C':
            pile->refb=C;
            break;
        case 'G':
            pile->refb=G;
            break;
        case 'T':
            pile->refb=T;
            break;   
        case 'N':
            pile->refb=N;
            break; 
        default:
            break;
    }
    
    int alignLen = (int)strlen(alignString);
    
    
    for (i=0; i<alignLen; i++) {
        alignChar = alignString[i];
        switch (alignChar) {
            case '.': // forward reference
                pile->depthF[ref]++;
                break;
            case ',': //reverse reference
                pile->depthR[ref]++;
                break;
                
                
            case 'A': 
                pile->depthF[0]++;
                break;
            case 'C':
                pile->depthF[1]++;
                break;
            case 'G':
                pile->depthF[2]++;
                break;
            case 'T':
                pile->depthF[3]++;
                break;
            case 'N':
                pile->depthF[4]++;
                break;
                
                // reverse read different from reference
            case 'a':
                pile->depthR[0]++;
                break;
            case 'c':
                pile->depthR[1]++;
                break;
            case 'g':
                pile->depthR[2]++;
                break;
            case 't':
                pile->depthR[3]++;
                break;
            case 'n':
                pile->depthR[4]++;
                break;
                
            case '+': //skip indels
                i=i+2;
                break;
                
            default:
                break;
        }
    }
    
    // accumulate the forward and reverse depths
    for(i = 0; i<5; i++) {
        pile->depthT[i] = pile->depthF[i]+pile->depthR[i];
    }
    
    return;
}

// get the majority vote base
static int get_majvot (const int *depth)
{
    int i;
    
    int majVot = 0, majDepth=0;
    for (i = 0; i<5; i++) {
        if (depth[i] > majDepth) {
            majVot = i;
            majDepth = depth[i];
        }
    }
    return majVot;
}


// export the data in depth chart format
static void export_depthchart(const pile_t *pile, const int numLines)
{
    char *chr="chr1";
    char baseChar[5]="ACGTN";
    int i;
    
    // Export header information
    fprintf(stdout, "count\tchr\tloc\tsubref\trefb");
    fprintf(stdout, "\t");
    
    fprintf(stdout, "A_All\tC_All\tG_All\tT_All\tN_All\tTot_All\tmatch_All\tMajVot_All");
    fprintf(stdout, "\t");
    
    fprintf(stdout, "A_F\tC_F\tG_F\tT_F\tN_F\tTot_F\tmatch_F\tMajVot_F");
    fprintf(stdout, "\t");
    
    fprintf(stdout, "A_R\tC_R\tG_R\tT_R\tN_R\tTot_R\tmatch_R\tMajVot_R");
    fprintf(stdout, "\n");
    
    for (i=0; i<numLines; i++) {
        fprintf(stdout, "%d\t%s\t%d\t%s\t%c", i+1, chr, pile[i].pos, chr, baseChar[pile[i].refb]);
        fprintf(stdout, "\t");
        fprintf(stdout, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%c", pile[i].depthT[A], pile[i].depthT[C], 
                pile[i].depthT[G], pile[i].depthT[T], pile[i].depthT[N],
                sum_depth(pile[i].depthT), pile[i].depthT[pile[i].refb], baseChar[get_majvot(pile[i].depthT)]);
        fprintf(stdout, "\t");
        fprintf(stdout, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%c", pile[i].depthF[A], pile[i].depthF[C], 
                pile[i].depthF[G], pile[i].depthF[T], pile[i].depthF[N],
                sum_depth(pile[i].depthF), pile[i].depthF[pile[i].refb], baseChar[get_majvot(pile[i].depthF)]);
        fprintf(stdout, "\t");
        fprintf(stdout, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%c", pile[i].depthR[A], pile[i].depthR[C], 
                pile[i].depthR[G], pile[i].depthR[T], pile[i].depthR[N],
                sum_depth(pile[i].depthR), pile[i].depthR[pile[i].refb], baseChar[get_majvot(pile[i].depthR)]);
        fprintf(stdout, "\n");
    }
}

int main (int argc, char *argv[])
{
    char lineBuf[MAXLINE];
    char refName[40], refb, alignString[10000], qualityString[10000];
    int pos, totDepth;
    int i;
    FILE *plfid;
    
    if (argc == 1) {
        fprintf(stderr, "Usage: pileup2dc [in.pileup]\n");
        return 1;
    }
    
    
    // open the pileup file for reading
    // die if error opening or reading file
    plfid = fopen(argv[1], "rd");
    
    
    // get the number of positions in the pileup
    int numLines = get_filelength(plfid);
    
    // allocate space for the depth chart
    pile_t pile[numLines];
    
    // read the pileup file and store in the structure
    for (i=0; i<numLines; i++) {
        fgets(lineBuf, MAXLINE, plfid);
        sscanf(lineBuf, "%s\t%d\t%c\t%d\t%s\t%s", refName, &pos, &refb, 
               &totDepth, alignString, qualityString);
        
        pile[i].pos = pos;
        parse_align_string ( toupper(refb), alignString, &pile[i] );
    }
    
    fclose(plfid);
    
    // export the pileup structure in depth chart format
    export_depthchart ( pile, numLines);
    
    return 0;
}