/************************************************************************************\
*                                                                                    *
* Copyright (c) 2014, Dr. Eugene W. Myers (EWM). All rights reserved.                *
*                                                                                    *
* Redistribution and use in source and binary forms, with or without modification,   *
* are permitted provided that the following conditions are met:                      *
*                                                                                    *
*  · Redistributions of source code must retain the above copyright notice, this     *
*    list of conditions and the following disclaimer.                                *
*                                                                                    *
*  · Redistributions in binary form must reproduce the above copyright notice, this  *
*    list of conditions and the following disclaimer in the documentation and/or     *
*    other materials provided with the distribution.                                 *
*                                                                                    *
*  · The name of EWM may not be used to endorse or promote products derived from     *
*    this software without specific prior written permission.                        *
*                                                                                    *
* THIS SOFTWARE IS PROVIDED BY EWM ”AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES,    *
* INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND       *
* FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL EWM BE LIABLE   *
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES *
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS  *
* OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY      *
* THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING     *
* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN  *
* IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                                      *
*                                                                                    *
* For any issues regarding this software and its use, contact EWM at:                *
*                                                                                    *
*   Eugene W. Myers Jr.                                                              *
*   Bautzner Str. 122e                                                               *
*   01099 Dresden                                                                    *
*   GERMANY                                                                          *
*   Email: gene.myers@gmail.com                                                      *
*                                                                                    *
\************************************************************************************/

/*******************************************************************************************
 *
 *  Utility for writing the overlaps in a .las file into Graphical Fragment Assembly format
 *  This program is derived from LAshow by EWM
 *
 *  Author:  Jared Simpson
 *  Date  :  September 2014
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <assert.h>

#include "DB.h"
#include "align.h"

static char *Usage[] =
    { "[-coU] [-(a|r):<db>] [-i<int(4)>] [-w<int(100)>] [-b<int(10)>] ",
      "       <align:las> [ <reads:range> ... ]"
    };

void Print_Cigar(FILE *file, Alignment *align, int a_is_base_sequence)
{
  int i;

  // Get the start and endpoints of the local alignment on read A and B
  int abeg = align->path->abpos;
  int bbeg = align->path->bbpos;

  int aend = align->path->aepos;
  int bend = align->path->bepos;
  
  fprintf(stderr, "CIGAR -- LA coordinates A: [%d %d] B: [%d %d]\n", abeg, aend, bbeg, bend);
  
  // Grab the path trace from the alignment
  int *trace = align->path->trace;
  int tlen  = align->path->tlen;
  
  int acurr = abeg;
  int bcurr = bbeg;
  
  i = 0;
  while(i < tlen)
  { int md = 0;

    // Calculate the number of matched characters before this gap
    if(trace[i] < 0)
      md = -trace[i] - acurr - 1;
    else
      md =  trace[i] - bcurr - 1;

    // Handle the case where we have to introduce multiple gaps at the
    // same position
    int j = i;
    while(j < tlen && trace[i] == trace[++j]) {}

    // Update current positions
    acurr += md;
    bcurr += md;

    // Skip the gap bases in the appropriate read
    char type;
    if(trace[i] < 0)
    { bcurr += j - i;
      type = a_is_base_sequence ? 'I' : 'D';
    } 
    else
    { acurr += j - i;
      type = a_is_base_sequence ? 'D' : 'I';
    }

    fprintf(file, "%dM%d%c", md, j - i, type);

#ifdef DEBUG_CIGAR
    fprintf(stderr, "  %3d %dM G: %d a: %d b: %d\n", trace[i], md, j-i, acurr, bcurr);
    int k;
    for(k = 0; k < md; ++k)
      fprintf(stderr, "%c", "ACGT"[align->aseq[acurr + k]]);
    fprintf(stderr, "\n");
    
    for(k = 0; k < md; ++k)
      fprintf(stderr, "%c", "ACGT"[align->bseq[bcurr + k]]);
    fprintf(stderr, "\n");
#endif

    i = j;
  }

  // Print the trailing matches
  int md = aend - acurr;
  assert(bend - bcurr == md);
  fprintf(file, "%dM", md);
}

int main(int argc, char *argv[])
{ HITS_DB   _db,  *db  = &_db; 
  Overlap   _ovl, *ovl = &_ovl;
  Alignment _aln, *aln = &_aln;

  FILE   *input;
  int64   novl;
  int     tspace, tbytes, small;

  int     ALIGN = 0;

  //  Process options

  { int    i, j;
    int    flags[128];

    ARG_INIT("LA2gfa")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            break;
          case 'a':
            ALIGN = 1;
            if (argv[i][2] != ':')
              { fprintf(stderr,"%s: Unrecognizable option %s\n",Prog_Name,argv[i]);
                exit (1);
              }
            if (Open_DB(argv[i]+3,db))
              exit (1);
            Trim_DB(db);
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    if (argc <= 1)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        exit (1);
      }
    (void)flags; //shhh warnings
  }
  

  //  Initiate file reading and read (novl, tspace) header
  
  { char  *over, *pwd, *root;

    pwd   = PathTo(argv[1]);
    root  = Root(argv[1],".las");
    over  = Catenate(pwd,"/",root,".las");
    input = Fopen(over,"r");
    if (input == NULL)
      exit (1);

    if (fread(&novl,sizeof(int64),1,input) != 1)
      SYSTEM_ERROR
    if (fread(&tspace,sizeof(int),1,input) != 1)
      SYSTEM_ERROR

    if (tspace <= TRACE_XOVR)
      { small  = 1;
        tbytes = sizeof(uint8);
      }
    else
      { small  = 0;
        tbytes = sizeof(uint16);
      }

    free(pwd);
    free(root);
    free(over);
  }

  // Write the header
  char* version_string = "VN:Z:1.0";
  fprintf(stdout, "H\t%s\n", version_string);

  // Read the DB and write out the segment records
  { int j;
    char* seq_buffer = New_Read_Buffer(db);
    char name_buffer[MAX_NAME];

    for (j = 0; j < db->oreads; j++)
      { Load_Read(db,j,seq_buffer,2);

        snprintf(name_buffer, MAX_NAME, "%d", j);
        fprintf(stdout, "S\t%s\t%s\n", name_buffer, seq_buffer);
      }
    free(seq_buffer - 1); // See New_Read_Buffer in DB.h
  }

  //  Read the overlap file and display selected records
  
  { int        j;
    uint16    *trace;
    Work_Data *work;
    int        tmax;

    work = New_Work_Data();
    aln->path = &(ovl->path);
    aln->aseq = New_Read_Buffer(db);
    aln->bseq = New_Read_Buffer(db);

    tmax  = 1000;
    trace = (uint16 *) Malloc(sizeof(uint16)*tmax,"Allocating trace vector");
    if (trace == NULL)
      exit (1);

    //  For each overlap record do
    for (j = 0; j < novl; j++)

       //  Read it in
      { Read_Overlap(input,ovl);
        if (ovl->path.tlen > tmax)
          { tmax = ((int) 1.2*ovl->path.tlen) + 100;
            trace = (uint16 *) Realloc(trace,sizeof(uint16)*tmax,"Allocating trace vector");
            if (trace == NULL)
              exit (1);
          }

        ovl->path.trace = (void *) trace;
        Read_Trace(input,ovl,tbytes);
        
        if (small)
          Decompress_TraceTo16(ovl);

        aln->alen  = ovl->alen;
        aln->blen  = ovl->blen;
        aln->flags = ovl->flags;
        Load_Read(db,ovl->aread,aln->aseq,0);
        Load_Read(db,ovl->bread,aln->bseq,0);
        
        if (COMP(aln->flags))
          Complement_Seq(aln->bseq);

        int abeg = aln->path->abpos;
        int bbeg = aln->path->bbpos;

        int aend = aln->path->aepos;
        int bend = aln->path->bepos;
        
        int alen = aln->alen;
        int blen = aln->blen;
        
        int a_left_match = abeg == 0;
        int a_right_match = aend == alen;
        
        int b_left_match = bbeg == 0;
        int b_right_match = bend == blen;

        char a_name[MAX_NAME];
        char b_name[MAX_NAME];

        // Lazily use read index as ID right now
        snprintf(a_name, MAX_NAME, "%d", ovl->aread);
        snprintf(b_name, MAX_NAME, "%d", ovl->bread);

        fprintf(stderr, "%s %s %d %d %d %d\n", a_name, b_name, abeg, aend, bbeg, bend);
        fprintf(stderr, "%s %s %d %d %d %d\n", a_name, b_name, a_left_match, a_right_match, b_left_match, b_right_match);
        
        // Check for local local alignments and skip
        if(!a_left_match && !a_right_match && !b_left_match && !b_right_match)
          { fprintf(stderr, "Skipping local alignment\n");
            continue;
          }
        
        // When printing containments, we might switch
        // which read is the the first sequence printed.
        // The Print_CIGAR function needs to know where read
        // deletions refer to and which read insertions refer to.
        int a_is_base_sequence = 1;

        // Handle containnment records
        if( (a_left_match && a_right_match) ||
            (b_left_match && b_right_match))
          { 
            fprintf(stderr, "Containment\n");
  
            char* container_name;
            char* contained_name;

            // Determine the record that is contained within the other record
            int position = 0;
            if(a_left_match)
              { contained_name = a_name;
                container_name = b_name;
                position = aln->path->bbpos;
                a_is_base_sequence = 0;
              }
            else
              { contained_name = b_name;
                contained_name = a_name;
                position = aln->path->abpos;
                a_is_base_sequence = 1;
              }

             char ori = COMP(aln->flags) ? '-' : '+';
             fprintf(stdout, "C\t%s\t%c\t%s\t%c\t%d\t", container_name, '+',
                                                        contained_name, ori,
                                                        position);
          }
        else
          { //
            // GFA link orientations:
            //
            // S1 -------------->
            // S2     --------------->
            // Record: S1 + S2 +

            // S1    -------------->
            // S2 ------------->
            // Record: S2 + S1 +

            // S1 -------------->
            // S2     <---------------
            // Record: S1 + S2 -

            // S1       <--------------
            // S2  -------------->
            // Record: S2 + S2 -
            char* first_name;
            char* second_name;

            char first_ori = '+';
            char second_ori = '+';

            // One of the two reads must have a suffix overlap
            // and one of the two must have a prefix overlap (after
            // accounting for reverse complementation)
            // We write the ID of the read with the suffix overlap
            // first.
            // The orientation flags are always set such that the "a"
            // sequence is on the + strand
            if(a_right_match)
              { assert(b_left_match);
                first_name = a_name;
                second_name = b_name;
                first_ori = '+';
                second_ori = COMP(aln->flags) ? '-' : '+';
                a_is_base_sequence = 1;
              } 
            else 
              { assert(b_right_match && a_left_match);
                first_name = b_name;
                second_name = a_name;
                first_ori = COMP(aln->flags) ? '-' : '+';
                second_ori = '+';
                a_is_base_sequence = 0;
              }  

            fprintf(stdout, "L\t%s\t%c\t%s\t%c\t", first_name,  first_ori,
                                                   second_name, second_ori);
          }

        // Use the partial trace to reconstruct an alignment between the reads
        Compute_Trace_PTS(aln,work,tspace);

        //  Calculate the CIGAR string from the trace and write it directly to stdout
        Print_Cigar(stdout,aln, a_is_base_sequence);
        fprintf(stdout, "\n");
      }

    free(trace);
    if (ALIGN)
      { free(aln->bseq-1);
        free(aln->aseq-1);
        Free_Work_Data(work);
        Close_DB(db);
      }
  }

  fclose(input);
  exit (0);
}

