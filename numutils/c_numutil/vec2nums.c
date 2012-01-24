/************************************************************************/
/* vec2nums.c - 							*/
/* This program converts real or complex vector files to simple number  */
/* list files, such as used by Matlab, by removing the header line and  */
/* index column.							*/
/* 									*/
/* Example Usage:							*/
/*	vec2nums vecfile  newfile  				  	*/
/* 								  	*/
/* File Types:								*/
/*      In:   Real or complex Vector.					*/
/*      Out:  (textual list of numbers)					*/
/*									*/
/* From:  www.atl.lmco.com/proj/csim/xgraph/numutil                     */
/*	  chein@atl.lmco.com						*/
/*									*/
/* Copyright (C) 2002 CHein						*/
/* This program is free software; you can redistribute it and/or modify	*/
/* it under the terms of the GNU General Public License as published by	*/
/* the Free Software Foundation, Version 2 of the GPL.			*/
/* This program is distributed in the hope that it will be useful,	*/
/* but WITHOUT ANY WARRANTY; without even the implied warranty of	*/
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	*/
/* GNU General Public License for more details.				*/
/* You should have received a copy of the GNU General Public License	*/
/* along with this program; if not, write to the Free Software		*/
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307.	*/
/************************************************************************/

#include "commonlibs.c"
 

main (argc, argv)
int argc;
char *argv[];
{
 char fname1[999], fname2[999]="", line[500], word[100], trsh[100];
 int X, Y, csw=0;
 double x,y;
 complex a, b;
 int i, j, k, arg;
 struct file_kind *infile=0;
 FILE *outfile=0;

 arg = 1;  k = 0;
 while (argc>arg)
 { /*arg*/

  if (argv[arg][0]=='-')
   { /*option*/
     {printf("ERROR vec2nums: Unknown command-line option %s\n", argv[arg]); exit(1);}
   } /*option*/
  else
   { /*file_name*/
    if (k==0)  
     { 
      strcpy(fname1,argv[arg]);
      infile = file_open_read( fname1 );
      if (infile->fpt == 0) 
       {printf("ERROR vec2nums: Cannot open input file '%s'.\n",fname1); exit(1);}
      k = k + 1;
     }
    else
    if (k==1)
     { strcpy(fname2,argv[arg]);  
       if (strcmp(fname1,fname2)==0) {printf("ERROR vec2nums: Can't over-write source.\n"); exit(1);}
       k = k + 1; 
     }
    else
     {printf("ERROR vec2nums: Too many file-names on command line.\n"); exit(1);}
   } /*file_name*/

  arg = arg + 1;
 } /*arg*/

 if (k==0) { printf("ERROR vec2nums: No input file on command line.\n"); exit(1); }

 if (k<2)  /* No output file specified, so make default name by appending suffix. */
  { 
    strcpy(fname2,fname1);
    i = strlen(fname2)-1;  while ((i>=0) && (fname2[i]!='.')) i = i-1;
    if (i<0) i = strlen(fname2);
    fname2[i] = '\0';
    strcat(fname2,".nums");
    if (strcmp(fname1,fname2)==0) {printf("ERROR vec2nums: Can't over-write source.\n"); exit(1);}
  }

 printf("	(Writing output to %s)\n", fname2 );
 outfile = fopen( fname2, "w" );
 if (outfile==0) {printf("ERROR vec2nums: Can't open output file '%s'.\n",fname2); exit(1);}

 k = 0;
 Read_Line( infile->fpt, line, 500 );
 while (! feof( infile->fpt ))
  {
   Next_Word(line,trsh,word," 	,()");
   if (word[0]!='\0')
    {
     if (infile->kind=='R')
      {
	if (infile->dim2>1) Next_Word(line,trsh,word," 	,()");
	Next_Word(line,trsh,word," 	,()");
       	fprintf(outfile,"%s\n", word);
      }
     else
      {
	fprintf(outfile,"%s	", word);
	Next_Word(line,trsh,word,"      ,()");
	fprintf(outfile,"%s\n", word);
      }
     k = k + 1;
    }
   i = i + 1;
   Read_Line( infile->fpt, line, 500 );
  }
 fclose(infile->fpt);
 fclose(outfile);
 printf("Read %d points\n", k);
}