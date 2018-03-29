/*!
 *    \file  main.c
 *   \brief  
 *
 *  <+DETAILED+>
 *
 *  \author  KIM Hyeok (kh), ekh0324@gmail.com
 *
 *  \internal
 *       Created:  2018년 03월 28일
 *      Revision:  none
 *      Compiler:  gcc
 *  Organization:  Konkuk University
 *     Copyright:  Copyright (c) 2018, KIM Hyeok
 *
 *  This source code is released for free distribution under the terms of the
 *  GNU General Public License as published by the Free Software Foundation.
 */



#include <assert.h>

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include	<unistd.h>
#include <stdbool.h>
#define M_PI       3.14159265358979323846
#include<math.h>
#include"snapshot.h"
#include <stdio.h>
#include <fcntl.h>	/* for open flags */
#include <limits.h>	/* for PATH_MAX */


#include	<stdlib.h>

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  
 * =====================================================================================
 */
static char template[] = "/tmp/myfileXXXXXX";
	int
main ( int argc, char *argv[] )
{
	char filename[100];
	strcpy( filename,argv[1]);
	Atomstream_list* streamlist;

	//		FILE* fp = fopen( filename ,"r");
	if( fp == NULL) {
		//			files_on [opt_num] = false;
		fprintf(stderr,"Can`t open file (%s)!!\n", filename);
		exit(EXIT_FAILURE);
	}


	Make_Info(fp);

	malloc_atomstreamlist(streamlist);
	fill_atomstreamlist(streamlist,fp);

	
	static char template[] = "/tmp/myfileXXXXXX";
	strcpy(fname, template);		/* Copy template */
	char* outname = tmpnam( NULL);
	printf("Temporary output filename : %s \n",outname);

	FILE* ofp = fopen(outname,"w");


	

	free_atomstreamlist(streamlist);
	fclose (fp);
	fclose(ofp);

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */


