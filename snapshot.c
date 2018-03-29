/*!
 *    \file  snapshot.c
 *   \brief  This file have not thread-safe properties.
 *
 *  <+DETAILED+>
 *
 *  \author  KIM Hyeok (kh), ekh0324@gmail.com
 *
 *  \internal
 *       Created:  2018년 03월 29일
 *      Revision:  none
 *      Compiler:  gcc
 *  Organization:  Konkuk University
 *     Copyright:  Copyright (c) 2018, KIM Hyeok
 *
 *  This source code is released for free distribution under the terms of the
 *  GNU General Public License as published by the Free Software Foundation.
 */


#include "snapshot.h"
#include <stdbool.h>

#define COPY(type) int copy_ ## type (type * a, type * b) {    \
	if (a==NULL || b==NULL) return 255; \
	memcpy(a,b,sizeof(type));           \
	return 0;                           \
}
void *safe_malloc(size_t n)
{
	void *p = malloc(n);
	if (p == NULL) {
		fprintf(stderr, "Fatal: failed to allocate %zu bytes.\n", n);
		abort();
	}
	return p;
}
void safe_free(void* ptr) 
{
	if(ptr != NULL) {free(ptr);}
}

COPY(atom);
COPY(Box3);

static const char delimeter[] = " ";
static const char s_atoms[]    = "ITEM: ATOMS id type xu yu zu";
static const char s_atoms_all[]    = "ITEM: ATOMS id type xu yu zu mux muy muz vx vy vz";
static const char s_atoms_dipole[]    = "ITEM: ATOMS id type xu yu zu mux muy muz";
static const char s_atoms_vel[]    = "ITEM: ATOMS id type xu yu zu vx vy vz";
static const char s_box_bounds[] = "ITEM: BOX BOUNDS pp pp pp";
static const char s_box_bounds_xy[] = "ITEM: BOX BOUNDS ff ff pp";
static const char s_box_bounds_xyz[] = "ITEM: BOX BOUNDS ff ff ff";
static const char s_box_bounds_z[] = "ITEM: BOX BOUNDS pp pp ff";
static const char s_n_atoms[] = "ITEM: NUMBER OF ATOMS";
static const char s_timestep[] = "ITEM: TIMESTEP";
static int i_atoms ;
static int i_atoms_all ;
static int i_atoms_dipole ;
static int i_atoms_vel ;
static int i_box_bounds ;
static int i_n_atoms ;
static int i_timestep ;

static unsigned int n_atoms;
static unsigned int n_snap;
static Box3 static_box ;
static int datatype;

char * strtok_null_safe ( char * str, const char * delimiters );
int read_dump_OnlyCheck( FILE* fp) {
	i_atoms = strlen(s_atoms);
	i_atoms_all = strlen(s_atoms_all);
	i_atoms_dipole = strlen(s_atoms_dipole);
	i_atoms_vel = strlen(s_atoms_vel);
	i_box_bounds = strlen(s_box_bounds);
	i_n_atoms = strlen(s_n_atoms);
	i_timestep = strlen(s_timestep);

	bigint timestep;
//	int n_atoms;
	real xlow,xhigh,ylow,yhigh,zlow,zhigh;
	error_code =0;
/* 	atom_column atom_ids = { 
 * 		.id = -1 ,.type=-1,
 * 	 .x= -1,.y=-1,.z=-1,
 * 	 .mux=-1,.muy=-1,.muz=-1,
 * 	 .vx = -1,.vy=-1,.vz =-1
 * 	};
 */
	int id,type;
	real xu,yu,zu;
	real mux,muy,muz;
	real vx,vy,vz;
	int i;
	atom* p_atom;
#define FAIL false
#define SUCCESS n_atoms
	read_lines(1,fp);
	if( strncmp(s_timestep,line,i_timestep) !=0) {
		return FAIL;
	}

	read_lines(1,fp);
	timestep = atol(line);
#ifndef NDEBUG
	fprintf(stderr,"atol(line) = %ld\n"
			"atoi(line) = %d\n"
			, timestep,atoi(line));
#endif

	read_lines(1,fp);
	if( strncmp(s_n_atoms,line,i_n_atoms) !=0)
		return FAIL;

	read_lines(1,fp);
	n_atoms = atoi(line);
//	int  pbc[3]={0,0,0};
	int  pbcTYPE=0;
	read_lines(1,fp);
	if( strncmp(s_box_bounds,line,i_box_bounds) ==0) {
		pbcTYPE = PBC_X|PBC_Y|PBC_Z;
//		pbc[0] = 1; pbc[1]=1; pbc[2]=1;
	}
	else if( strncmp(s_box_bounds_z,line,i_box_bounds) ==0) {
		pbcTYPE = PBC_X|PBC_Y;
//		pbc[0] = 1; pbc[1]=1; pbc[2]=0;
	}
	else if( strncmp(s_box_bounds_xy,line,i_box_bounds) ==0) {
		pbcTYPE = PBC_Z;
//		pbc[0] = 0; pbc[1]=0; pbc[2]=1;
	}
	else if( strncmp(s_box_bounds_xyz,line,i_box_bounds) ==0) {
		 pbcTYPE =0;
	}
	else
		return FAIL;

	read_lines(1,fp); 
	xlow = atof(strtok_null_safe(line,delimeter));
	xhigh = atof(strtok_null_safe(NULL,delimeter));
	read_lines(1,fp); 
	ylow = atof(strtok_null_safe(line,delimeter));
	yhigh = atof(strtok_null_safe(NULL,delimeter));
	read_lines(1,fp); 
	zlow = atof(strtok_null_safe(line,delimeter));
	zhigh = atof(strtok_null_safe(NULL,delimeter));

	
	Snapshot *snap   = new_Snapshot(timestep,n_atoms);
	Box3 box = {xlow,xhigh,ylow,yhigh,zlow,zhigh, pbcTYPE};
//	Box3 box = {xlow,xhigh,ylow,yhigh,zlow,zhigh, {pbc[0],pbc[1],pbc[2]}};
	copy_Box3( &snap->box,  &box);
	copy_Box3( &static_box,  &box);

//	error( (char*)s_timestep);
//	fprintf(stderr,"%ld\n", timestep);
//	error((char*)s_n_atoms);
//	fprintf(stderr,"%d\n", snap->n_atoms);
//	error((char*)s_box_bounds);
//	fprintf(stderr,"%f %f\n", snap->box.xlow,snap->box.xhigh);
//	fprintf(stderr,"%f %f\n", snap->box.ylow,snap->box.yhigh);
//	fprintf(stderr,"%f %f\n", snap->box.zlow,snap->box.zhigh);


	read_lines(1,fp);

	if( strncmp(s_atoms,line,i_atoms) ==0) {
		datatype=0;
		for (i=0; i<n_atoms; i++){
			read_lines(1,fp);
			/*	fprintf(stderr,"%d %d %f %f %f %f %f %f\n", 
					id,type,
					xu,yu,zu,
					mux,muy,muz);*/
			if (error_code == 1) return FAIL;
		}
	}
	else if (strncmp(s_atoms_vel,line,i_atoms_vel) ==0 ) {
		datatype= FLAG_VEL ;
		for (i=0; i<n_atoms; i++){
			read_lines(1,fp);
			if (error_code == 1) return FAIL;
			/*	fprintf(stderr,"%d %d %f %f %f %f %f %f\n", 
					id,type,
					xu,yu,zu,
					mux,muy,muz);*/
		}
	}
	else if (strncmp(s_atoms_dipole,line,i_atoms_dipole) ==0 ) {
		datatype= FLAG_DIPOLE ;
		for (i=0; i<n_atoms; i++){
			read_lines(1,fp);
			if (error_code == 1) return FAIL;
			/*	fprintf(stderr,"%d %d %f %f %f %f %f %f\n", 
					id,type,
					xu,yu,zu,
					mux,muy,muz);*/
		}
	}
	else if (strncmp(s_atoms_all,line,i_atoms_all) ==0 ) {
		datatype= FLAG_VEL|FLAG_DIPOLE ;
		for (i=0; i<n_atoms; i++){
			read_lines(1,fp);
			if (error_code == 1) return FAIL;
			/*	fprintf(stderr,"%d %d %f %f %f %f %f %f\n", 
					id,type,
					xu,yu,zu,
					mux,muy,muz);*/
		}
	}
	else  {
		return FAIL;
	}

	return SUCCESS;
}
#define GET_COMMON \
			read_lines(1,fp);\
			if (error_code == 1) return FAIL;\
			id =atoi(strtok_null_safe(line, delimeter));\
			type = atoi(strtok_null_safe(NULL,delimeter));\
			xu   = atof(strtok_null_safe(NULL,delimeter));\
			yu   = atof(strtok_null_safe(NULL,delimeter));\
			zu   = atof(strtok_null_safe(NULL,delimeter));
#define GET_DIPOLE\
			mux   = atof(strtok_null_safe(NULL,delimeter));\
			muy   = atof(strtok_null_safe(NULL,delimeter));\
			muz   = atof(strtok_null_safe(NULL,delimeter));
#define GET_VELOCITY\
			vx   = atof(strtok_null_safe(NULL,delimeter));\
			vy   = atof(strtok_null_safe(NULL,delimeter));\
			vz   = atof(strtok_null_safe(NULL,delimeter));
Snapshot* read_dump( FILE* fp) {
	bigint timestep;
//	int n_atoms;
	real xlow,xhigh,ylow,yhigh,zlow,zhigh;
	int id,type;
	real xu,yu,zu,mux,muy,muz;
	real vx,vy,vz;
	int i;
	atom* p_atom;
	error_code =0;

	read_lines(1,fp);
	if( strncmp(s_timestep,line,i_timestep) !=0) {
		return (Snapshot*)(error("not ITEM: TIMESTEP"));
	}

	read_lines(1,fp);
	timestep = atol(line);
#ifndef NDEBUG
	fprintf(stderr,"atol(line) = %ld\n"
			"atoi(line) = %d\n"
			, timestep,atoi(line));
#endif

	read_lines(1,fp);
	if( strncmp(s_n_atoms,line,i_n_atoms) !=0)
		return (Snapshot*)(error("not ITEM: NUMBER OF ATOMS"));

	read_lines(1,fp);
	n_atoms = atoi(line);
	int  pbcTYPE=0;
	read_lines(1,fp);
	if( strncmp(s_box_bounds,line,i_box_bounds) ==0) {
		pbcTYPE = PBC_X|PBC_Y|PBC_Z;
//		pbc[0] = 1; pbc[1]=1; pbc[2]=1;
	}
	else if( strncmp(s_box_bounds_z,line,i_box_bounds) ==0) {
		pbcTYPE = PBC_X|PBC_Y;
//		pbc[0] = 1; pbc[1]=1; pbc[2]=0;
	}
	else if( strncmp(s_box_bounds_xy,line,i_box_bounds) ==0) {
		pbcTYPE = PBC_Z;
//		pbc[0] = 0; pbc[1]=0; pbc[2]=1;
	}
	else if( strncmp(s_box_bounds_xyz,line,i_box_bounds) ==0) {
		 pbcTYPE =0;
	}
	else
		return (Snapshot*)(error("not ITEM: BOX BOUNDS pp pp pp(ff)"));

	read_lines(1,fp); 
	xlow = atof(strtok_null_safe(line,delimeter));
	xhigh = atof(strtok_null_safe(NULL,delimeter));
	read_lines(1,fp); 
	ylow = atof(strtok_null_safe(line,delimeter));
	yhigh = atof(strtok_null_safe(NULL,delimeter));
	read_lines(1,fp); 
	zlow = atof(strtok_null_safe(line,delimeter));
	zhigh = atof(strtok_null_safe(NULL,delimeter));

	
	Snapshot *snap   = new_Snapshot(timestep,n_atoms);
//	Box3 box = {xlow,xhigh,ylow,yhigh,zlow,zhigh, {pbc[0],pbc[1],pbc[2]}};
	Box3 box = {xlow,xhigh,ylow,yhigh,zlow,zhigh, pbcTYPE};
	copy_Box3( &snap->box,  &box);
	copy_Box3( &static_box,  &box);

	error( (char*)s_timestep);
//	fprintf(stderr,"%ld\n", timestep);
	error((char*)s_n_atoms);
//	fprintf(stderr,"%d\n", snap->n_atoms);
	error((char*)s_box_bounds);
//	fprintf(stderr,"%f %f\n", snap->box.xlow,snap->box.xhigh);
//	fprintf(stderr,"%f %f\n", snap->box.ylow,snap->box.yhigh);
//	fprintf(stderr,"%f %f\n", snap->box.zlow,snap->box.zhigh);

	read_lines(1,fp);
	
	if( strncmp(s_atoms_dipole,line,i_atoms_dipole) ==0) {
#ifndef NDEBUG
		fprintf(stderr,"s_atoms_dipole style snapshot\n");
#endif
		for (i=0; i<n_atoms; i++){
			GET_COMMON;
			GET_DIPOLE;
			p_atom = &(snap->atoms[i]);
			make_atom_dipole( p_atom,id,type,xu,yu,zu,mux,muy,muz);
		}
	}
	else if (strncmp(s_atoms_vel,line,i_atoms_vel) ==0 ) {
#ifndef NDEBUG
		fprintf(stderr,"s_atoms_vel style snapshot\n");
#endif
		for (i=0; i<n_atoms; i++){
			GET_COMMON ; 
			GET_VELOCITY ;

			p_atom = &(snap->atoms[i]);
			make_atom_vel( p_atom,id,type,xu,yu,zu,vx,vy,vz);
			/*	fprintf(stderr,"%d %d %f %f %f %f %f %f\n", 
					id,type,
					xu,yu,zu,
					mux,muy,muz);*/
		}

	}
	else if (strncmp(s_atoms_all,line,i_atoms_all) ==0 ) {
#ifndef NDEBUG
		fprintf(stderr,"s_atoms_all style snapshot\n");
#endif
		for (i=0; i<n_atoms; i++){
			GET_COMMON;
			GET_DIPOLE;
			GET_VELOCITY;

			p_atom = &(snap->atoms[i]);
			make_atom_all( p_atom,id,type,xu,yu,zu,mux,muy,muz,vx,vy,vz);
		}

	}
	else if (strncmp(s_atoms,line,i_atoms) ==0 ) {
		for (i=0; i<n_atoms; i++){
			GET_COMMON ;

			p_atom = &(snap->atoms[i]);
			make_atom( p_atom,id,type,xu,yu,zu);
			/*	fprintf(stderr,"%d %d %f %f %f %f %f %f\n", 
					id,type,
					xu,yu,zu,
					mux,muy,muz);*/
		}

	}
	else  {
		return (Snapshot*)(error("not ITEM: ATOMS id type xu yu zu mux muy muz"));
	}

	error_code =0;
	return snap;
}
void* error( char string[MAXLINE] ) {
#ifndef NDEBUG
	fputs( string, stderr );
	fputs( "\n", stderr );
#endif
	return NULL;
	//	exit(1);
}
char * strtok_null_safe ( char * str, const char * delimiters ) {

	char* ret_str = strtok(str,delimiters);
	if (ret_str == NULL ) {
		error_code = 1 ; 
		return "0";
	}
	return  ret_str;
}

void read_lines(int n,FILE* fp)  // from lammps reader_native.cpp
{
	char *eof;int i;
	for (i = 0; i < n; i++) eof = fgets(line,MAXLINE,fp);
	if (eof == NULL) {
		error_code = 1;
		error("Unexpected end of dump file");
	}
}


Snapshot* new_Snapshot(bigint timestep, int n) {
	Snapshot* snap = (Snapshot*) malloc(sizeof(Snapshot));
	snap->timestep = timestep;
	snap-> n_atoms = n;
	snap->atoms = (atom*) malloc(sizeof(atom)*n);
	return snap;
}


void free_Snapshot(Snapshot* snap) {
 if(snap->atoms !=NULL)
	 free(snap->atoms);
 if(snap !=NULL)
	 free(snap);
}

int make_atom(atom* col,int id, int type, 
		real x, real y, real z) {
	real mu1;
	if (col ==NULL)
		return 255;
	col->id=id; col->type=type;
	col->x=x;col->y=y;col->z=z;
	col->atomType = 0;
	return 0;
}
int make_atom_dipole(atom* col,int id, int type, 
		real x, real y, real z,
		real mux,real muy, real muz) {
	real mu1;
	if (col ==NULL)
		return 255;
	col->id=id; col->type=type;
	col->x=x;col->y=y;col->z=z;
	mu1 = sqrt(mux*mux+muy*muy+muz*muz);
	col->mu1 = mu1;

	if (mu1>0.001) { 
		col->mux=mux/mu1;col->muy=muy/mu1;col->muz=muz/mu1;
	}
	col->atomType = FLAG_DIPOLE;
	return 0;
}
int make_atom_vel(atom* col,int id, int type, 
		real x, real y, real z,
		real vx,real vy, real vz) {
	real mu1;
	if (col ==NULL)
		return 255;
	col->id=id; col->type=type;
	col->x=x;col->y=y;col->z=z;
	col->vx=vx; 
	col->vy=vy;
	col->vz=vz;
	col->atomType = FLAG_VEL;
	return 0;
}
int make_atom_all(atom* col,int id, int type, 
		real x, real y, real z,
		real mux,real muy, real muz,
		real vx, real vy, real vz) {
	real mu1;
	if (col ==NULL)
		return 255;
	col->id=id; col->type=type;
	col->x=x;col->y=y;col->z=z;
	mu1 = sqrt(mux*mux+muy*muy+muz*muz);
	col->mu1 = mu1;
	col->vx=vx; 
	col->vy=vy;
	col->vz=vz;
	col->atomType = FLAG_VEL|FLAG_DIPOLE;

	if (mu1>0.001) { 
		col->mux=mux/mu1;col->muy=muy/mu1;col->muz=muz/mu1;
	}
	return 0;
}
int dump_stream(atomstream* stream, FILE* fp, int nTime, int n_atoms,int s_id, int s_type) 
{
	int id,type;
	real xu,yu,zu,mux,muy,muz;
	real vx,vy,vz;
	atom* p_atom;
	fseek(fp, 0, SEEK_SET); // 
	int i;
	
	for( i=0; i<nTime; i++) {
		read_lines(9,fp);

		if( strncmp(s_atoms_dipole,line,i_atoms_dipole) ==0) {
			for (i=0; i<n_atoms; i++){
				read_lines(1,fp);
				if (error_code == 1) return FAIL;
				id =atoi(strtok_null_safe(line, delimeter));
				type = atoi(strtok_null_safe(NULL,delimeter));
				if ( s_id==id && s_type == type) {
					xu   = atof(strtok_null_safe(NULL,delimeter));
					yu   = atof(strtok_null_safe(NULL,delimeter));
					zu   = atof(strtok_null_safe(NULL,delimeter));
					mux   = atof(strtok_null_safe(NULL,delimeter));
					muy   = atof(strtok_null_safe(NULL,delimeter));
					muz   = atof(strtok_null_safe(NULL,delimeter));

					stream->x[i] = xu;
					stream->y[i] = yu;
					stream->z[i] = zu;
					stream->mux[i] = mux;
					stream->muy[i] = muy;
					stream->muz[i] = muz;
					stream->atomType = FLAG_DIPOLE;
				}

			}
		}
		else if (strncmp(s_atoms_vel,line,i_atoms_vel) ==0 ) {
			for (i=0; i<n_atoms; i++){
				read_lines(1,fp);
				if (error_code == 1) return FAIL;
				id =atoi(strtok_null_safe(line, delimeter));
				type = atoi(strtok_null_safe(NULL,delimeter));
				if ( s_id==id && s_type == type) {
					xu   = atof(strtok_null_safe(NULL,delimeter));
					yu   = atof(strtok_null_safe(NULL,delimeter));
					zu   = atof(strtok_null_safe(NULL,delimeter));
					vx   = atof(strtok_null_safe(NULL,delimeter));
					vy   = atof(strtok_null_safe(NULL,delimeter));
					vz   = atof(strtok_null_safe(NULL,delimeter));

					stream->x[i] = xu;
					stream->y[i] = yu;
					stream->z[i] = zu;
					stream->vx[i] = vx;
					stream->vy[i] = vy;
					stream->vz[i] = vz;
					stream->atomType = FLAG_VEL;
				}

			}
		}
		else if (strncmp(s_atoms_all,line,i_atoms_all) ==0 ) {
			for (i=0; i<n_atoms; i++){
				read_lines(1,fp);
			if (error_code == 1) return FAIL;
				id =atoi(strtok_null_safe(line, delimeter));
				type = atoi(strtok_null_safe(NULL,delimeter));
				if ( s_id==id && s_type == type) {
					xu   = atof(strtok_null_safe(NULL,delimeter));
					yu   = atof(strtok_null_safe(NULL,delimeter));
					zu   = atof(strtok_null_safe(NULL,delimeter));
					mux   = atof(strtok_null_safe(NULL,delimeter));
					muy   = atof(strtok_null_safe(NULL,delimeter));
					muz   = atof(strtok_null_safe(NULL,delimeter));
					vx   = atof(strtok_null_safe(NULL,delimeter));
					vy   = atof(strtok_null_safe(NULL,delimeter));
					vz   = atof(strtok_null_safe(NULL,delimeter));

					stream->x[i] = xu;
					stream->y[i] = yu;
					stream->z[i] = zu;
					stream->vx[i] = vx;
					stream->vy[i] = vy;
					stream->vz[i] = vz;
					stream->mux[i] = mux;
					stream->muy[i] = muy;
					stream->muz[i] = muz;
					stream->atomType = FLAG_VEL|FLAG_DIPOLE;
				}

			}
		}
		else if (strncmp(s_atoms,line,i_atoms) ==0 ) {
			for (i=0; i<n_atoms; i++){
				read_lines(1,fp);
			if (error_code == 1) return FAIL;
				id =atoi(strtok_null_safe(line, delimeter));
				type = atoi(strtok_null_safe(NULL,delimeter));
				if ( s_id==id && s_type == type) {
					xu   = atof(strtok_null_safe(NULL,delimeter));
					yu   = atof(strtok_null_safe(NULL,delimeter));
					zu   = atof(strtok_null_safe(NULL,delimeter));
					stream->x[i] = xu;
					stream->y[i] = yu;
					stream->z[i] = zu;
					stream->atomType = 0;
				}
			}
		}
	}
	return 0;
}
int malloc_stream( atomstream* stream) 
{
	stream = (atomstream*) safe_malloc( sizeof(atomstream));
	stream->x = (real*) safe_malloc( nTime* sizeof(real));
	stream->y = (real*) safe_malloc( nTime* sizeof(real));
	stream->z = (real*) safe_malloc( nTime* sizeof(real));
	if ((datatype & FLAG_DIPOLE) == FLAG_DIPOLE) {
		stream->mux = (real*) safe_malloc( nTime* sizeof(real));
		stream->muy = (real*) safe_malloc( nTime* sizeof(real));
		stream->muz = (real*) safe_malloc( nTime* sizeof(real));
	}
	if ((datatype & FLAG_VELOCITY) == FLAG_VELOCITY) {
		stream->vx = (real*) safe_malloc( nTime* sizeof(real));
		stream->vy = (real*) safe_malloc( nTime* sizeof(real));
		stream->vz = (real*) safe_malloc( nTime* sizeof(real));
	}
	stream->nTime = nTime;
}
int free_stream( atomstream* stream)
{
	safe_free(&(stream->x));
	safe_free(&(stream->y));
	safe_free(&(stream->z));
	safe_free(&(stream->mux));
	safe_free(&(stream->muy));
	safe_free(&(stream->muz));
	safe_free(&(stream->mux));
	safe_free(&(stream->muy));
	safe_free(&(stream->muz));

	safe_free(stream);
}
int free_atomstreamlist(Atomstream_list * list)
{
	int i;
	int n_atoms = list->n_atoms;
	
	for ( i = 0; i < n_atoms; i += 1 ) {
		free_stream(&list->streams[i]);
	}
	safe_free(list);
}
Atomstream_list *  new_atomstreamlist() 
{
	Atomstream_list * list;
	malloc_atomstreamlist(list);
	return list;
}
#define SET_STREAM_COMMON\
				stream[in]->x[it] = xu;\
				stream[in]->y[it] = yu;\
				stream[in]->z[it] = zu;
#define SET_STREAM_VELOCITY\
				stream[in]->vx[it] = vx;\
				stream[in]->vy[it] = vy;\
				stream[in]->vz[it] = vz;
#define SET_STREAM_DIPOLE\
				stream[in]->mux[it] = mux;\
				stream[in]->muy[it] = muy;\
				stream[in]->muz[it] = muz;

int malloc_atomstreamlist(Atomstream_list *list )
{
	list = (Atomstream_list*) safe_malloc(sizeof(Atomstream_list));
	list->n_atoms= n_atoms;
	list->streams = (atomstream*) safe_malloc(
			n_atoms*(sizeof(atomstream)));
	
	for ( int i = 0; i < n_atoms; i += 1 ) {
		malloc_stream(& list->streams[i]);
	}
	copy_Box3(&list->box, &static_box);
	return 0;
}
int fill_atomstreamlist(Atomstream_list* list, FILE*)
{
	int id,type;
	real xu,yu,zu,mux,muy,muz;
	real vx,vy,vz;
	atom* p_atom;
	atomstream* stream = list->streams;
	rewind(fp);
	int it,in;
	
	for( it=0; it<nTitme; it++) {
		read_lines(9,fp);
		if( strncmp(s_atoms_dipole,line,i_atoms_dipole) ==0) {
			for (in=0; in<n_atoms; in++){
				GET_COMMON;
				GET_DIPOLE;

				SET_STREAM_COMMON;
				SET_STREAM_DIPOLE;
				stream[in]->atomType = FLAG_DIPOLE;

			}
		}
		else if (strncmp(s_atoms_vel,line,i_atoms_vel) ==0 ) {
			for (in=0; in<n_atoms; in++){
				GET_COMMON;
				GET_DIPOLE;
				GET_VELOCITY;

				SET_STREAM_COMMON;
				SET_STREAM_DIPOLE;
				SET_STREAM_VELOCITY;
				stream[in]->atomType = FLAG_VEL;

			}
		}
		else if (strncmp(s_atoms_all,line,i_atoms_all) ==0 ) {
			for (in=0; in<n_atoms; in++){
				GET_COMMON;
				GET_DIPOLE;
				GET_VELOCITY;

				SET_STREAM_COMMON;
				SET_STREAM_DIPOLE;
				SET_STREAM_VELOCITY;
				stream[in]->atomType = FLAG_VEL|FLAG_DIPOLE;
			}
		}
		else if (strncmp(s_atoms,line,i_atoms) ==0 ) {
			for (in=0; in<n_atoms; in++){
				GET_COMMON;

				SET_STREAM_COMMON;
				stream[in]->atomType = 0;
			}
		}
	}
	return 0;
}

int Make_Info(FILE* fp) 
{
	n_snap = 0;	

	rewind(fp);
	while(1) {
		bool check =	read_dump_OnlyCheck(fp);
		if (check == false )
			break;
		n_snap++; 
	}
	return 0;
}
unsigned int get_number_of_snapshots()
{
	return n_snap;
}
unsigned int get_number_of_atoms()
{
	return n_atoms;
}
Box3 get_box() 
{
	return static_box;
}
void fwrite_streamlist(FILE * fp,Atomstream_list* streamlist)
{
	fwrite(&streamlist->n_atoms,sizeof(n_atoms),1,fp);
	fwrite(&streamlist->nTime,sizeof(n_atoms),1,fp);
	fwrite(&streamlist->atomType,sizeof(n_atoms),1,fp);
	atomstream* stream = streamlist->streams;
	
	for ( in = 0; in < n_atoms; in += 1 ) {
		fwrite(&stream[in].x, sizeof(real), nTime,fp);
	}
	for ( in = 0; in < n_atoms; in += 1 ) {
		fwrite(&stream[in].y, sizeof(real), nTime,fp);
	}
	for ( in = 0; in < n_atoms; in += 1 ) {
		fwrite(&stream[in].z, sizeof(real), nTime,fp);
	}
	if ((datatype & FLAG_DIPOLE) == FLAG_DIPOLE) {
		for ( in = 0; in < n_atoms; in += 1 ) {
			fwrite(&stream[in].mux, sizeof(real), nTime,fp);
		}
		for ( in = 0; in < n_atoms; in += 1 ) {
			fwrite(&stream[in].muy, sizeof(real), nTime,fp);
		}
		for ( in = 0; in < n_atoms; in += 1 ) {
			fwrite(&stream[in].muz, sizeof(real), nTime,fp);
		}
	}
	if ((datatype & FLAG_VELOCITY) == FLAG_VELOCITY) {
		for ( in = 0; in < n_atoms; in += 1 ) {
			fwrite(&stream[in].vx, sizeof(real), nTime,fp);
		}
		for ( in = 0; in < n_atoms; in += 1 ) {
			fwrite(&stream[in].vy, sizeof(real), nTime,fp);
		}
		for ( in = 0; in < n_atoms; in += 1 ) {
			fwrite(&stream[in].vz, sizeof(real), nTime,fp);
		}
	}
	
}
Atomstream_list* fread_streamlist(FILE * fp)
{
	Atomstream_list* streamlist ;

	
	fread(&n_atoms,sizeof(n_atoms),1,fp);
	fread(&nTime,sizeof(n_atoms),1,fp);
	fread(&atomType,sizeof(n_atoms),1,fp);

	malloc_atomstreamlist(streamlist);

	atomstream* stream = streamlist->streams;
	
	for ( in = 0; in < n_atoms; in += 1 ) {
		fread(&stream[in].x, sizeof(real), nTime,fp);
	}
	for ( in = 0; in < n_atoms; in += 1 ) {
		fread(&stream[in].y, sizeof(real), nTime,fp);
	}
	for ( in = 0; in < n_atoms; in += 1 ) {
		fread(&stream[in].z, sizeof(real), nTime,fp);
	}
	if ((datatype & FLAG_DIPOLE) == FLAG_DIPOLE) {
		for ( in = 0; in < n_atoms; in += 1 ) {
			fread(&stream[in].mux, sizeof(real), nTime,fp);
		}
		for ( in = 0; in < n_atoms; in += 1 ) {
			fread(&stream[in].muy, sizeof(real), nTime,fp);
		}
		for ( in = 0; in < n_atoms; in += 1 ) {
			fread(&stream[in].muz, sizeof(real), nTime,fp);
		}
	}
	if ((datatype & FLAG_VELOCITY) == FLAG_VELOCITY) {
		for ( in = 0; in < n_atoms; in += 1 ) {
			fread(&stream[in].vx, sizeof(real), nTime,fp);
		}
		for ( in = 0; in < n_atoms; in += 1 ) {
			fread(&stream[in].vy, sizeof(real), nTime,fp);
		}
		for ( in = 0; in < n_atoms; in += 1 ) {
			fread(&stream[in].vz, sizeof(real), nTime,fp);
		}
	}
}
#undef FAIL
#undef SUCCESS
