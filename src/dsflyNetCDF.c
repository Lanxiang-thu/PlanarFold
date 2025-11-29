/* Copyright 2019 University Corporation for Atmospheric
   Research/Unidata.  See COPYRIGHT file for conditions of use. */
#include <stdio.h>
#include <string.h>
//#include <iostream>
#include "netcdf.h"
#include "dsfly.h"
 
/* Dimensions. */
#define NSPATIAL 3          // spatial
#define NCELLSPATIAL 3      // cell_spatial
#define NLABLE 5            // cell_spatial
#define NCELLANGULAR 3      // cell_spatial
/// #define NATOM
/// #define NFRAME

/* Names of things. */
#define FRAME_NAME "frame"
#define ATOM_NAME "atom"
#define TIME_NAME "time"
#define LABLE_NAME "lable"
#define SPAT_NAME "spatial"
#define CRD_NAME "coordinates"
#define VEL_NAME "velocities"
#define CSPAT_NAME "cell_spatial"
#define CANGU_NAME "cell_angular"
#define CLENS_NAME "cell_lengths"
#define CANGLS_NAME "cell_angles"

/* Handle errors by printing an error message  */
#define ERR(e) {printf("Error: %s\n", nc_strerror(e));}

//----------- C ---------- R e s t a r t --------------------------------------
void createRst(int natom, ncfile *ncf){
    int spat_dimid, cSpat_dimid, lable_dimid, cAngu_dimid;
    int atom_dimid; // change according to 'natom'
    int dimid[2];   // coordinates, velocities, cell_angular
    int retval;
    int NATM = natom;
    int time_dimid; //???

    printf("Entered createRst!!\n");
    /* Define the file. */
    if ((retval = nc_create(ncf->Rst_filename, NC_64BIT_OFFSET, &ncf->ncid_Rst)))   // NC_64BIT_OFFSET or NC_CLOBBER
        {ERR(retval); printf("Bug1: createRst\n");}

    /* Define the dimensions. */
    printf("Define the dimensions...\n");
    // spatial(3)
    if ((retval = nc_def_dim(ncf->ncid_Rst, SPAT_NAME, NSPATIAL, &spat_dimid)))
        {ERR(retval); printf("Bug2: createRst\n");}
    // atom(natom) !!!!!!!!
    if ((retval = nc_def_dim(ncf->ncid_Rst, ATOM_NAME, NATM, &atom_dimid)))
        {ERR(retval); printf("Bug3: createRst\n");}
    // cell_spatial(3)
    if ((retval = nc_def_dim(ncf->ncid_Rst, CSPAT_NAME, NCELLSPATIAL, &cSpat_dimid)))
        {ERR(retval); printf("Bug4: createRst\n");}
    // label(5)
    if ((retval = nc_def_dim(ncf->ncid_Rst, LABLE_NAME, NLABLE, &lable_dimid)))
        {ERR(retval); printf("Bug5: createRst\n");}
    // cell_angular(3)
    if ((retval = nc_def_dim(ncf->ncid_Rst, CANGU_NAME, NCELLANGULAR, &cAngu_dimid)))
        {ERR(retval); printf("Bug6: createRst\n");}

    /* Define the variables.*/
    printf("Define the variables...\n");
    //  float64 time()
    if ((retval = nc_def_var(ncf->ncid_Rst, TIME_NAME, NC_DOUBLE, 0, &time_dimid, &ncf->varid_time_Rst)))       // ???
        {ERR(retval); printf("Bug7: createRst\n");}
    //  |S1 spatial(spatial)
    if ((retval = nc_def_var(ncf->ncid_Rst, SPAT_NAME, NC_CHAR, 1, &spat_dimid, &ncf->varid_spatial_Rst)))
        {ERR(retval); printf("Bug8: createRst\n");}
    //  float64 coordinates(atom,spatial)
    dimid[0] = atom_dimid;
    dimid[1] = spat_dimid;
    if ((retval = nc_def_var(ncf->ncid_Rst, CRD_NAME, NC_DOUBLE, 2, dimid, &ncf->varid_crd_Rst)))
        {ERR(retval); printf("Bug9: createRst\n");}
    //  float64 velocities(atom,spatial)
    dimid[0] = atom_dimid;
    dimid[1] = spat_dimid;
    if ((retval = nc_def_var(ncf->ncid_Rst, VEL_NAME, NC_DOUBLE, 2, dimid, &ncf->varid_vel_Rst)))
        {ERR(retval); printf("Bug10: createRst\n");}
    //  |S1 cell_spatial(cell_spatial)
    if ((retval = nc_def_var(ncf->ncid_Rst, CSPAT_NAME, NC_CHAR, 1, &cSpat_dimid, &ncf->varid_cellSpatial_Rst)))
        {ERR(retval); printf("Bug11: createRst\n");}
    //  |S1 cell_angular(cell_angular,label)
    dimid[0] = cAngu_dimid;
    dimid[1] = lable_dimid;
    if ((retval = nc_def_var(ncf->ncid_Rst, CANGU_NAME, NC_CHAR, 2, dimid, &ncf->varid_cellAngular_Rst)))
        {ERR(retval); printf("Bug12: createRst\n");}
    //  float64 cell_lengths(cell_spatial)
    if ((retval = nc_def_var(ncf->ncid_Rst, CLENS_NAME, NC_DOUBLE, 1, &cSpat_dimid, &ncf->varid_cellLengths_Rst)))       // ???
        {ERR(retval); printf("Bug13: createRst\n");}
    //  float64 cell_angles(cell_angular)
    if ((retval = nc_def_var(ncf->ncid_Rst, CANGLS_NAME, NC_DOUBLE, 1, &cAngu_dimid, &ncf->varid_cellAngles_Rst)))       // ???
        {ERR(retval); printf("Bug14: createRst\n");}

    /* Write Global attributes */
    printf("Write Global attributes...\n");
    char title[13];
    strcpy(title,"default_name");
    if ((retval = nc_put_att_text(ncf->ncid_Rst, NC_GLOBAL, "title", strlen(title), title)))
        {ERR(retval); printf("Bug19: createRst\n");}
    if ((retval = nc_put_att_text(ncf->ncid_Rst, NC_GLOBAL, "application", 5, "AMBER")))
        {ERR(retval); printf("Bug20: createRst\n");}
    if ((retval = nc_put_att_text(ncf->ncid_Rst, NC_GLOBAL, "program", 9, "damselfly")))
        {ERR(retval); printf("Bug21: createRst\n");}
    if ((retval = nc_put_att_text(ncf->ncid_Rst, NC_GLOBAL, "programVersion", 3, "1.0")))
        {ERR(retval); printf("Bug22: createRst\n");}
    if ((retval = nc_put_att_text(ncf->ncid_Rst, NC_GLOBAL, "Conventions", 12, "AMBERRESTART")))    // must be 'AMBERRESTART' for rst
        {ERR(retval); printf("Bug23: createRst\n");}
    if ((retval = nc_put_att_text(ncf->ncid_Rst, NC_GLOBAL, "ConventionVersion", 3, "1.0")))
        {ERR(retval); printf("Bug24: createRst\n");}

    /* End define mode. */
    if ((retval = nc_enddef(ncf->ncid_Rst)))
        {ERR(retval); printf("Bug15: createRst\n");}


    //------------------------------------------------------
    /* Write spatial, cell_spatial, and cell_angular */
    printf("Write spatial, cell_spatial, and cell_angular...\n");
    size_t start, count, start2[2], count2[2];
    char xyz[3];
    char abc[15] = {'a', 'l', 'p', 'h', 'a',
                    'b', 'e', 't', 'a', ' ',
                    'g', 'a', 'm', 'm', 'a' };
    // spatial
    start = 0;
    count = 3;
    xyz[0] = 'x'; xyz[1] = 'y'; xyz[2] = 'z';
    if ((retval = nc_put_vara_text(ncf->ncid_Rst, ncf->varid_spatial_Rst, &start, &count, xyz)))
        {ERR(retval); printf("Bug16: createRst\n");}
    // cell_spatial
    xyz[0] = 'a'; xyz[1] = 'b'; xyz[2] = 'c';
    if ((retval = nc_put_vara_text(ncf->ncid_Rst, ncf->varid_cellSpatial_Rst, &start, &count, xyz)))
        {ERR(retval); printf("Bug17: createRst\n");}
    // cell_angular
    start2[0] = 0;
    count2[0] = 3;
    start2[1] = 0;
    count2[1] = 5;
    if ((retval = nc_put_vara_text(ncf->ncid_Rst, ncf->varid_cellAngular_Rst, start2, count2, abc)))
        {ERR(retval); printf("Bug18: createRst\n");}
}

//----------- C ---------- T r a j e c t o r y---------------------------------
void createTrj(int natom, int nframe, ncfile *ncf){
    int spat_dimid, cSpat_dimid, lable_dimid, cAngu_dimid;
    int atom_dimid, frame_dimid; // change according to 'natom' and 'nframe'
    int dimid2[2], dimid3[3];   // coordinates, velocities, cell_angular
    int retval;
    int NATM = natom, NFRAME = nframe;

    /* Define the file. */
    if ((retval = nc_create(ncf->Trj_filename, NC_64BIT_OFFSET, &ncf->ncid_Trj)))   // NC_64BIT_OFFSET or NC_CLOBBER
        ERR(retval);

    /* Define the dimensions. */
    // frame(nframe) !!!!!!!!
    if ((retval = nc_def_dim(ncf->ncid_Trj, FRAME_NAME, NFRAME, &frame_dimid)))
        ERR(retval);
    // spatial(3)
    if ((retval = nc_def_dim(ncf->ncid_Trj, SPAT_NAME, NSPATIAL, &spat_dimid)))
        ERR(retval);
    // atom(natom) !!!!!!!!
    if ((retval = nc_def_dim(ncf->ncid_Trj, ATOM_NAME, NATM, &atom_dimid)))
        ERR(retval);
    // cell_spatial(3)
    if ((retval = nc_def_dim(ncf->ncid_Trj, CSPAT_NAME, NCELLSPATIAL, &cSpat_dimid)))
        ERR(retval);
    // label(5)
    if ((retval = nc_def_dim(ncf->ncid_Trj, LABLE_NAME, NLABLE, &lable_dimid)))
        ERR(retval);
    // cell_angular(3)
    if ((retval = nc_def_dim(ncf->ncid_Trj, CANGU_NAME, NCELLANGULAR, &cAngu_dimid)))
        ERR(retval);

    /* Define the variables.*/
    //  float32 time(frame)
    if ((retval = nc_def_var(ncf->ncid_Trj, TIME_NAME, NC_FLOAT, 1, &frame_dimid, &ncf->varid_time_Trj)))       // ???
        ERR(retval);
    //  |S1 spatial(spatial)
    if ((retval = nc_def_var(ncf->ncid_Trj, SPAT_NAME, NC_CHAR, 1, &spat_dimid, &ncf->varid_spatial_Trj)))
        ERR(retval);
    //  float32 coordinates(frame,atom,spatial)
    dimid3[0] = frame_dimid;
    dimid3[1] = atom_dimid;
    dimid3[2] = spat_dimid;
    if ((retval = nc_def_var(ncf->ncid_Trj, CRD_NAME, NC_FLOAT, 3, dimid3, &ncf->varid_crd_Trj)))
        ERR(retval);
    //  |S1 cell_spatial(cell_spatial)
    if ((retval = nc_def_var(ncf->ncid_Trj, CSPAT_NAME, NC_CHAR, 1, &cSpat_dimid, &ncf->varid_cellSpatial_Trj)))
        ERR(retval);
    //  |S1 cell_angular(cell_angular,label)
    dimid2[0] = cAngu_dimid;
    dimid2[1] = lable_dimid;
    if ((retval = nc_def_var(ncf->ncid_Trj, CANGU_NAME, NC_CHAR, 2, dimid2, &ncf->varid_cellAngular_Trj)))
        ERR(retval);
    //  float64 cell_lengths(frame,cell_spatial)
    dimid2[0] = frame_dimid;
    dimid2[1] = cSpat_dimid;
    if ((retval = nc_def_var(ncf->ncid_Trj, CLENS_NAME, NC_DOUBLE, 2, dimid2, &ncf->varid_cellLengths_Trj)))       // ???
        ERR(retval);
    //  float64 cell_angles(frame,cell_angular)
    dimid2[0] = frame_dimid;
    dimid2[1] = cAngu_dimid;
    if ((retval = nc_def_var(ncf->ncid_Trj, CANGLS_NAME, NC_DOUBLE, 2, dimid2, &ncf->varid_cellAngles_Trj)))       // ???
        ERR(retval);
    
    /* Write Global attributes */
    char title[13];
    strcpy(title,"default_name");
    if ((retval = nc_put_att_text(ncf->ncid_Trj, NC_GLOBAL, "title", strlen(title), title)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncf->ncid_Trj, NC_GLOBAL, "application", 5, "AMBER")))
        ERR(retval);
    if ((retval = nc_put_att_text(ncf->ncid_Trj, NC_GLOBAL, "program", 9, "damselfly")))
        ERR(retval);
    if ((retval = nc_put_att_text(ncf->ncid_Trj, NC_GLOBAL, "programVersion", 3, "1.0")))
        ERR(retval);
    if ((retval = nc_put_att_text(ncf->ncid_Trj, NC_GLOBAL, "Conventions", 5, "AMBER")))    // must be 'AMBER' for trj
        ERR(retval);
    if ((retval = nc_put_att_text(ncf->ncid_Trj, NC_GLOBAL, "ConventionVersion", 3, "1.0")))
        ERR(retval);
    
    /* End define mode. */
    if ((retval = nc_enddef(ncf->ncid_Trj)))
        ERR(retval);

    //------------------------------------------------------
    /* Write spatial, cell_spatial, and cell_angular */
    size_t start, count, start2[2], count2[2];
    char xyz[3];
    char abc[15] = {'a', 'l', 'p', 'h', 'a',
                    'b', 'e', 't', 'a', ' ',
                    'g', 'a', 'm', 'm', 'a' };
    // spatial
    start = 0;
    count = 3;
    xyz[0] = 'x'; xyz[1] = 'y'; xyz[2] = 'z';
    if ((retval = nc_put_vara_text(ncf->ncid_Trj, ncf->varid_spatial_Trj, &start, &count, xyz)))
        ERR(retval);
    // cell_spatial
    xyz[0] = 'a'; xyz[1] = 'b'; xyz[2] = 'c';
    if ((retval = nc_put_vara_text(ncf->ncid_Trj, ncf->varid_cellSpatial_Trj, &start, &count, xyz)))
        ERR(retval);
    // cell_angular
    start2[0] = 0;
    count2[0] = 3;
    start2[1] = 0;
    count2[1] = 5;
    if ((retval = nc_put_vara_text(ncf->ncid_Trj, ncf->varid_cellAngular_Trj, start2, count2, abc)))
        ERR(retval);
}

//----------- W ---------- R e s t a r t --------------------------------------
void writeRst(int nframe, state *sdt, ncfile *ncf){
    int i, retval, natom3;
    double angles[3];
    size_t start, count;
    size_t start2[2], count2[2];

    // convert velocity to half-unit
    natom3 = 3 * sdt->natom;
    /*
    for (i = 0; i < natom3; i++){
        sdt->vel[i] /= 20.4548282808730;
    }
    */
    // write float64 coordinates(atom,spatial) and float64 velocities(atom,spatial)
    start2[0] = 0;
    count2[0] = sdt->natom;
    start2[1] = 0;
    count2[1] = 3;
    if ((retval = nc_put_vara_double(ncf->ncid_Rst, ncf->varid_crd_Rst, start2,
                count2, &sdt->crd[0]))) ERR(retval);
    if ((retval = nc_put_vara_double(ncf->ncid_Rst, ncf->varid_vel_Rst, start2,
                count2, &sdt->vel[0]))) ERR(retval);

    // write float64 cell_lengths(cell_spatial) and float64 cell_angles(cell_angular)
    start = 0;
    count = 3;
    if ((retval = nc_put_vara_double(ncf->ncid_Rst, ncf->varid_cellLengths_Rst, &start,
                &count, &sdt->box[0]))) ERR(retval);
    // Convert from radian to degree
    for (i = 0; i < 3; i++){
        angles[i] = sdt->box[i+3] * 180 / PI;
    }
    if ((retval = nc_put_vara_double(ncf->ncid_Rst, ncf->varid_cellAngles_Rst, &start,
                &count, &angles[0]))) ERR(retval);
}

//----------- W ----------- T r a j e c t o r y--------------------------------
void writeTrj(int currFrame, state *sdt, ncfile *ncf, float *crdF){
    int i, retval, natom3;
    double angles[3];
    size_t start2[2], count2[2];
    size_t start3[3], count3[3];

    // Convert crd from double to single precision for storage
    natom3 = 3 * sdt->natom;
    for (i = 0; i < natom3; i++){
        crdF[i] = (float)(sdt->crd[i]);
    }

    // write float32 coordinates(frame,atom,spatial)
    start3[0] = currFrame;
    count3[0] = 1;   // just write this one frame         // '1' or 'currFrame+1' ?
    start3[1] = 0;
    count3[1] = sdt->natom;
    start3[2] = 0;
    count3[2] = 3;
    if ((retval = nc_put_vara_float(ncf->ncid_Trj, ncf->varid_crd_Trj, start3,
                count3, &crdF[0]))) ERR(retval);

    // write float64 cell_lengths(frame,cell_spatial) and float64 cell_angles(frame,cell_angular)
    start2[0] = currFrame;
    count2[0] = 1;
    start2[1] = 0;
    count2[1] = 3;
    if ((retval = nc_put_vara_double(ncf->ncid_Trj, ncf->varid_cellLengths_Trj, start2,
                count2, &sdt->box[0]))) ERR(retval);
    // Convert from radian to degree
    for (i = 0; i < 3; i++){
        angles[i] = sdt->box[i+3] * 180 / PI;
    }
    if ((retval = nc_put_vara_double(ncf->ncid_Trj, ncf->varid_cellAngles_Trj, start2,
                count2, &angles[0]))) ERR(retval);
}

//----------- C l o s e -- R e s t a r t / T r a j e c t o r y ----------------
void closeRstTrj(ncfile *ncf)
{
   /* Close the file. */
    int retval;
    printf("about to close NC files...\n");
    if ((retval = nc_close(ncf->ncid_Rst)))
        ERR(retval);
    printf("Rst file is correctly closed!\n");
    if ((retval = nc_close(ncf->ncid_Trj)))
        ERR(retval);
    printf("Trj file is correctly closed!\n");
}
