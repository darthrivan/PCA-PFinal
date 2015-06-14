/*
This file is part of ftdock, a program for rigid-body protein-protein docking 
Copyright (C) 1997-2000 Gidon Moont

Biomolecular Modelling Laboratory
Imperial Cancer Research Fund
44 Lincoln's Inn Fields
London WC2A 3PX

+44 (0)20 7269 3348
http://www.bmm.icnet.uk/

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/

#include "structures.h"

#include <sys/times.h>
#include <unistd.h>
#define BUFFER_SIZE (1<<20)


void print_electric_grid ( fftw_real *grid, int grid_size)
{
  int i, diff;

  for (i=0; i<(grid_size * grid_size * ( 2 * ( grid_size / 2 + 1 ) ) * sizeof( fftw_real ))-3; i+=4) 
       write (1, (void *)((char *)grid+i), sizeof(int));  
  if ((diff=i-(grid_size * grid_size * ( 2 * ( grid_size / 2 + 1 ) ) * sizeof( fftw_real )))>0)
       write(1,(void *)((char *)grid+i),diff);

}

int main( int argc , char *argv[] ) {

  /* index counters */

  int i ;

  /* Command line options */

  char    *output_file_name ;
  char    *static_file_name ;
  char    *mobile_file_name ;
  int   global_grid_size ;
  int   angle_step ;
  float   surface ;
  float   internal_value ;
  int   electrostatics ;
  int   keep_per_rotation ;
  int     kept_scores ;
  int   rescue ;
  int   calculate ;
  float   reverse_calculated_one_span ;

  char    *default_global_grid_size ;
  char    *default_angle_step ;
  char    *default_surface ;
  char    *default_internal_value ;
  char    *default_electrostatics ;
  char    *default_keep_per_rotation ;

  /* File stuff */

  FILE    *ftdock_file ;
  char      out_buffer[BUFFER_SIZE];
  char    line_buffer[100] ;
  int   id , id2 , SCscore ;
  float   RPscore ;
  int   x , y , z , z_twist , theta , phi ;

  /* Angles stuff */

  struct Angle  Angles ;
  int   first_rotation , rotation ;

  /* Structures */

  struct Structure  Static_Structure , Mobile_Structure ;
  struct Structure  Origin_Static_Structure , Origin_Mobile_Structure ;
  struct Structure  Rotated_at_Origin_Mobile_Structure ;

  /* Co-ordinates */

  int   xyz , fx , fy , fz , fxyz ;

  /* Grid stuff */

  float   grid_span , one_span ;

  fftw_real *static_grid ;
  fftw_real *mobile_grid ;
  fftw_real *convoluted_grid ;

  fftw_real *static_elec_grid = ( void * ) 0 ;
  fftw_real *mobile_elec_grid = ( void * ) 0 ;
  fftw_real *convoluted_elec_grid = ( void * ) 0 ;

  /* FFTW stuff */

  //rfftwnd_plan  p , pinv ;
  fftwf_plan  plan_stat, plan_stat_el, plan_mob, plan_mob_el;
  fftwf_plan pinv_multi, pinv_multi_el;

  fftwf_complex  *static_fsg ;
  fftwf_complex  *mobile_fsg ;
  fftwf_complex  *multiple_fsg ;

  fftwf_complex  *static_elec_fsg = ( void * ) 0 ;
  fftwf_complex  *mobile_elec_fsg = ( void * ) 0 ;
  fftwf_complex  *multiple_elec_fsg = ( void * ) 0 ;

  /* Scores */

  struct Score  *Scores ;
  float   max_es_value ;

  /* Timing */
  struct tms start, end;
  
/************/

  /* Its nice to tell people what going on straight away */

  setvbuf( stdout , out_buffer, _IOFBF , BUFFER_SIZE ) ;


  printf( "\n          3D-Dock Suite (March 2001)\n" ) ;
  printf( "          Copyright (C) 1997-2000 Gidon Moont\n" ) ;
  printf( "          This program comes with ABSOLUTELY NO WARRANTY\n" ) ;
  printf( "          for details see license. This program is free software,\n"); 
  printf( "          and you may redistribute it under certain conditions.\n\n"); 

  printf( "          Biomolecular Modelling Laboratory\n" ) ;
  printf( "          Imperial Cancer Research Fund\n" ) ;
  printf( "          44 Lincoln's Inn Fields\n" ) ;
  printf( "          London WC2A 3PX\n" ) ;
  printf( "          +44 (0)20 7269 3348\n" ) ;
  printf( "          http://www.bmm.icnet.uk/\n\n" ) ;


  printf( "Starting FTDock (v2.0) global search program\n" ) ;


/************/

  /* Memory allocation */

  if( ( ( output_file_name  = ( char * ) malloc ( 500 * sizeof( char ) ) ) == NULL ) ||
      ( ( static_file_name  = ( char * ) malloc ( 500 * sizeof( char ) ) ) == NULL ) ||
      ( ( mobile_file_name  = ( char * ) malloc ( 500 * sizeof( char ) ) ) == NULL ) ) {
    GENERAL_MEMORY_PROBLEM 
  }

/************/

  /* Command Line defaults */

  strcpy( output_file_name , "ftdock_global.dat" ) ;
  strcpy( static_file_name , " --static file name was not provided--" ) ;
  strcpy( mobile_file_name , " --mobile file name was not provided--" ) ;
  global_grid_size = 128 ;
  angle_step = 12 ;
  surface = 1.3 ;
  internal_value = -15 ;
  electrostatics = 1 ;
  keep_per_rotation = 3 ;
  rescue = 0 ;
  calculate = 1 ;
  reverse_calculated_one_span = 0.7 ;

  default_global_grid_size = "(default calculated)" ;
  default_angle_step = "(default)" ;
  default_surface = "(default)" ;
  default_internal_value = "(default)" ;
  default_electrostatics = "(default)" ;
  default_keep_per_rotation = "(default)" ;

  /* Command Line parse */

  for( i = 1 ; i < argc ; i ++ ) {

    if( strcmp( argv[i] , "-out" ) == 0 ) {
      i ++ ;
      if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
        printf( "Bad command line\n" ) ;
        exit( EXIT_FAILURE ) ;
      }
      strcpy( output_file_name , argv[i] ) ;
    } else {
      if( strcmp( argv[i] , "-static" ) == 0 ) {
        i ++ ;
        if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
          printf( "Bad command line\n" ) ;
          exit( EXIT_FAILURE ) ;
        }
        strcpy( static_file_name , argv[i] ) ;
      } else {
        if( strcmp( argv[i] , "-mobile" ) == 0 ) {
          i ++ ;
          if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
            printf( "Bad command line\n" ) ;
            exit( EXIT_FAILURE ) ;
          }
          strcpy( mobile_file_name , argv[i] ) ;
        } else {
          if( strcmp( argv[i] , "-grid" ) == 0 ) {
            i ++ ;
            if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
              printf( "Bad command line\n" ) ;
              exit( EXIT_FAILURE ) ;
            }
            sscanf( argv[i] , "%d" , &global_grid_size ) ;
            if( ( global_grid_size % 2 ) != 0 ) {
              printf( "Grid size must be even\n" ) ;
              exit( EXIT_FAILURE ) ;
            }
            default_global_grid_size = "(user defined)" ;
            calculate = 0 ;
          } else {
            if( strcmp( argv[i] , "-angle_step" ) == 0 ) {
              i ++ ;
              if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
                printf( "Bad command line\n" ) ;
                exit( EXIT_FAILURE ) ;
              }
              sscanf( argv[i] , "%d" , &angle_step ) ;
              default_angle_step = "(user defined)" ;
            } else {
              if( strcmp( argv[i] , "-surface" ) == 0 ) {
                i ++ ;
                if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
                  printf( "Bad command line\n" ) ;
                  exit( EXIT_FAILURE ) ;
                }
                sscanf( argv[i] , "%f" , &surface ) ;
                default_surface = "(user defined)" ;
              } else {
                if( strcmp( argv[i] , "-internal" ) == 0 ) {
                  i ++ ;
                  if( i == argc ) {
                    printf( "Bad command line\n" ) ;
                    exit( EXIT_FAILURE ) ;
                  }
                  sscanf( argv[i] , "%f" , &internal_value ) ;
                  default_internal_value = "(user defined)" ;
                } else {
                  if( strcmp( argv[i] , "-noelec" ) == 0 ) {
                    electrostatics = 0 ;
                    default_electrostatics = "(user defined)" ;
                  } else {
                    if( strcmp( argv[i] , "-keep" ) == 0 ) {
                      i ++ ;
                      if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
                        printf( "Bad command line\n" ) ;
                        exit( EXIT_FAILURE ) ;
                      }
                      sscanf( argv[i] , "%d" , &keep_per_rotation ) ;
                      default_keep_per_rotation = "(user defined)" ;
                    } else {
                      if( strcmp( argv[i] , "-rescue" ) == 0 ) {
                        rescue = 1 ;
                      } else {
                        if( strcmp( argv[i] , "-calculate_grid" ) == 0 ) {
                          i ++ ;
                          if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
                            printf( "Bad command line\n" ) ;
                            exit( EXIT_FAILURE ) ;
                          }
                          calculate = 1 ;
                          default_global_grid_size = "(user defined calculated)" ;
                          sscanf( argv[i] , "%f" , &reverse_calculated_one_span ) ;
                        } else {
                          printf( "Bad command line\n" ) ;
                          exit( EXIT_FAILURE ) ;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

  }

/************/

  /* Rescue option */

  if( rescue == 1 ) {

    printf( "RESCUE mode\n" ) ;

    if( ( ftdock_file = fopen( "scratch_parameters.dat" , "r" ) ) == NULL ) {
      printf( "Could not open scratch_parameters.dat for reading.\nDying\n\n" ) ;
      exit( EXIT_FAILURE ) ;
    }

    calculate = 0 ;

    default_global_grid_size = "(read from rescue file)" ;
    default_angle_step = "(read from rescue file)" ;
    default_surface = "(read from rescue file)" ;
    default_internal_value = "(read from rescue file)" ;
    default_electrostatics = "(read from rescue file)" ;

    while( fgets( line_buffer , 99 , ftdock_file ) ) {

      if( strncmp( line_buffer , "Static molecule" , 15 ) == 0 ) sscanf( line_buffer , "Static molecule :: %s" , static_file_name ) ;
      if( strncmp( line_buffer , "Mobile molecule" , 15 ) == 0 ) sscanf( line_buffer , "Mobile molecule :: %s" , mobile_file_name ) ;
      if( strncmp( line_buffer , "Output file name" , 16 ) == 0 ) sscanf( line_buffer , "Output file name :: %s" , output_file_name ) ;
      if( strncmp( line_buffer , "Global grid size" , 16 ) == 0 ) sscanf( line_buffer , "Global grid size :: %d" , &global_grid_size ) ;
      if( strncmp( line_buffer , "Global search angle step" , 24 ) == 0 ) sscanf( line_buffer , "Global search angle step :: %d" , &angle_step ) ;
      if( strncmp( line_buffer , "Global surface thickness" , 24 ) == 0 ) sscanf( line_buffer , "Global surface thickness :: %f" , &surface ) ;
      if( strncmp( line_buffer , "Global internal deterrent value" , 31 ) == 0 ) sscanf( line_buffer , "Global internal deterrent value :: %f" , &internal_value ) ;
      if( strncmp( line_buffer , "Electrostatics                     ::     on" , 44 ) == 0 ) electrostatics = 1 ;    
      if( strncmp( line_buffer , "Electrostatics                     ::    off" , 44 ) == 0 ) electrostatics = 0 ;    
      if( strncmp( line_buffer , "Global keep per rotation" , 25 ) == 0 ) sscanf( line_buffer , "Global keep per rotation :: %d" , &keep_per_rotation ) ;

    }

    fclose( ftdock_file ) ;

    if( ( ftdock_file = fopen( "scratch_scores.dat" , "r" ) ) == NULL ) {
      printf( "Could not open scratch_scores.dat for reading.\nDying\n\n" ) ;
      exit( EXIT_FAILURE ) ;
    }

    fgets( line_buffer , 99 , ftdock_file ) ;

    while( fgets( line_buffer , 99 , ftdock_file ) ) {

      sscanf( line_buffer , "G_DATA %d " , &first_rotation ) ;

    }

    fclose( ftdock_file ) ;

    first_rotation ++ ;

    printf( "Will be starting from rotation %d\n" , first_rotation ) ;

/************/

  } else {

    first_rotation = 1 ;

  }

/************/

  /* Do these things first so that bad inputs will be caught soonest */

  /* Read in Structures from pdb files */
  Static_Structure = read_pdb_to_structure( static_file_name ) ;
  Mobile_Structure = read_pdb_to_structure( mobile_file_name ) ;

  if( Mobile_Structure.length > Static_Structure.length ) {
    printf( "WARNING\n" ) ;
    printf( "The mobile molecule has more residues than the static\n" ) ;
    printf( "Are you sure you have the correct molecules?\n" ) ;
    printf( "Continuing anyway\n" ) ;
  }
  
/************/

  /* Get angles */
  Angles = generate_global_angles( angle_step ) ;

  printf( "Total number of rotations is %d\n" , Angles.n ) ;

/************/

  /* Assign charges */

  if( electrostatics == 1 ) {
    printf( "Assigning charges\n" ) ;
    assign_charges( Static_Structure ) ;
    assign_charges( Mobile_Structure ) ;

  /************/

    /* Store new structures centered on Origin */

    Origin_Static_Structure = translate_structure_onto_origin( Static_Structure ) ;
    Origin_Mobile_Structure = translate_structure_onto_origin( Mobile_Structure ) ;

    /* Free some memory */

    for( i = 1 ; i <= Static_Structure.length ; i ++ ) {
      free( Static_Structure.Residue[i].Atom ) ;
    }
    free( Static_Structure.Residue ) ;

    for( i = 1 ; i <= Mobile_Structure.length ; i ++ ) {
      free( Mobile_Structure.Residue[i].Atom ) ;
    }
    free( Mobile_Structure.Residue ) ;

  /************/

    /* Calculate Grid stuff */

    grid_span = total_span_of_structures( Origin_Static_Structure , Origin_Mobile_Structure ) ;

    if( calculate == 1 ) {
      printf( "Using automatic calculation for grid size\n" ) ;
      global_grid_size = (int)( grid_span / reverse_calculated_one_span ) ;
      if( ( global_grid_size % 2 ) != 0 ) global_grid_size ++ ;
    }

    one_span = grid_span / (float)global_grid_size ;

    printf( "Span = %.3f angstroms\n" , grid_span ) ;
    printf( "Grid size = %d\n" , global_grid_size ) ;
    printf( "Each Grid cube = %.5f angstroms\n" , one_span ) ;

  /************/

    /* Memory Allocation */

    if( ( Scores = ( struct Score * ) fftwf_malloc ( ( keep_per_rotation + 2 ) * sizeof( struct Score ) ) ) == NULL ) {
      GENERAL_MEMORY_PROBLEM
    }

    if(
      ( ( static_grid = ( fftw_real * ) fftwf_malloc
       ( global_grid_size * global_grid_size * ( 2 * ( global_grid_size / 2 + 1 ) ) * sizeof( fftw_real ) ) ) == NULL )
      ||
      ( ( mobile_grid = ( fftw_real * ) fftwf_malloc
       ( global_grid_size * global_grid_size * ( 2 * ( global_grid_size / 2 + 1 ) ) * sizeof( fftw_real ) ) ) == NULL )
      ||
      ( ( convoluted_grid = ( fftw_real * ) fftwf_malloc
       ( global_grid_size * global_grid_size * ( 2 * ( global_grid_size / 2 + 1 ) ) * sizeof( fftw_real ) ) ) == NULL )
      ) {
      printf( "Not enough memory for surface grids\nUse (sensible) smaller grid size\nDying\n\n" ) ;
      exit( EXIT_FAILURE ) ;
    }

    static_fsg = ( fftwf_complex * ) static_grid ;
    mobile_fsg = ( fftwf_complex * ) mobile_grid ;
    multiple_fsg = ( fftwf_complex * ) convoluted_grid ;

    if(
      ( ( static_elec_grid = ( fftw_real * ) fftwf_malloc
       ( global_grid_size * global_grid_size * ( 2 * ( global_grid_size / 2 + 1 ) ) * sizeof( fftw_real ) ) ) == NULL )
      ||
      ( ( mobile_elec_grid = ( fftw_real * ) fftwf_malloc
       ( global_grid_size * global_grid_size * ( 2 * ( global_grid_size / 2 + 1 ) ) * sizeof( fftw_real ) ) ) == NULL )
      ||
      ( ( convoluted_elec_grid = ( fftw_real * ) fftwf_malloc
       ( global_grid_size * global_grid_size * ( 2 * ( global_grid_size / 2 + 1 ) ) * sizeof( fftw_real ) ) ) == NULL )
      ) {
      printf( "Not enough memory for electrostatic grids\nSwitch off electrostatics or use (sensible) smaller grid size\nDying\n\n" ) ;
      exit( EXIT_FAILURE ) ;
    } else {
       // all ok 
      printf( "Electrostatics are on\n" ) ;
    }

    static_elec_fsg = ( fftwf_complex * ) static_elec_grid ;
    mobile_elec_fsg = ( fftwf_complex * ) mobile_elec_grid ;
    multiple_elec_fsg = ( fftwf_complex * ) convoluted_elec_grid ;

  /************/

     // Create FFTW plans 
    printf( "Creating plans\n" ) ;
    plan_stat       = fftwf_plan_dft_r2c_3d( global_grid_size , global_grid_size , global_grid_size ,
                                 static_grid , (fftwf_complex *)static_grid, FFTW_MEASURE ) ;
    plan_stat_el    = fftwf_plan_dft_r2c_3d( global_grid_size , global_grid_size , global_grid_size ,
                                 static_elec_grid , (fftwf_complex *)static_elec_grid, FFTW_MEASURE ) ;
                     
    plan_mob       = fftwf_plan_dft_r2c_3d( global_grid_size , global_grid_size , global_grid_size ,
                                 mobile_grid , (fftwf_complex *)mobile_grid, FFTW_MEASURE ) ;
    plan_mob_el    = fftwf_plan_dft_r2c_3d( global_grid_size , global_grid_size , global_grid_size ,
                                 mobile_elec_grid , (fftwf_complex *)mobile_elec_grid, FFTW_MEASURE ) ;
                     
    pinv_multi    = fftwf_plan_dft_c2r_3d( global_grid_size , global_grid_size , global_grid_size ,
                                 multiple_fsg , (fftw_real *)multiple_fsg, FFTW_MEASURE ) ;
    pinv_multi_el = fftwf_plan_dft_c2r_3d( global_grid_size , global_grid_size , global_grid_size ,
                                 multiple_elec_fsg , (fftw_real *)multiple_elec_fsg, FFTW_MEASURE ) ;
  /************/

    //printf( "PCA TIMING SHOULD start here\n");
    
      //TIMING
      //struct tms start, end;

    if (times(&start) == (clock_t)-1) exit(0);
      
    printf( "Setting up Static Structure\n" ) ;

    /* Discretise and surface the Static Structure (need do only once) */
    discretise_structure( Origin_Static_Structure , grid_span , global_grid_size , static_grid ) ;
    printf( "  surfacing grid\n" ) ;
    surface_grid( grid_span , global_grid_size , static_grid , surface , internal_value ) ;

    /* Calculate electic field at all grid nodes (need do only once) */
    electric_field( Origin_Static_Structure , grid_span , global_grid_size , static_elec_grid ) ;
    electric_field_zero_core( global_grid_size , static_elec_grid , static_grid , internal_value ) ;

    /* Fourier Transform the static grids (need do only once) */
    printf( "  one time forward FFT calculations\n" ) ;
    fftwf_execute(plan_stat);
    fftwf_execute(plan_stat_el);

    printf( "  done\n" ) ;

  /************/

    /* Store paramaters in case of rescue */

    if( ( ftdock_file = fopen( "scratch_parameters.dat" , "w" ) ) == NULL ) {
      printf( "Could not open scratch_parameters.dat for writing.\nDying\n\n" ) ;
      exit( EXIT_FAILURE ) ;
    }

    fprintf( ftdock_file, "\nGlobal Scan\n" ) ;

    fprintf( ftdock_file, "\nCommand line controllable values\n" ) ;
    fprintf( ftdock_file, "Static molecule                    :: %s\n" , static_file_name ) ;
    fprintf( ftdock_file, "Mobile molecule                    :: %s\n" , mobile_file_name ) ;
    fprintf( ftdock_file, "Output file name                   :: %s\n" , output_file_name ) ;
    fprintf( ftdock_file, "\n" ) ;
    fprintf( ftdock_file, "Global grid size                   :: %6d      %s\n" , global_grid_size , default_global_grid_size ) ;
    fprintf( ftdock_file, "Global search angle step           :: %6d      %s\n" , angle_step , default_angle_step ) ;
    fprintf( ftdock_file, "Global surface thickness           :: %9.2f   %s\n" , surface , default_surface ) ;
    fprintf( ftdock_file, "Global internal deterrent value    :: %9.2f   %s\n" , internal_value , default_internal_value ) ;
    fprintf( ftdock_file, "Electrostatics                     ::     on      %s\n" , default_electrostatics ) ;
    fprintf( ftdock_file, "Global keep per rotation           :: %6d      %s\n" , keep_per_rotation , default_keep_per_rotation ) ;

    fprintf( ftdock_file, "\nCalculated values\n" ) ;
    fprintf( ftdock_file, "Global rotations                   :: %6d\n" , Angles.n ) ;
    fprintf( ftdock_file, "Global total span (angstroms)      :: %10.3f\n" , grid_span ) ;
    fprintf( ftdock_file, "Global grid cell span (angstroms)  :: %10.3f\n" , one_span ) ;

    fclose( ftdock_file ) ;

  /************/

    /* Main program loop */

    max_es_value = 0 ;

    printf( "Starting main loop through the rotations\n" ) ;

    for( rotation = first_rotation ; rotation <= 20/*Angles.n*/ ; rotation ++ ) {
      /* Rotate Mobile Structure */
      Rotated_at_Origin_Mobile_Structure =
       rotate_structure( Origin_Mobile_Structure , (int)Angles.z_twist[rotation] , (int)Angles.theta[rotation] , (int)Angles.phi[rotation] ) ;

      /* Discretise the rotated Mobile Structure */
      discretise_structure( Rotated_at_Origin_Mobile_Structure , grid_span , global_grid_size , mobile_grid ) ;

      /* Electic point charge approximation onto grid calculations ( quicker than filed calculations by a long way! ) */
      electric_point_charge( Rotated_at_Origin_Mobile_Structure , grid_span , global_grid_size , mobile_elec_grid ) ;

      /* Forward Fourier Transforms */
      fftwf_execute(plan_mob);
      fftwf_execute(plan_mob_el);

  /************/

      /* Do convolution of the two sets of grids
         convolution is equivalent to multiplication of the complex conjugate of one
         fourier grid with other (raw) one
         hence the sign changes from a normal complex number multiplication
      */

      for( fx = 0 ; fx < global_grid_size ; fx ++ ) {
        for( fy = 0 ; fy < global_grid_size ; fy ++ ) {
          for( fz = 0 ; fz < global_grid_size/2 + 1 ; fz ++ ) {

            fxyz = fz + ( global_grid_size/2 + 1 ) * ( fy + global_grid_size * fx ) ;
            multiple_fsg[fxyz][0] =
             static_fsg[fxyz][0] * mobile_fsg[fxyz][0] + static_fsg[fxyz][1] * mobile_fsg[fxyz][1] ;
            multiple_fsg[fxyz][1] =
             static_fsg[fxyz][1] * mobile_fsg[fxyz][0] - static_fsg[fxyz][0] * mobile_fsg[fxyz][1] ;           
            multiple_elec_fsg[fxyz][0] =
             static_elec_fsg[fxyz][0] * mobile_elec_fsg[fxyz][0] + static_elec_fsg[fxyz][1] * mobile_elec_fsg[fxyz][1] ;
            multiple_elec_fsg[fxyz][1] =
             static_elec_fsg[fxyz][1] * mobile_elec_fsg[fxyz][0] - static_elec_fsg[fxyz][0] * mobile_elec_fsg[fxyz][1] ;
          }
        }
      }

      /* Reverse Fourier Transform */
      fftwf_execute(pinv_multi);
      fftwf_execute(pinv_multi_el);

  /************/

      /* Get best scores */

      for( i = 0 ; i < keep_per_rotation ; i ++ ) {

        Scores[i].score = 0 ;
        Scores[i].rpscore = 0.0 ;
        Scores[i].coord[1] = 0 ;
        Scores[i].coord[2] = 0 ;
        Scores[i].coord[3] = 0 ;

      }

      for( x = 0 ; x < global_grid_size ; x ++ ) {
        fx = x ;
        if( fx > ( global_grid_size / 2 ) ) fx -= global_grid_size ;

        for( y = 0 ; y < global_grid_size ; y ++ ) {
          fy = y ;
          if( fy > ( global_grid_size / 2 ) ) fy -= global_grid_size ;

          for( z = 0 ; z < global_grid_size ; z ++ ) {
            fz = z ;
            if( fz > ( global_grid_size / 2 ) ) fz -= global_grid_size ;

            xyz = z + ( 2 * ( global_grid_size / 2 + 1 ) ) * ( y + global_grid_size * x ) ;

            if ( convoluted_elec_grid[xyz] < 0 ) {

              /* Scale factor from FFTs */
              convoluted_grid[xyz] /= ( global_grid_size * global_grid_size * global_grid_size ) ;

              if( (int)convoluted_grid[xyz] > Scores[keep_per_rotation-1].score ) {

                i = keep_per_rotation - 2 ;

                while( ( (int)convoluted_grid[xyz] > Scores[i].score ) && ( i >= 0 ) ) {
                  Scores[i+1].score    = Scores[i].score ;
                  Scores[i+1].rpscore  = Scores[i].rpscore ;
                  Scores[i+1].coord[1] = Scores[i].coord[1] ;
                  Scores[i+1].coord[2] = Scores[i].coord[2] ;
                  Scores[i+1].coord[3] = Scores[i].coord[3] ;
                  i -- ;
                }

                Scores[i+1].score    = (int)convoluted_grid[xyz] ;
                Scores[i+1].rpscore  = (float)convoluted_elec_grid[xyz] ;
                Scores[i+1].coord[1] = fx ;
                Scores[i+1].coord[2] = fy ;
                Scores[i+1].coord[3] = fz ;

              }

            }

          }
        }
      }

      for( i = 0 ; i < keep_per_rotation ; i ++ ) {

        max_es_value = min( max_es_value , Scores[i].rpscore ) ;

        fprintf( stdout, "G_DATA %6d   %6d    %7d       %4d %4d %4d      %4d%4d%4d\n" ,
                  rotation , 0 , Scores[i].score , Scores[i].coord[1] , Scores[i].coord[2] , Scores[i].coord[3 ] ,
                   Angles.z_twist[rotation] , Angles.theta[rotation]  , Angles.phi[rotation] ) ;

      }


      /* Free some memory */
      for( i = 1 ; i <= Rotated_at_Origin_Mobile_Structure.length ; i ++ ) {
        free( Rotated_at_Origin_Mobile_Structure.Residue[i].Atom ) ;
      }
      free( Rotated_at_Origin_Mobile_Structure.Residue ) ;
    }

    /* Finished main loop */
  }
  else {

  /************/

    /* Store new structures centered on Origin */

    Origin_Static_Structure = translate_structure_onto_origin( Static_Structure ) ;
    Origin_Mobile_Structure = translate_structure_onto_origin( Mobile_Structure ) ;

    /* Free some memory */

    for( i = 1 ; i <= Static_Structure.length ; i ++ ) {
      free( Static_Structure.Residue[i].Atom ) ;
    }
    free( Static_Structure.Residue ) ;

    for( i = 1 ; i <= Mobile_Structure.length ; i ++ ) {
      free( Mobile_Structure.Residue[i].Atom ) ;
    }
    free( Mobile_Structure.Residue ) ;

  /************/

    /* Calculate Grid stuff */

    grid_span = total_span_of_structures( Origin_Static_Structure , Origin_Mobile_Structure ) ;

    if( calculate == 1 ) {
      printf( "Using automatic calculation for grid size\n" ) ;
      global_grid_size = (int)( grid_span / reverse_calculated_one_span ) ;
      if( ( global_grid_size % 2 ) != 0 ) global_grid_size ++ ;
    }

    one_span = grid_span / (float)global_grid_size ;

    printf( "Span = %.3f angstroms\n" , grid_span ) ;
    printf( "Grid size = %d\n" , global_grid_size ) ;
    printf( "Each Grid cube = %.5f angstroms\n" , one_span ) ;

  /************/

    /* Memory Allocation */

    if( ( Scores = ( struct Score * ) fftwf_malloc ( ( keep_per_rotation + 2 ) * sizeof( struct Score ) ) ) == NULL ) {
      GENERAL_MEMORY_PROBLEM
    }

    if(
      ( ( static_grid = ( fftw_real * ) fftwf_malloc
       ( global_grid_size * global_grid_size * ( 2 * ( global_grid_size / 2 + 1 ) ) * sizeof( fftw_real ) ) ) == NULL )
      ||
      ( ( mobile_grid = ( fftw_real * ) fftwf_malloc
       ( global_grid_size * global_grid_size * ( 2 * ( global_grid_size / 2 + 1 ) ) * sizeof( fftw_real ) ) ) == NULL )
      ||
      ( ( convoluted_grid = ( fftw_real * ) fftwf_malloc
       ( global_grid_size * global_grid_size * ( 2 * ( global_grid_size / 2 + 1 ) ) * sizeof( fftw_real ) ) ) == NULL )
      ) {
      printf( "Not enough memory for surface grids\nUse (sensible) smaller grid size\nDying\n\n" ) ;
      exit( EXIT_FAILURE ) ;
    }

    static_fsg = ( fftwf_complex * ) static_grid ;
    mobile_fsg = ( fftwf_complex * ) mobile_grid ;
    multiple_fsg = ( fftwf_complex * ) convoluted_grid ;

  /************/

     // Create FFTW plans 
    printf( "Creating plans\n" ) ;
    plan_stat       = fftwf_plan_dft_r2c_3d( global_grid_size , global_grid_size , global_grid_size ,
                                 static_grid , (fftwf_complex *)static_grid, FFTW_MEASURE ) ;
    plan_stat_el    = fftwf_plan_dft_r2c_3d( global_grid_size , global_grid_size , global_grid_size ,
                                 static_elec_grid , (fftwf_complex *)static_elec_grid, FFTW_MEASURE ) ;
               
    plan_mob       = fftwf_plan_dft_r2c_3d( global_grid_size , global_grid_size , global_grid_size ,
                                 mobile_grid , (fftwf_complex *)mobile_grid, FFTW_MEASURE ) ;
    plan_mob_el    = fftwf_plan_dft_r2c_3d( global_grid_size , global_grid_size , global_grid_size ,
                                 mobile_elec_grid , (fftwf_complex *)mobile_elec_grid, FFTW_MEASURE ) ;
               
    pinv_multi    = fftwf_plan_dft_c2r_3d( global_grid_size , global_grid_size , global_grid_size ,
                                 multiple_fsg , (fftw_real *)multiple_fsg, FFTW_MEASURE ) ;
    pinv_multi_el = fftwf_plan_dft_c2r_3d( global_grid_size , global_grid_size , global_grid_size ,
                                 multiple_elec_fsg , (fftw_real *)multiple_elec_fsg, FFTW_MEASURE ) ;
  /************/

    //printf( "PCA TIMING SHOULD start here\n");
    
    //TIMING
    //struct tms start, end;

      if (times(&start) == (clock_t)-1) exit(0);
      
    printf( "Setting up Static Structure\n" ) ;

    /* Discretise and surface the Static Structure (need do only once) */
    discretise_structure( Origin_Static_Structure , grid_span , global_grid_size , static_grid ) ;
    printf( "  surfacing grid\n" ) ;
    surface_grid( grid_span , global_grid_size , static_grid , surface , internal_value ) ;

    /* Fourier Transform the static grids (need do only once) */
    printf( "  one time forward FFT calculations\n" ) ;
    fftwf_execute(plan_stat);
    printf( "  done\n" ) ;

  /************/

    /* Store paramaters in case of rescue */

    if( ( ftdock_file = fopen( "scratch_parameters.dat" , "w" ) ) == NULL ) {
      printf( "Could not open scratch_parameters.dat for writing.\nDying\n\n" ) ;
      exit( EXIT_FAILURE ) ;
    }

    fprintf( ftdock_file, "\nGlobal Scan\n" ) ;

    fprintf( ftdock_file, "\nCommand line controllable values\n" ) ;
    fprintf( ftdock_file, "Static molecule                    :: %s\n" , static_file_name ) ;
    fprintf( ftdock_file, "Mobile molecule                    :: %s\n" , mobile_file_name ) ;
    fprintf( ftdock_file, "Output file name                   :: %s\n" , output_file_name ) ;
    fprintf( ftdock_file, "\n" ) ;
    fprintf( ftdock_file, "Global grid size                   :: %6d      %s\n" , global_grid_size , default_global_grid_size ) ;
    fprintf( ftdock_file, "Global search angle step           :: %6d      %s\n" , angle_step , default_angle_step ) ;
    fprintf( ftdock_file, "Global surface thickness           :: %9.2f   %s\n" , surface , default_surface ) ;
    fprintf( ftdock_file, "Global internal deterrent value    :: %9.2f   %s\n" , internal_value , default_internal_value ) ;
    fprintf( ftdock_file, "Electrostatics                     ::    off      %s\n" , default_electrostatics ) ;
    fprintf( ftdock_file, "Global keep per rotation           :: %6d      %s\n" , keep_per_rotation , default_keep_per_rotation ) ;

    fprintf( ftdock_file, "\nCalculated values\n" ) ;
    fprintf( ftdock_file, "Global rotations                   :: %6d\n" , Angles.n ) ;
    fprintf( ftdock_file, "Global total span (angstroms)      :: %10.3f\n" , grid_span ) ;
    fprintf( ftdock_file, "Global grid cell span (angstroms)  :: %10.3f\n" , one_span ) ;

    fclose( ftdock_file ) ;

  /************/

    /* Main program loop */

    max_es_value = 0 ;

    printf( "Starting main loop through the rotations\n" ) ;
    for( rotation = first_rotation ; rotation <= 20/*Angles.n*/ ; rotation ++ ) {
      /* Rotate Mobile Structure */
      Rotated_at_Origin_Mobile_Structure =
       rotate_structure( Origin_Mobile_Structure , (int)Angles.z_twist[rotation] , (int)Angles.theta[rotation] , (int)Angles.phi[rotation] ) ;

      /* Discretise the rotated Mobile Structure */
      discretise_structure( Rotated_at_Origin_Mobile_Structure , grid_span , global_grid_size , mobile_grid ) ;

      /* Electic point charge approximation onto grid calculations ( quicker than filed calculations by a long way! ) */
      /* Forward Fourier Transforms */
      fftwf_execute(plan_mob);

  /************/

      /* Do convolution of the two sets of grids
         convolution is equivalent to multiplication of the complex conjugate of one
         fourier grid with other (raw) one
         hence the sign changes from a normal complex number multiplication
      */

      for( fx = 0 ; fx < global_grid_size ; fx ++ ) {
        for( fy = 0 ; fy < global_grid_size ; fy ++ ) {
          for( fz = 0 ; fz < global_grid_size/2 + 1 ; fz ++ ) {

            fxyz = fz + ( global_grid_size/2 + 1 ) * ( fy + global_grid_size * fx ) ;
            multiple_fsg[fxyz][0] =
             static_fsg[fxyz][0] * mobile_fsg[fxyz][0] + static_fsg[fxyz][1] * mobile_fsg[fxyz][1] ;
            multiple_fsg[fxyz][1] =
             static_fsg[fxyz][1] * mobile_fsg[fxyz][0] - static_fsg[fxyz][0] * mobile_fsg[fxyz][1] ;
          }
        }
      }

      /* Reverse Fourier Transform */
      fftwf_execute(pinv_multi);

  /************/

      /* Get best scores */

      for( i = 0 ; i < keep_per_rotation ; i ++ ) {

        Scores[i].score = 0 ;
        Scores[i].rpscore = 0.0 ;
        Scores[i].coord[1] = 0 ;
        Scores[i].coord[2] = 0 ;
        Scores[i].coord[3] = 0 ;

      }

      for( x = 0 ; x < global_grid_size ; x ++ ) {
        fx = x ;
        if( fx > ( global_grid_size / 2 ) ) fx -= global_grid_size ;

        for( y = 0 ; y < global_grid_size ; y ++ ) {
          fy = y ;
          if( fy > ( global_grid_size / 2 ) ) fy -= global_grid_size ;

          for( z = 0 ; z < global_grid_size ; z ++ ) {
            fz = z ;
            if( fz > ( global_grid_size / 2 ) ) fz -= global_grid_size ;

            xyz = z + ( 2 * ( global_grid_size / 2 + 1 ) ) * ( y + global_grid_size * x ) ;
            /* Scale factor from FFTs */
            if( (int)convoluted_grid[xyz] != 0 ) {
              convoluted_grid[xyz] /= ( global_grid_size * global_grid_size * global_grid_size ) ;
            }

            if( (int)convoluted_grid[xyz] > Scores[keep_per_rotation-1].score ) {

              i = keep_per_rotation - 2 ;

              while( ( (int)convoluted_grid[xyz] > Scores[i].score ) && ( i >= 0 ) ) {
                Scores[i+1].score    = Scores[i].score ;
                Scores[i+1].rpscore  = Scores[i].rpscore ;
                Scores[i+1].coord[1] = Scores[i].coord[1] ;
                Scores[i+1].coord[2] = Scores[i].coord[2] ;
                Scores[i+1].coord[3] = Scores[i].coord[3] ;
                i -- ;
              }

              Scores[i+1].score    = (int)convoluted_grid[xyz] ;
              Scores[i+1].rpscore  = (float)0 ;
              Scores[i+1].coord[1] = fx ;
              Scores[i+1].coord[2] = fy ;
              Scores[i+1].coord[3] = fz ;

            }
          }
        }
      }

      for( i = 0 ; i < keep_per_rotation ; i ++ ) {

        max_es_value = min( max_es_value , Scores[i].rpscore ) ;

        fprintf( stdout, "G_DATA %6d   %6d    %7d       %4d %4d %4d      %4d%4d%4d\n" ,
                  rotation , 0 , Scores[i].score , Scores[i].coord[1] , Scores[i].coord[2] , Scores[i].coord[3 ] ,
                   Angles.z_twist[rotation] , Angles.theta[rotation]  , Angles.phi[rotation] ) ;

      }


      /* Free some memory */
      for( i = 1 ; i <= Rotated_at_Origin_Mobile_Structure.length ; i ++ ) {
        free( Rotated_at_Origin_Mobile_Structure.Residue[i].Atom ) ;
      }
      free( Rotated_at_Origin_Mobile_Structure.Residue ) ;
    }
    // fclose( ftdock_file ) ;

    /* Finished main loop */

  /************/

  }
  /* Free the memory */
  fftwf_destroy_plan( plan_stat ) ;
  fftwf_destroy_plan( plan_stat_el ) ;
  fftwf_destroy_plan( plan_mob ) ;
  fftwf_destroy_plan( plan_mob_el ) ;
  fftwf_destroy_plan( pinv_multi ) ;
  fftwf_destroy_plan( pinv_multi_el ) ;

  fftwf_free( static_grid ) ;
  fftwf_free( mobile_grid ) ;
  fftwf_free( convoluted_grid ) ;

  for( i = 1 ; i <= Origin_Static_Structure.length ; i ++ ) {
    free( Origin_Static_Structure.Residue[i].Atom ) ;
  }
  free( Origin_Static_Structure.Residue ) ;

  for( i = 1 ; i <= Origin_Mobile_Structure.length ; i ++ ) {
    free( Origin_Mobile_Structure.Residue[i].Atom ) ;
  }
  free( Origin_Mobile_Structure.Residue ) ;

      /* PCA: Finishing programm here*/ 
      if (times(&end) == (clock_t)-1) exit(0);
      // fprintf(stderr, "\n Timing amb crida times: user %f segons, system: %f segons\n", (float)(end.tms_utime-start.tms_utime)/sysconf(_SC_CLK_TCK), (float)(end.tms_stime-start.tms_stime)/sysconf(_SC_CLK_TCK));
      fprintf(stderr, "User:%f\nSystem:%f\nElapsed:%f\n", (float)(end.tms_utime-start.tms_utime)/sysconf(_SC_CLK_TCK), (float)(end.tms_stime-start.tms_stime)/sysconf(_SC_CLK_TCK),
        (float)((end.tms_utime-start.tms_utime)+(end.tms_stime-start.tms_stime))/sysconf(_SC_CLK_TCK));

    //  printf("PCA TIMING SHOULD stop here\n");
    //  printf("PCA STOPS HERE\n");
      return 0;
}
