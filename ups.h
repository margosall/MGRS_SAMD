#include "polarst.h"

#define UPS_NO_ERROR                0x0000
#define UPS_LAT_ERROR               0x0001
#define UPS_LON_ERROR               0x0002
#define UPS_HEMISPHERE_ERROR        0x0004
#define UPS_EASTING_ERROR           0x0008
#define UPS_NORTHING_ERROR          0x0010
#define UPS_A_ERROR                 0x0020
#define UPS_INV_F_ERROR             0x0040

#define PI       3.14159265358979323e0  /* PI     */
#define PI_OVER    (PI/2.0e0)           /* PI over 2 */
#define MAX_LAT    ((PI * 90)/180.0)    /* 90 degrees in radians */
#define MAX_ORIGIN_LAT ((81.114528 * PI) / 180.0)
#define MIN_NORTH_LAT (83.5*PI/180.0)
#define MIN_SOUTH_LAT (-79.5*PI/180.0)
#define MIN_EAST_NORTH 0
#define MAX_EAST_NORTH 4000000

/* Ellipsoid Parameters, default to WGS 84  */
static double UPS_a = 6378137.0;          /* Semi-major axis of ellipsoid in meters   */
static double UPS_f = 1 / 298.257223563;  /* Flattening of ellipsoid  */
const double UPS_False_Easting = 2000000;
const double UPS_False_Northing = 2000000;
static double UPS_Origin_Latitude = MAX_ORIGIN_LAT;  /*set default = North Hemisphere */
static double UPS_Origin_Longitude = 0.0;

long Set_UPS_Parameters( double a,
                         double f)
{
/*
 * The function SET_UPS_PARAMETERS receives the ellipsoid parameters and sets
 * the corresponding state variables. If any errors occur, the error code(s)
 * are returned by the function, otherwise UPS_NO_ERROR is returned.
 *
 *   a     : Semi-major axis of ellipsoid in meters (input)
 *   f     : Flattening of ellipsoid                (input)
 */

  double inv_f = 1 / f;
  long Error_Code = UPS_NO_ERROR;

  if (a <= 0.0)
  { /* Semi-major axis must be greater than zero */
    Error_Code |= UPS_A_ERROR;
  }
  if ((inv_f < 250) || (inv_f > 350))
  { /* Inverse flattening must be between 250 and 350 */
    Error_Code |= UPS_INV_F_ERROR;
  }

  if (!Error_Code)
  { /* no errors */
    UPS_a = a;
    UPS_f = f;
  }
  return (Error_Code);
}  /* END of Set_UPS_Parameters  */
long Convert_Geodetic_To_UPS ( double Latitude,
                               double Longitude,
                               char   *Hemisphere,
                               double *Easting,
                               double *Northing)
{
/*
 *  The function Convert_Geodetic_To_UPS converts geodetic (latitude and
 *  longitude) coordinates to UPS (hemisphere, easting, and northing)
 *  coordinates, according to the current ellipsoid parameters. If any 
 *  errors occur, the error code(s) are returned by the function, 
 *  otherwide UPS_NO_ERROR is returned.
 *
 *    Latitude      : Latitude in radians                       (input)
 *    Longitude     : Longitude in radians                      (input)
 *    Hemisphere    : Hemisphere either 'N' or 'S'              (output)
 *    Easting       : Easting/X in meters                       (output)
 *    Northing      : Northing/Y in meters                      (output)
 */

  double tempEasting, tempNorthing;
  long Error_Code = UPS_NO_ERROR;

  if ((Latitude < -MAX_LAT) || (Latitude > MAX_LAT))
  {   /* latitude out of range */
    Error_Code |= UPS_LAT_ERROR;
  }
  if ((Latitude < 0) && (Latitude > MIN_SOUTH_LAT))
    Error_Code |= UPS_LAT_ERROR;
  if ((Latitude >= 0) && (Latitude < MIN_NORTH_LAT))
    Error_Code |= UPS_LAT_ERROR;
  if ((Longitude < -PI) || (Longitude > (2 * PI)))
  {  /* slam out of range */
    Error_Code |= UPS_LON_ERROR;
  }

  if (!Error_Code)
  {  /* no errors */
    if (Latitude < 0)
    {
      UPS_Origin_Latitude = -MAX_ORIGIN_LAT; 
      *Hemisphere = 'S';
    }
    else
    {
      UPS_Origin_Latitude = MAX_ORIGIN_LAT; 
      *Hemisphere = 'N';
    }


    Set_Polar_Stereographic_Parameters( UPS_a,
                                        UPS_f,
                                        UPS_Origin_Latitude,
                                        UPS_Origin_Longitude,
                                        UPS_False_Easting,
                                        UPS_False_Northing);

    Convert_Geodetic_To_Polar_Stereographic(Latitude,
                                            Longitude,
                                            &tempEasting,
                                            &tempNorthing);

    *Easting = tempEasting;
    *Northing = tempNorthing;
  }  /*  END of if(!Error_Code)   */

  return Error_Code;
}  /* END OF Convert_Geodetic_To_UPS  */

