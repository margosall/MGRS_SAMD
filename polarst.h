#define POLAR_NO_ERROR                0x0000
#define POLAR_LAT_ERROR               0x0001
#define POLAR_LON_ERROR               0x0002
#define POLAR_ORIGIN_LAT_ERROR        0x0004
#define POLAR_ORIGIN_LON_ERROR        0x0008
#define POLAR_EASTING_ERROR        0x0010
#define POLAR_NORTHING_ERROR      0x0020
#define POLAR_A_ERROR                 0x0040
#define POLAR_INV_F_ERROR             0x0080
#define POLAR_RADIUS_ERROR            0x0100

#define PI           3.14159265358979323e0       /* PI     */
#define PI_OVER_2    (PI / 2.0)           
#define TWO_PI       (2.0 * PI)
#define POLAR_POW(EsSin)     pow((1.0 - EsSin) / (1.0 + EsSin), es_OVER_2)

/************************************************************************/
/*                           GLOBAL DECLARATIONS
 *
 */

const double PI_Over_4 = (PI / 4.0);

/* Ellipsoid Parameters, default to WGS 84  */
static double Polar_a = 6378137.0;                    /* Semi-major axis of ellipsoid in meters  */
static double Polar_f = 1 / 298.257223563;            /* Flattening of ellipsoid  */
static double es = 0.08181919084262188000;            /* Eccentricity of ellipsoid    */
static double es_OVER_2 = .040909595421311;           /* es / 2.0 */
static double Southern_Hemisphere = 0;                /* Flag variable */
static double tc = 1.0;
static double e4 = 1.0033565552493;
static double Polar_a_mc = 6378137.0;                 /* Polar_a * mc */
static double two_Polar_a = 12756274.0;               /* 2.0 * Polar_a */

/* Polar Stereographic projection Parameters */
static double Polar_Origin_Lat = ((PI * 90) / 180);   /* Latitude of origin in radians */
static double Polar_Origin_Long = 0.0;                /* Longitude of origin in radians */
static double Polar_False_Easting = 0.0;              /* False easting in meters */
static double Polar_False_Northing = 0.0;             /* False northing in meters */

/* Maximum variance for easting and northing values for WGS 84. */
static double Polar_Delta_Easting = 12713601.0;
static double Polar_Delta_Northing = 12713601.0;
long Convert_Geodetic_To_Polar_Stereographic (double Latitude,
                                              double Longitude,
                                              double *Easting,
                                              double *Northing)

{  /* BEGIN Convert_Geodetic_To_Polar_Stereographic */

/*
 * The function Convert_Geodetic_To_Polar_Stereographic converts geodetic
 * coordinates (latitude and longitude) to Polar Stereographic coordinates
 * (easting and northing), according to the current ellipsoid
 * and Polar Stereographic projection parameters. If any errors occur, error
 * code(s) are returned by the function, otherwise POLAR_NO_ERROR is returned.
 *
 *    Latitude   :  Latitude, in radians                      (input)
 *    Longitude  :  Longitude, in radians                     (input)
 *    Easting    :  Easting (X), in meters                    (output)
 *    Northing   :  Northing (Y), in meters                   (output)
 */

  double dlam;
  double slat;
  double essin;
  double t;
  double rho;
  double pow_es;
  long Error_Code = POLAR_NO_ERROR;

  if ((Latitude < -PI_OVER_2) || (Latitude > PI_OVER_2))
  {   /* Latitude out of range */
    Error_Code |= POLAR_LAT_ERROR;
  }
  if ((Latitude < 0) && (Southern_Hemisphere == 0))
  {   /* Latitude and Origin Latitude in different hemispheres */
    Error_Code |= POLAR_LAT_ERROR;
  }
  if ((Latitude > 0) && (Southern_Hemisphere == 1))
  {   /* Latitude and Origin Latitude in different hemispheres */
    Error_Code |= POLAR_LAT_ERROR;
  }
  if ((Longitude < -PI) || (Longitude > TWO_PI))
  {  /* Longitude out of range */
    Error_Code |= POLAR_LON_ERROR;
  }


  if (!Error_Code)
  {  /* no errors */

    if (fabs(fabs(Latitude) - PI_OVER_2) < 1.0e-10)
    {
      *Easting = Polar_False_Easting;
      *Northing = Polar_False_Northing;
    }
    else
    {
      if (Southern_Hemisphere != 0)
      {
        Longitude *= -1.0;
        Latitude *= -1.0;
      }
      dlam = Longitude - Polar_Origin_Long;
      if (dlam > PI)
      {
        dlam -= TWO_PI;
      }
      if (dlam < -PI)
      {
        dlam += TWO_PI;
      }
      slat = sin(Latitude);
      essin = es * slat;
      pow_es = POLAR_POW(essin);
      t = tan(PI_Over_4 - Latitude / 2.0) / pow_es;

      if (fabs(fabs(Polar_Origin_Lat) - PI_OVER_2) > 1.0e-10)
        rho = Polar_a_mc * t / tc;
      else
        rho = two_Polar_a * t / e4;


      if (Southern_Hemisphere != 0)
      {
        *Easting = -(rho * sin(dlam) - Polar_False_Easting);
     //   *Easting *= -1.0;
        *Northing = rho * cos(dlam) + Polar_False_Northing;
      }
      else
      {
        *Easting = rho * sin(dlam) + Polar_False_Easting;
        *Northing = -rho * cos(dlam) + Polar_False_Northing;
      }

    }
  }
  return (Error_Code);
} /* END OF Convert_Geodetic_To_Polar_Stereographic */

long Set_Polar_Stereographic_Parameters (double a,
                                         double f,
                                         double Latitude_of_True_Scale,
                                         double Longitude_Down_from_Pole,
                                         double False_Easting,
                                         double False_Northing)

{  /* BEGIN Set_Polar_Stereographic_Parameters   */
/*  
 *  The function Set_Polar_Stereographic_Parameters receives the ellipsoid
 *  parameters and Polar Stereograpic projection parameters as inputs, and
 *  sets the corresponding state variables.  If any errors occur, error
 *  code(s) are returned by the function, otherwise POLAR_NO_ERROR is returned.
 *
 *  a                : Semi-major axis of ellipsoid, in meters         (input)
 *  f                : Flattening of ellipsoid                         (input)
 *  Latitude_of_True_Scale  : Latitude of true scale, in radians       (input)
 *  Longitude_Down_from_Pole : Longitude down from pole, in radians    (input)
 *  False_Easting    : Easting (X) at center of projection, in meters  (input)
 *  False_Northing   : Northing (Y) at center of projection, in meters (input)
 */

  double es2;
  double slat, clat;
  double essin;
  double one_PLUS_es, one_MINUS_es;
  double pow_es;
  double temp, temp_northing;
  double inv_f = 1 / f;
  double mc;                    
//  const double  epsilon = 1.0e-2;
  long Error_Code = POLAR_NO_ERROR;

  if (a <= 0.0)
  { /* Semi-major axis must be greater than zero */
    Error_Code |= POLAR_A_ERROR;
  }
  if ((inv_f < 250) || (inv_f > 350))
  { /* Inverse flattening must be between 250 and 350 */
    Error_Code |= POLAR_INV_F_ERROR;
  }
  if ((Latitude_of_True_Scale < -PI_OVER_2) || (Latitude_of_True_Scale > PI_OVER_2))
  { /* Origin Latitude out of range */
    Error_Code |= POLAR_ORIGIN_LAT_ERROR;
  }
  if ((Longitude_Down_from_Pole < -PI) || (Longitude_Down_from_Pole > TWO_PI))
  { /* Origin Longitude out of range */
    Error_Code |= POLAR_ORIGIN_LON_ERROR;
  }

  if (!Error_Code)
  { /* no errors */

    Polar_a = a;
    two_Polar_a = 2.0 * Polar_a;
    Polar_f = f;

    if (Longitude_Down_from_Pole > PI)
      Longitude_Down_from_Pole -= TWO_PI;
    if (Latitude_of_True_Scale < 0)
    {
      Southern_Hemisphere = 1;
      Polar_Origin_Lat = -Latitude_of_True_Scale;
      Polar_Origin_Long = -Longitude_Down_from_Pole;
    }
    else
    {
      Southern_Hemisphere = 0;
      Polar_Origin_Lat = Latitude_of_True_Scale;
      Polar_Origin_Long = Longitude_Down_from_Pole;
    }
    Polar_False_Easting = False_Easting;
    Polar_False_Northing = False_Northing;

    es2 = 2 * Polar_f - Polar_f * Polar_f;
    es = sqrt(es2);
    es_OVER_2 = es / 2.0;

    if (fabs(fabs(Polar_Origin_Lat) - PI_OVER_2) > 1.0e-10)
    {
      slat = sin(Polar_Origin_Lat);
      essin = es * slat;
      pow_es = POLAR_POW(essin);
      clat = cos(Polar_Origin_Lat);
      mc = clat / sqrt(1.0 - essin * essin);
      Polar_a_mc = Polar_a * mc;
      tc = tan(PI_Over_4 - Polar_Origin_Lat / 2.0) / pow_es;
    }
    else
    {
      one_PLUS_es = 1.0 + es;
      one_MINUS_es = 1.0 - es;
      e4 = sqrt(pow(one_PLUS_es, one_PLUS_es) * pow(one_MINUS_es, one_MINUS_es));
    }

    /* Calculate Radius */
    Convert_Geodetic_To_Polar_Stereographic(0, Longitude_Down_from_Pole, 
                                            &temp, &temp_northing);

    Polar_Delta_Northing = temp_northing;
    if(Polar_False_Northing)
      Polar_Delta_Northing -= Polar_False_Northing;
    if (Polar_Delta_Northing < 0)
      Polar_Delta_Northing = -Polar_Delta_Northing;
    Polar_Delta_Northing *= 1.01;

    Polar_Delta_Easting = Polar_Delta_Northing;

  /*  Polar_Delta_Easting = temp_northing;
    if(Polar_False_Easting)
      Polar_Delta_Easting -= Polar_False_Easting;
    if (Polar_Delta_Easting < 0)
      Polar_Delta_Easting = -Polar_Delta_Easting;
    Polar_Delta_Easting *= 1.01;*/
  }

  return (Error_Code);
} /* END OF Set_Polar_Stereographic_Parameters */


