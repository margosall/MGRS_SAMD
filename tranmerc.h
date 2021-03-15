#include <math.h>
#include "pi.h"

#define TRANMERC_NO_ERROR           0x0000
#define TRANMERC_LAT_ERROR          0x0001
#define TRANMERC_LON_ERROR          0x0002
#define TRANMERC_EASTING_ERROR      0x0004
#define TRANMERC_NORTHING_ERROR     0x0008
#define TRANMERC_ORIGIN_LAT_ERROR   0x0010
#define TRANMERC_CENT_MER_ERROR     0x0020
#define TRANMERC_A_ERROR            0x0040
#define TRANMERC_INV_F_ERROR        0x0080
#define TRANMERC_SCALE_FACTOR_ERROR 0x0100
#define TRANMERC_LON_WARNING        0x0200

#define MAX_DELTA_LONG  ((PI * 90)/180.0)       /* 90 degrees in radians */
#define MIN_SCALE_FACTOR  0.3
#define MAX_SCALE_FACTOR  3.0

#define SPHTMD(Latitude) ((double) (TranMerc_ap * Latitude \
      - TranMerc_bp * sin(2.e0 * Latitude) + TranMerc_cp * sin(4.e0 * Latitude) \
      - TranMerc_dp * sin(6.e0 * Latitude) + TranMerc_ep * sin(8.e0 * Latitude) ) )

#define SPHSN(Latitude) ((double) (TranMerc_a / sqrt( 1.e0 - TranMerc_es * \
      pow(sin(Latitude), 2))))

#define SPHSR(Latitude) ((double) (TranMerc_a * (1.e0 - TranMerc_es) / \
    pow(DENOM(Latitude), 3)))

#define DENOM(Latitude) ((double) (sqrt(1.e0 - TranMerc_es * pow(sin(Latitude),2))))


/**************************************************************************/
/*                               GLOBAL DECLARATIONS
 *
 */

/* Ellipsoid Parameters, default to WGS 84  */
static double TranMerc_a = 6378137.0;              /* Semi-major axis of ellipsoid in meters */
static double TranMerc_f = 1 / 298.257223563;      /* Flattening of ellipsoid  */
static double TranMerc_es = 0.0066943799901413800; /* Eccentricity (0.08181919084262188000) squared */
static double TranMerc_ebs = 0.0067394967565869;   /* Second Eccentricity squared */

/* Transverse_Mercator projection Parameters */
static double TranMerc_Origin_Lat = 0.0;           /* Latitude of origin in radians */
static double TranMerc_Origin_Long = 0.0;          /* Longitude of origin in radians */
static double TranMerc_False_Northing = 0.0;       /* False northing in meters */
static double TranMerc_False_Easting = 0.0;        /* False easting in meters */
static double TranMerc_Scale_Factor = 1.0;         /* Scale factor  */

/* Isometeric to geodetic latitude parameters, default to WGS 84 */
static double TranMerc_ap = 6367449.1458008;
static double TranMerc_bp = 16038.508696861;
static double TranMerc_cp = 16.832613334334;
static double TranMerc_dp = 0.021984404273757;
static double TranMerc_ep = 3.1148371319283e-005;

/* Maximum variance for easting and northing values for WGS 84. */
static double TranMerc_Delta_Easting = 40000000.0;
static double TranMerc_Delta_Northing = 40000000.0;

long Convert_Geodetic_To_Transverse_Mercator (double Latitude,
                                              double Longitude,
                                              double *Easting,
                                              double *Northing)

{      /* BEGIN Convert_Geodetic_To_Transverse_Mercator */

  /*
   * The function Convert_Geodetic_To_Transverse_Mercator converts geodetic
   * (latitude and longitude) coordinates to Transverse Mercator projection
   * (easting and northing) coordinates, according to the current ellipsoid
   * and Transverse Mercator projection coordinates.  If any errors occur, the
   * error code(s) are returned by the function, otherwise TRANMERC_NO_ERROR is
   * returned.
   *
   *    Latitude      : Latitude in radians                         (input)
   *    Longitude     : Longitude in radians                        (input)
   *    Easting       : Easting/X in meters                         (output)
   *    Northing      : Northing/Y in meters                        (output)
   */

  double c;       /* Cosine of latitude                          */
  double c2;
  double c3;
  double c5;
  double c7;
  double dlam;    /* Delta longitude - Difference in Longitude       */
  double eta;     /* constant - TranMerc_ebs *c *c                   */
  double eta2;
  double eta3;
  double eta4;
  double s;       /* Sine of latitude                        */
  double sn;      /* Radius of curvature in the prime vertical       */
  double t;       /* Tangent of latitude                             */
  double tan2;
  double tan3;
  double tan4;
  double tan5;
  double tan6;
  double t1;      /* Term in coordinate conversion formula - GP to Y */
  double t2;      /* Term in coordinate conversion formula - GP to Y */
  double t3;      /* Term in coordinate conversion formula - GP to Y */
  double t4;      /* Term in coordinate conversion formula - GP to Y */
  double t5;      /* Term in coordinate conversion formula - GP to Y */
  double t6;      /* Term in coordinate conversion formula - GP to Y */
  double t7;      /* Term in coordinate conversion formula - GP to Y */
  double t8;      /* Term in coordinate conversion formula - GP to Y */
  double t9;      /* Term in coordinate conversion formula - GP to Y */
  double tmd;     /* True Meridional distance                        */
  double tmdo;    /* True Meridional distance for latitude of origin */
  long    Error_Code = TRANMERC_NO_ERROR;
  double temp_Origin;
  double temp_Long;

  if ((Latitude < -MAX_LAT) || (Latitude > MAX_LAT))
  {  /* Latitude out of range */
    Error_Code|= TRANMERC_LAT_ERROR;
  }
  if (Longitude > PI)
    Longitude -= (2 * PI);
  if ((Longitude < (TranMerc_Origin_Long - MAX_DELTA_LONG))
      || (Longitude > (TranMerc_Origin_Long + MAX_DELTA_LONG)))
  {
    if (Longitude < 0)
      temp_Long = Longitude + 2 * PI;
    else
      temp_Long = Longitude;
    if (TranMerc_Origin_Long < 0)
      temp_Origin = TranMerc_Origin_Long + 2 * PI;
    else
      temp_Origin = TranMerc_Origin_Long;
    if ((temp_Long < (temp_Origin - MAX_DELTA_LONG))
        || (temp_Long > (temp_Origin + MAX_DELTA_LONG)))
      Error_Code|= TRANMERC_LON_ERROR;
  }
  if (!Error_Code)
  { /* no errors */

    /* 
     *  Delta Longitude
     */
    dlam = Longitude - TranMerc_Origin_Long;

    if (fabs(dlam) > (9.0 * PI / 180))
    { /* Distortion will result if Longitude is more than 9 degrees from the Central Meridian */
      Error_Code |= TRANMERC_LON_WARNING;
    }

    if (dlam > PI)
      dlam -= (2 * PI);
    if (dlam < -PI)
      dlam += (2 * PI);
    if (fabs(dlam) < 2.e-10)
      dlam = 0.0;

    s = sin(Latitude);
    c = cos(Latitude);
    c2 = c * c;
    c3 = c2 * c;
    c5 = c3 * c2;
    c7 = c5 * c2;
    t = tan (Latitude);
    tan2 = t * t;
    tan3 = tan2 * t;
    tan4 = tan3 * t;
    tan5 = tan4 * t;
    tan6 = tan5 * t;
    eta = TranMerc_ebs * c2;
    eta2 = eta * eta;
    eta3 = eta2 * eta;
    eta4 = eta3 * eta;

    /* radius of curvature in prime vertical */
    sn = SPHSN(Latitude);

    /* True Meridianal Distances */
    tmd = SPHTMD(Latitude);

    /*  Origin  */
    tmdo = SPHTMD (TranMerc_Origin_Lat);

    /* northing */
    t1 = (tmd - tmdo) * TranMerc_Scale_Factor;
    t2 = sn * s * c * TranMerc_Scale_Factor/ 2.e0;
    t3 = sn * s * c3 * TranMerc_Scale_Factor * (5.e0 - tan2 + 9.e0 * eta 
                                                + 4.e0 * eta2) /24.e0; 

    t4 = sn * s * c5 * TranMerc_Scale_Factor * (61.e0 - 58.e0 * tan2
                                                + tan4 + 270.e0 * eta - 330.e0 * tan2 * eta + 445.e0 * eta2
                                                + 324.e0 * eta3 -680.e0 * tan2 * eta2 + 88.e0 * eta4 
                                                -600.e0 * tan2 * eta3 - 192.e0 * tan2 * eta4) / 720.e0;

    t5 = sn * s * c7 * TranMerc_Scale_Factor * (1385.e0 - 3111.e0 * 
                                                tan2 + 543.e0 * tan4 - tan6) / 40320.e0;

    *Northing = TranMerc_False_Northing + t1 + pow(dlam,2.e0) * t2
                + pow(dlam,4.e0) * t3 + pow(dlam,6.e0) * t4
                + pow(dlam,8.e0) * t5; 

    /* Easting */
    t6 = sn * c * TranMerc_Scale_Factor;
    t7 = sn * c3 * TranMerc_Scale_Factor * (1.e0 - tan2 + eta ) /6.e0;
    t8 = sn * c5 * TranMerc_Scale_Factor * (5.e0 - 18.e0 * tan2 + tan4
                                            + 14.e0 * eta - 58.e0 * tan2 * eta + 13.e0 * eta2 + 4.e0 * eta3 
                                            - 64.e0 * tan2 * eta2 - 24.e0 * tan2 * eta3 )/ 120.e0;
    t9 = sn * c7 * TranMerc_Scale_Factor * ( 61.e0 - 479.e0 * tan2
                                             + 179.e0 * tan4 - tan6 ) /5040.e0;

    *Easting = TranMerc_False_Easting + dlam * t6 + pow(dlam,3.e0) * t7 
               + pow(dlam,5.e0) * t8 + pow(dlam,7.e0) * t9;
  }
  return (Error_Code);
} /* END OF Convert_Geodetic_To_Transverse_Mercator */

long Set_Transverse_Mercator_Parameters(double a,
                                        double f,
                                        double Origin_Latitude,
                                        double Central_Meridian,
                                        double False_Easting,
                                        double False_Northing,
                                        double Scale_Factor)

{ /* BEGIN Set_Tranverse_Mercator_Parameters */
  /*
   * The function Set_Tranverse_Mercator_Parameters receives the ellipsoid
   * parameters and Tranverse Mercator projection parameters as inputs, and
   * sets the corresponding state variables. If any errors occur, the error
   * code(s) are returned by the function, otherwise TRANMERC_NO_ERROR is
   * returned.
   *
   *    a                 : Semi-major axis of ellipsoid, in meters    (input)
   *    f                 : Flattening of ellipsoid                     (input)
   *    Origin_Latitude   : Latitude in radians at the origin of the   (input)
   *                         projection
   *    Central_Meridian  : Longitude in radians at the center of the  (input)
   *                         projection
   *    False_Easting     : Easting/X at the center of the projection  (input)
   *    False_Northing    : Northing/Y at the center of the projection (input)
   *    Scale_Factor      : Projection scale factor                    (input) 
   */

  double tn;        /* True Meridianal distance constant  */
  double tn2;
  double tn3;
  double tn4;
  double tn5;
  double dummy_northing;
  double TranMerc_b; /* Semi-minor axis of ellipsoid, in meters */
  double inv_f = 1 / f;
  long Error_Code = TRANMERC_NO_ERROR;

  if (a <= 0.0)
  { /* Semi-major axis must be greater than zero */
    Error_Code |= TRANMERC_A_ERROR;
  }
  if ((inv_f < 250) || (inv_f > 350))
  { /* Inverse flattening must be between 250 and 350 */
    Error_Code |= TRANMERC_INV_F_ERROR;
  }
  if ((Origin_Latitude < -PI_OVER_2) || (Origin_Latitude > PI_OVER_2))
  { /* origin latitude out of range */
    Error_Code |= TRANMERC_ORIGIN_LAT_ERROR;
  }
  if ((Central_Meridian < -PI) || (Central_Meridian > (2*PI)))
  { /* origin longitude out of range */
    Error_Code |= TRANMERC_CENT_MER_ERROR;
  }
  if ((Scale_Factor < MIN_SCALE_FACTOR) || (Scale_Factor > MAX_SCALE_FACTOR))
  {
    Error_Code |= TRANMERC_SCALE_FACTOR_ERROR;
  }
  if (!Error_Code)
  { /* no errors */
    TranMerc_a = a;
    TranMerc_f = f;
    TranMerc_Origin_Lat = Origin_Latitude;
    if (Central_Meridian > PI)
      Central_Meridian -= (2*PI);
    TranMerc_Origin_Long = Central_Meridian;
    TranMerc_False_Northing = False_Northing;
    TranMerc_False_Easting = False_Easting; 
    TranMerc_Scale_Factor = Scale_Factor;

    /* Eccentricity Squared */
    TranMerc_es = 2 * TranMerc_f - TranMerc_f * TranMerc_f;
    /* Second Eccentricity Squared */
    TranMerc_ebs = (1 / (1 - TranMerc_es)) - 1;

    TranMerc_b = TranMerc_a * (1 - TranMerc_f);    
    /*True meridianal constants  */
    tn = (TranMerc_a - TranMerc_b) / (TranMerc_a + TranMerc_b);
    tn2 = tn * tn;
    tn3 = tn2 * tn;
    tn4 = tn3 * tn;
    tn5 = tn4 * tn;

    TranMerc_ap = TranMerc_a * (1.e0 - tn + 5.e0 * (tn2 - tn3)/4.e0
                                + 81.e0 * (tn4 - tn5)/64.e0 );
    TranMerc_bp = 3.e0 * TranMerc_a * (tn - tn2 + 7.e0 * (tn3 - tn4)
                                       /8.e0 + 55.e0 * tn5/64.e0 )/2.e0;
    TranMerc_cp = 15.e0 * TranMerc_a * (tn2 - tn3 + 3.e0 * (tn4 - tn5 )/4.e0) /16.0;
    TranMerc_dp = 35.e0 * TranMerc_a * (tn3 - tn4 + 11.e0 * tn5 / 16.e0) / 48.e0;
    TranMerc_ep = 315.e0 * TranMerc_a * (tn4 - tn5) / 512.e0;
    Convert_Geodetic_To_Transverse_Mercator(MAX_LAT,
                                            MAX_DELTA_LONG + Central_Meridian,
                                            &TranMerc_Delta_Easting,
                                            &TranMerc_Delta_Northing);
    Convert_Geodetic_To_Transverse_Mercator(0,
                                            MAX_DELTA_LONG + Central_Meridian,
                                            &TranMerc_Delta_Easting,
                                            &dummy_northing);
    TranMerc_Delta_Northing++;
    TranMerc_Delta_Easting++;

  } /* END OF if(!Error_Code) */
  return (Error_Code);
}  /* END of Set_Transverse_Mercator_Parameters  */
