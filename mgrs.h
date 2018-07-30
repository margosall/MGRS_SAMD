#include "ups.h"
#include "utm.h"
#define PI180 0.01745329251
#define DEG_TO_RAD       0.017453292519943295 /* PI/180                      */
#define RAD_TO_DEG       57.29577951308232087 /* 180/PI                      */
#define LETTER_A               0   /* ARRAY INDEX FOR LETTER A               */
#define LETTER_B               1   /* ARRAY INDEX FOR LETTER B               */
#define LETTER_C               2   /* ARRAY INDEX FOR LETTER C               */
#define LETTER_D               3   /* ARRAY INDEX FOR LETTER D               */
#define LETTER_E               4   /* ARRAY INDEX FOR LETTER E               */
#define LETTER_F               5   /* ARRAY INDEX FOR LETTER F               */
#define LETTER_G               6   /* ARRAY INDEX FOR LETTER G               */
#define LETTER_H               7   /* ARRAY INDEX FOR LETTER H               */
#define LETTER_I               8   /* ARRAY INDEX FOR LETTER I               */
#define LETTER_J               9   /* ARRAY INDEX FOR LETTER J               */
#define LETTER_K              10   /* ARRAY INDEX FOR LETTER K               */
#define LETTER_L              11   /* ARRAY INDEX FOR LETTER L               */
#define LETTER_M              12   /* ARRAY INDEX FOR LETTER M               */
#define LETTER_N              13   /* ARRAY INDEX FOR LETTER N               */
#define LETTER_O              14   /* ARRAY INDEX FOR LETTER O               */
#define LETTER_P              15   /* ARRAY INDEX FOR LETTER P               */
#define LETTER_Q              16   /* ARRAY INDEX FOR LETTER Q               */
#define LETTER_R              17   /* ARRAY INDEX FOR LETTER R               */
#define LETTER_S              18   /* ARRAY INDEX FOR LETTER S               */
#define LETTER_T              19   /* ARRAY INDEX FOR LETTER T               */
#define LETTER_U              20   /* ARRAY INDEX FOR LETTER U               */
#define LETTER_V              21   /* ARRAY INDEX FOR LETTER V               */
#define LETTER_W              22   /* ARRAY INDEX FOR LETTER W               */
#define LETTER_X              23   /* ARRAY INDEX FOR LETTER X               */
#define LETTER_Y              24   /* ARRAY INDEX FOR LETTER Y               */
#define LETTER_Z              25   /* ARRAY INDEX FOR LETTER Z               */
#define MGRS_LETTERS            3  /* NUMBER OF LETTERS IN MGRS              */
#define ONEHT          100000.e0    /* ONE HUNDRED THOUSAND                  */
#define TWOMIL        2000000.e0    /* TWO MILLION                           */
#define TRUE                      1  /* CONSTANT VALUE FOR TRUE VALUE  */
#define FALSE                     0  /* CONSTANT VALUE FOR FALSE VALUE */
#define PI    3.14159265358979323e0  /* PI                             */
#define PI_OVER_2  (PI / 2.0e0)


#define MIN_EASTING  100000
#define MAX_EASTING  900000
#define MIN_NORTHING 0
#define MAX_NORTHING 10000000
#define MAX_PRECISION           5   /* Maximum precision of easting & northing */
#define MIN_UTM_LAT      ( (-80 * PI) / 180.0 ) /* -80 degrees in radians    */
#define MAX_UTM_LAT      ( (84 * PI) / 180.0 )  /* 84 degrees in radians     */

#define MIN_EAST_NORTH 0
#define MAX_EAST_NORTH 4000000
/* Ellipsoid parameters, default to WGS 84 */
double MGRS_a = 6378137.0;    /* Semi-major axis of ellipsoid in meters */
double MGRS_f = 1 / 298.257223563; /* Flattening of ellipsoid           */
char   MGRS_Ellipsoid_Code[3] = {'W','E',0};

/*
 *    CLARKE_1866 : Ellipsoid code for CLARKE_1866
 *    CLARKE_1880 : Ellipsoid code for CLARKE_1880
 *    BESSEL_1841 : Ellipsoid code for BESSEL_1841
 *    BESSEL_1841_NAMIBIA : Ellipsoid code for BESSEL 1841 (NAMIBIA)
 */
const char* CLARKE_1866 = "CC";
const char* CLARKE_1880 = "CD";
const char* BESSEL_1841 = "BR";
const char* BESSEL_1841_NAMIBIA = "BN";

#define MGRS_NO_ERROR                0x0000
#define MGRS_LAT_ERROR               0x0001
#define MGRS_LON_ERROR               0x0002
#define MGRS_STRING_ERROR            0x0004
#define MGRS_PRECISION_ERROR         0x0008
#define MGRS_A_ERROR                 0x0010
#define MGRS_INV_F_ERROR             0x0020
#define MGRS_EASTING_ERROR           0x0040
#define MGRS_NORTHING_ERROR          0x0080
#define MGRS_ZONE_ERROR              0x0100
#define MGRS_HEMISPHERE_ERROR        0x0200
#define MGRS_LAT_WARNING             0x0400

typedef struct Latitude_Band_Value
{
  long letter;            /* letter representing latitude band  */
  double min_northing;    /* minimum northing for latitude band */
  double north;           /* upper latitude for latitude band   */
  double south;           /* lower latitude for latitude band   */
  double northing_offset; /* latitude band northing offset      */
} Latitude_Band;

static const Latitude_Band Latitude_Band_Table[20] =
  {{LETTER_C, 1100000.0, -72.0, -80.5, 0.0},
  {LETTER_D, 2000000.0, -64.0, -72.0, 2000000.0},
  {LETTER_E, 2800000.0, -56.0, -64.0, 2000000.0},
  {LETTER_F, 3700000.0, -48.0, -56.0, 2000000.0},
  {LETTER_G, 4600000.0, -40.0, -48.0, 4000000.0},
  {LETTER_H, 5500000.0, -32.0, -40.0, 4000000.0},
  {LETTER_J, 6400000.0, -24.0, -32.0, 6000000.0},
  {LETTER_K, 7300000.0, -16.0, -24.0, 6000000.0},
  {LETTER_L, 8200000.0, -8.0, -16.0, 8000000.0},
  {LETTER_M, 9100000.0, 0.0, -8.0, 8000000.0},
  {LETTER_N, 0.0, 8.0, 0.0, 0.0},
  {LETTER_P, 800000.0, 16.0, 8.0, 0.0},
  {LETTER_Q, 1700000.0, 24.0, 16.0, 0.0},
  {LETTER_R, 2600000.0, 32.0, 24.0, 2000000.0},
  {LETTER_S, 3500000.0, 40.0, 32.0, 2000000.0},
  {LETTER_T, 4400000.0, 48.0, 40.0, 4000000.0},
  {LETTER_U, 5300000.0, 56.0, 48.0, 4000000.0},
  {LETTER_V, 6200000.0, 64.0, 56.0, 6000000.0},
  {LETTER_W, 7000000.0, 72.0, 64.0, 6000000.0},
  {LETTER_X, 7900000.0, 84.5, 72.0, 6000000.0}};
  
typedef struct UPS_Constant_Value
{
  long letter;            /* letter representing latitude band      */
  long ltr2_low_value;    /* 2nd letter range - low number         */
  long ltr2_high_value;   /* 2nd letter range - high number          */
  long ltr3_high_value;   /* 3rd letter range - high number (UPS)   */
  double false_easting;   /* False easting based on 2nd letter      */
  double false_northing;  /* False northing based on 3rd letter     */
} UPS_Constant;

static const UPS_Constant UPS_Constant_Table[4] =
  {{LETTER_A, LETTER_J, LETTER_Z, LETTER_Z, 800000.0, 800000.0},
  {LETTER_B, LETTER_A, LETTER_R, LETTER_Z, 2000000.0, 800000.0},
  {LETTER_Y, LETTER_J, LETTER_Z, LETTER_P, 800000.0, 1300000.0},
  {LETTER_Z, LETTER_A, LETTER_J, LETTER_P, 2000000.0, 1300000.0}};
  
long Make_MGRS_String (char* MGRS,
                       long Zone,
                       int Letters[MGRS_LETTERS],
                       double Easting,
                       double Northing,
                       long Precision)
/*
 * The function Make_MGRS_String constructs an MGRS string
 * from its component parts.
 *
 *   MGRS           : MGRS coordinate string          (output)
 *   Zone           : UTM Zone                        (input)
 *   Letters        : MGRS coordinate string letters  (input)
 *   Easting        : Easting value                   (input)
 *   Northing       : Northing value                  (input)
 *   Precision      : Precision level of MGRS string  (input)
 */
{ /* Make_MGRS_String */
  long i;
  long j;
  double divisor;
  long east;
  long north;
  char alphabet[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  long error_code = MGRS_NO_ERROR;

  i = 0;
  if (Zone)
    i = sprintf (MGRS+i,"%2.2ld",Zone);
  else
    strncpy(MGRS, "  ", 2);  // 2 spaces

  for (j=0;j<3;j++)
    MGRS[i++] = alphabet[Letters[j]];
  divisor = pow (10.0, (5 - Precision));
  Easting = fmod (Easting, 100000.0);
  if (Easting >= 99999.5)
    Easting = 99999.0;
  east = (long)(Easting/divisor);
  i += sprintf (MGRS+i, "%*.*ld", (int)Precision, (int)Precision, east);
  Northing = fmod (Northing, 100000.0);
  if (Northing >= 99999.5)
    Northing = 99999.0;
  north = (long)(Northing/divisor);
  i += sprintf (MGRS+i, "%*.*ld", (int)Precision, (int)Precision, north);
  return (error_code);
} /* Make_MGRS_String */

void Get_Grid_Values (long zone,
                      long* ltr2_low_value,
                      long* ltr2_high_value,
                      double *pattern_offset)
/*
 * The function getGridValues sets the letter range used for
 * the 2nd letter in the MGRS coordinate string, based on the set
 * number of the utm zone. It also sets the pattern offset using a
 * value of A for the second letter of the grid square, based on
 * the grid pattern and set number of the utm zone.
 *
 *    zone            : Zone number             (input)
 *    ltr2_low_value  : 2nd letter low number   (output)
 *    ltr2_high_value : 2nd letter high number  (output)
 *    pattern_offset  : Pattern offset          (output)
 */
{ /* BEGIN Get_Grid_Values */
  long set_number;    /* Set number (1-6) based on UTM zone number */
  long aa_pattern;    /* Pattern based on ellipsoid code */

  set_number = zone % 6;

  if (!set_number)
    set_number = 6;

  if (!strcmp(MGRS_Ellipsoid_Code,CLARKE_1866) || !strcmp(MGRS_Ellipsoid_Code, CLARKE_1880) ||
      !strcmp(MGRS_Ellipsoid_Code,BESSEL_1841) || !strcmp(MGRS_Ellipsoid_Code,BESSEL_1841_NAMIBIA))
    aa_pattern = FALSE;
  else
    aa_pattern = TRUE;

  if ((set_number == 1) || (set_number == 4))
  {
    *ltr2_low_value = LETTER_A;
    *ltr2_high_value = LETTER_H;
  }
  else if ((set_number == 2) || (set_number == 5))
  {
    *ltr2_low_value = LETTER_J;
    *ltr2_high_value = LETTER_R;
  }
  else if ((set_number == 3) || (set_number == 6))
  {
    *ltr2_low_value = LETTER_S;
    *ltr2_high_value = LETTER_Z;
  }

  /* False northing at A for second letter of grid square */
  if (aa_pattern)
  {
    if ((set_number % 2) ==  0)
      *pattern_offset = 500000.0;
    else
      *pattern_offset = 0.0;
  }
  else
  {
    if ((set_number % 2) == 0)
      *pattern_offset =  1500000.0;
    else
      *pattern_offset = 1000000.00;
  }
} /* END OF Get_Grid_Values */

long Get_Latitude_Letter(double latitude, int* letter)
/*
 * The function Get_Latitude_Letter receives a latitude value
 * and uses the Latitude_Band_Table to determine the latitude band
 * letter for that latitude.
 *
 *   latitude   : Latitude              (input)
 *   letter     : Latitude band letter  (output)
 */
{ /* Get_Latitude_Letter */
  double temp = 0.0;
  long error_code = MGRS_NO_ERROR;
  double lat_deg = latitude * RAD_TO_DEG;

  if (lat_deg >= 72 && lat_deg < 84.5)
    *letter = LETTER_X;
  else if (lat_deg > -80.5 && lat_deg < 72)
  {
    temp = ((latitude + (80.0 * DEG_TO_RAD)) / (8.0 * DEG_TO_RAD)) + 1.0e-12;
    *letter = Latitude_Band_Table[(int)temp].letter;
  }
  else
    error_code |= MGRS_LAT_ERROR;

  return error_code;
} /* Get_Latitude_Letter */

long UTM_To_MGRS (long Zone,
                  char Hemisphere,
                  double Longitude,
                  double Latitude,
                  double Easting,
                  double Northing,
                  long Precision,
                  char *MGRS)
/*
 * The function UTM_To_MGRS calculates an MGRS coordinate string
 * based on the zone, latitude, easting and northing.
 *
 *    Zone      : Zone number             (input)
 *    Hemisphere: Hemisphere              (input)
 *    Longitude : Longitude in radians    (input)
 *    Latitude  : Latitude in radians     (input)
 *    Easting   : Easting                 (input)
 *    Northing  : Northing                (input)
 *    Precision : Precision               (input)
 *    MGRS      : MGRS coordinate string  (output)
 */
{ /* BEGIN UTM_To_MGRS */
  double pattern_offset;      /* Northing offset for 3rd letter               */
  double grid_easting;        /* Easting used to derive 2nd letter of MGRS   */
  double grid_northing;       /* Northing used to derive 3rd letter of MGRS  */
  long ltr2_low_value;        /* 2nd letter range - low number               */
  long ltr2_high_value;       /* 2nd letter range - high number              */
  int letters[MGRS_LETTERS];  /* Number location of 3 letters in alphabet    */
  long temp_error_code = MGRS_NO_ERROR;
  long error_code = MGRS_NO_ERROR;


  /* Special check for rounding to (truncated) eastern edge of zone 31V */
  if ((Zone == 31) && (((Latitude >= 56.0 * DEG_TO_RAD) && (Latitude < 64.0 * DEG_TO_RAD)) && ((Longitude >= 3.0 * DEG_TO_RAD) || (Easting >= 500000.0))))
  { /* Reconvert to UTM zone 32 */
    Set_UTM_Parameters (MGRS_a, MGRS_f, 32);
    temp_error_code = Convert_Geodetic_To_UTM (Latitude, Longitude, &Zone, &Hemisphere, &Easting, &Northing);
    if(temp_error_code)
    {
      if(temp_error_code & UTM_LAT_ERROR)
        error_code |= MGRS_LAT_ERROR;
      if(temp_error_code & UTM_LON_ERROR)
        error_code |= MGRS_LON_ERROR;
      if(temp_error_code & UTM_ZONE_OVERRIDE_ERROR)
        error_code |= MGRS_ZONE_ERROR;
      if(temp_error_code & UTM_EASTING_ERROR)
        error_code |= MGRS_EASTING_ERROR;
      if(temp_error_code & UTM_NORTHING_ERROR)
        error_code |= MGRS_NORTHING_ERROR;

      return error_code;
    }
  }

  if( Latitude <= 0.0 && Northing == 1.0e7)
  {
    Latitude = 0.0;
    Northing = 0.0;
  }

  Get_Grid_Values(Zone, &ltr2_low_value, &ltr2_high_value, &pattern_offset);

  error_code = Get_Latitude_Letter(Latitude, &letters[0]);

  if (!error_code)
  {
    grid_northing = Northing;

    while (grid_northing >= TWOMIL)
    {
      grid_northing = grid_northing - TWOMIL;
    }
    grid_northing = grid_northing + pattern_offset;
    if(grid_northing >= TWOMIL)
      grid_northing = grid_northing - TWOMIL;

    letters[2] = (long)(grid_northing / ONEHT);
    if (letters[2] > LETTER_H)
      letters[2] = letters[2] + 1;

    if (letters[2] > LETTER_N)
      letters[2] = letters[2] + 1;

    grid_easting = Easting;
    if (((letters[0] == LETTER_V) && (Zone == 31)) && (grid_easting == 500000.0))
      grid_easting = grid_easting - 1.0; /* SUBTRACT 1 METER */

    letters[1] = ltr2_low_value + ((long)(grid_easting / ONEHT) -1);
    if ((ltr2_low_value == LETTER_J) && (letters[1] > LETTER_N))
      letters[1] = letters[1] + 1;

    Make_MGRS_String (MGRS, Zone, letters, grid_easting, Northing, Precision);
  }
  return error_code;
} /* END UTM_To_MGRS */
  
long Convert_UPS_To_MGRS (char   Hemisphere,
                          double Easting,
                          double Northing,
                          long   Precision,
                          char*  MGRS)
/*
 *  The function Convert_UPS_To_MGRS converts UPS (hemisphere, easting,
 *  and northing) coordinates to an MGRS coordinate string according to
 *  the current ellipsoid parameters.  If any errors occur, the error
 *  code(s) are returned by the function, otherwise UPS_NO_ERROR is
 *  returned.
 *
 *    Hemisphere    : Hemisphere either 'N' or 'S'     (input)
 *    Easting       : Easting/X in meters              (input)
 *    Northing      : Northing/Y in meters             (input)
 *    Precision     : Precision level of MGRS string   (input)
 *    MGRS          : MGRS coordinate string           (output)
 */
{ /* Convert_UPS_To_MGRS */
  double false_easting;       /* False easting for 2nd letter                 */
  double false_northing;      /* False northing for 3rd letter                */
  double grid_easting;        /* Easting used to derive 2nd letter of MGRS    */
  double grid_northing;       /* Northing used to derive 3rd letter of MGRS   */
  long ltr2_low_value;        /* 2nd letter range - low number                */
  int letters[MGRS_LETTERS];  /* Number location of 3 letters in alphabet     */
  int index = 0;
  long error_code = MGRS_NO_ERROR;

  if ((Hemisphere != 'N') && (Hemisphere != 'S'))
    error_code |= MGRS_HEMISPHERE_ERROR;
  if ((Easting < MIN_EAST_NORTH) || (Easting > MAX_EAST_NORTH))
    error_code |= MGRS_EASTING_ERROR;
  if ((Northing < MIN_EAST_NORTH) || (Northing > MAX_EAST_NORTH))
    error_code |= MGRS_NORTHING_ERROR;
  if ((Precision < 0) || (Precision > MAX_PRECISION))
    error_code |= MGRS_PRECISION_ERROR;
  if (!error_code)
  {

    if (Hemisphere == 'N')
    {
      if (Easting >= TWOMIL)
        letters[0] = LETTER_Z;
      else
        letters[0] = LETTER_Y;

      index = letters[0] - 22;
      ltr2_low_value = UPS_Constant_Table[index].ltr2_low_value;
      false_easting = UPS_Constant_Table[index].false_easting;
      false_northing = UPS_Constant_Table[index].false_northing;
    }
    else
    {
      if (Easting >= TWOMIL)
        letters[0] = LETTER_B;
      else
        letters[0] = LETTER_A;

      ltr2_low_value = UPS_Constant_Table[letters[0]].ltr2_low_value;
      false_easting = UPS_Constant_Table[letters[0]].false_easting;
      false_northing = UPS_Constant_Table[letters[0]].false_northing;
    }

    grid_northing = Northing;
    grid_northing = grid_northing - false_northing;
    letters[2] = (long)(grid_northing / ONEHT);

    if (letters[2] > LETTER_H)
      letters[2] = letters[2] + 1;

    if (letters[2] > LETTER_N)
      letters[2] = letters[2] + 1;

    grid_easting = Easting;
    grid_easting = grid_easting - false_easting;
    letters[1] = ltr2_low_value + ((long)(grid_easting / ONEHT));

    if (Easting < TWOMIL)
    {
      if (letters[1] > LETTER_L)
        letters[1] = letters[1] + 3;

      if (letters[1] > LETTER_U)
        letters[1] = letters[1] + 2;
    }
    else
    {
      if (letters[1] > LETTER_C)
        letters[1] = letters[1] + 2;

      if (letters[1] > LETTER_H)
        letters[1] = letters[1] + 1;

      if (letters[1] > LETTER_L)
        letters[1] = letters[1] + 3;
    }

    Make_MGRS_String (MGRS, 0, letters, Easting, Northing, Precision);
  }
  return (error_code);
} /* Convert_UPS_To_MGRS */

long Convert_Geodetic_To_MGRS (double Latitude,
                               double Longitude,
                               long Precision,
                               char* MGRS)
/*
 * The function Convert_Geodetic_To_MGRS converts Geodetic (latitude and
 * longitude) coordinates to an MGRS coordinate string, according to the
 * current ellipsoid parameters.  If any errors occur, the error code(s)
 * are returned by the function, otherwise MGRS_NO_ERROR is returned.
 *
 *    Latitude   : Latitude in radians              (input)
 *    Longitude  : Longitude in radians             (input)
 *    Precision  : Precision level of MGRS string   (input)
 *    MGRS       : MGRS coordinate string           (output)
 *
 */
{ /* Convert_Geodetic_To_MGRS */
  long zone;
  char hemisphere;
  double easting;
  double northing;
  long temp_error_code = MGRS_NO_ERROR;
  long error_code = MGRS_NO_ERROR;

  if ((Latitude < -PI_OVER_2) || (Latitude > PI_OVER_2))
  { /* Latitude out of range */
    error_code |= MGRS_LAT_ERROR;
  }
  if ((Longitude < -PI) || (Longitude > (2*PI)))
  { /* Longitude out of range */
    error_code |= MGRS_LON_ERROR;
  }
  if ((Precision < 0) || (Precision > MAX_PRECISION))
    error_code |= MGRS_PRECISION_ERROR;
  if (!error_code)
  {
    if ((Latitude < MIN_UTM_LAT) || (Latitude > MAX_UTM_LAT))
    {
      temp_error_code = Set_UPS_Parameters (MGRS_a, MGRS_f);
      if(!temp_error_code)
      {
        temp_error_code = Convert_Geodetic_To_UPS (Latitude, Longitude, &hemisphere, &easting, &northing);
        if(!temp_error_code)
        {
          error_code |= Convert_UPS_To_MGRS (hemisphere, easting, northing, Precision, MGRS);
        }
        else
        {
          if(temp_error_code & UPS_LAT_ERROR)
            error_code |= MGRS_LAT_ERROR;
          if(temp_error_code & UPS_LON_ERROR)
            error_code |= MGRS_LON_ERROR;
        }
      }
      else
      {
        if(temp_error_code & UPS_A_ERROR)
          error_code |= MGRS_A_ERROR;
        if(temp_error_code & UPS_INV_F_ERROR)
          error_code |= MGRS_INV_F_ERROR;
      }
    }
    else
    {
      temp_error_code = Set_UTM_Parameters (MGRS_a, MGRS_f, 0);
      if(!temp_error_code)
      {
        temp_error_code = Convert_Geodetic_To_UTM (Latitude, Longitude, &zone, &hemisphere, &easting, &northing);
        if(!temp_error_code)
          error_code |= UTM_To_MGRS (zone, hemisphere, Longitude, Latitude, easting, northing, Precision, MGRS);
        else
        {
          if(temp_error_code & UTM_LAT_ERROR)
            error_code |= MGRS_LAT_ERROR;
          if(temp_error_code & UTM_LON_ERROR)
            error_code |= MGRS_LON_ERROR;
          if(temp_error_code & UTM_ZONE_OVERRIDE_ERROR)
            error_code |= MGRS_ZONE_ERROR;
          if(temp_error_code & UTM_EASTING_ERROR)
            error_code |= MGRS_EASTING_ERROR;
          if(temp_error_code & UTM_NORTHING_ERROR)
            error_code |= MGRS_NORTHING_ERROR;
        }
      }
      else
      {
        if(temp_error_code & UTM_A_ERROR)
          error_code |= MGRS_A_ERROR;
        if(temp_error_code & UTM_INV_F_ERROR)
          error_code |= MGRS_INV_F_ERROR;
        if(temp_error_code & UTM_ZONE_OVERRIDE_ERROR)
          error_code |= MGRS_ZONE_ERROR;
      }
    }
  }
  return (error_code);
} /* Convert_Geodetic_To_MGRS */
