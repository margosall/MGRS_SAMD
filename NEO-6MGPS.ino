/*********************
 *10 to GPS Module TX*
 *09 to GPS Module RX*
 *********************/

//#include <SoftwareSerial.h>

#include "mgrs.h"
#include <TinyGPS.h>
#include <Wire.h>
#include <Adafruit_GFX.h>
#include <Adafruit_SSD1306.h>

#define OLED_RESET 4
Adafruit_SSD1306 display(OLED_RESET);

//SoftwareSerial mySerial(3, 4);
TinyGPS gps;

void gpsdump(TinyGPS &gps);
void printFloat(double f, int digits = 2);
char string[15];

void setup()  
{
  display.begin(SSD1306_SWITCHCAPVCC, 0x3C);
  //display.dim(1);

  display.clearDisplay();
  display.setTextSize(2);
  display.setTextColor(WHITE);
  display.setCursor(0,0);
  display.println("Loading.");
  SerialUSB.println(string);
  display.display();
  Serial1.begin(9600);
  //SerialUSB.begin(9600);
}

void loop() // run over and over
{
  bool newdata = false;
  unsigned long start = millis();
  // Every 5 seconds we print an update
  while (millis() - start < 5000) 
  {
    if (Serial1.available()) 
    
    {
      char c = Serial1.read();
      //SerialUSB.print(c);  // uncomment to see raw GPS data
      if (gps.encode(c)) 
      {
        newdata = true;
        break;  // uncomment to print new data immediately!
      }
    }
  }
  
  if (newdata) 
  {
    gpsdump(gps);
  }
  
}

void gpsdump(TinyGPS &gps)
{
  long lat, lon;
  float flat, flon;
  unsigned long age, date, time, chars;
  int year;
  int result;
  gps.get_datetime(&date, &time, &age);
  gps.f_get_position(&flat, &flon, &age);
  result = Convert_Geodetic_To_MGRS( double(flat) * PI180,  double(flon) * PI180, 5, string);
  int hours = time / 1000000;
  hours = (hours+3) % 24;
  int minutes = (time % 1000000) / 10000;
  int seconds = (time % 10000) / 100;
  int sats = gps.satellites();
  display.clearDisplay();
  display.setCursor(0,0);
  display.print(" ");
  if(hours < 10)
  {
    display.print(0);
  }
  display.print(hours); display.print(":");
  if(minutes < 10)
  {
    display.print(0);
  }
  display.print(minutes);display.print(":");
  if(seconds < 10)
  {
    display.print(0);
  }
  display.println(seconds);
  //display.print(flat,6); display.print(", "); display.println(flon,6);
  display.print("   ");
  for(int i = 0;i < 5;i++)
  {
    display.print(string[i]);
  }
  display.println();
  display.print("   ");
  for(int i = 5;i < 10;i++)
  {
    display.print(string[i]);
  }
  display.println();
  display.print("   ");
  for(int i = 10;i < 15;i++)
  {
    display.print(string[i]);
  }
  display.println();
  //SerialUSB.print(flat,6);
  //SerialUSB.print(", ");
  //SerialUSB.print(flon,6);
  //SerialUSB.print(" ");
  //SerialUSB.println(string);
  //display.setTextSize(1);

  /*display.print("Satellites: ");*/ 
  //display.setCursor(60,24);
  //display.println(sats);

  //display.setCursor(17,9);
  //display.setTextSize(2);

  display.display();

}
