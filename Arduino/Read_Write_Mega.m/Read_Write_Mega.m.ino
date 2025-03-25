/*
  Project: ECG Analyzer Application
  The University of Texas at Dallas
  Joe Epperson
  Shamman Noor Shoudha
  Jason Hoff
 */
#include <SoftwareSerial.h>  

const byte bluetoothTx = 10;  // TX-O pin of bluetooth mate, Arduino D2
const byte bluetoothRx = 11;  // RX-I pin of bluetooth mate, Arduino D3

SoftwareSerial bluetooth(bluetoothTx, bluetoothRx);

void setup()
{
  Serial.begin(9600);  // Begin the serial monitor at 9600bps
  bluetooth.begin(115200);  // The Bluetooth Mate defaults to 115200bps
  bluetooth.print("$");  // Print three times individually
  bluetooth.print("$");
  bluetooth.print("$");  // Enter command mode
  delay(100);  // Short delay, wait for CMD
  bluetooth.println("U,9600,N");  // Temporarily change the baudrate to 9600, no parity
  bluetooth.begin(9600);  // Start bluetooth serial at 9600
}

void loop()
{
  if(bluetooth.available())
  {
    Serial.print((char)bluetooth.read());  
  }
  if(Serial.available())
  {
    bluetooth.print((char)Serial.read());
  }
}
