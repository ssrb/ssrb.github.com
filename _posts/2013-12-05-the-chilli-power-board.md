---
layout: post
title: "The Chilli Power Board"
description: ""
category: "Home automation"
tags: [arduino, electronics]
assets: "assets/ChilliPowerBoard"
---
{% include JB/setup %}

## Hola, me dicen Dedoverde
I'm starting growing [Bhut Jolokia](http://en.wikipedia.org/wiki/Bhut_Jolokia) and other peppers out of seeds.
I would like to automate, monitor and document as much as possible the process so that I decided to build an overkill programable
power board out of stuff I ordered some time ago from [DX](http://dx.com/) as well as some scrap parts.
This includes:

* a [4x relay PCB](http://dx.com/p/4-channel-relay-module-extension-board-for-arduino-avr-arm-51-141651): ~8USD;
* a [tiny bluetooth <=> serial module](http://dx.com/p/jy-mcu-arduino-bluetooth-wireless-serial-port-module-104299): ~8USD;
* a [tiny ATmega328P based prototyping board](http://dx.com/p/diy-atmega328p-16mhz-electric-block-module-blue-172858) a.k.a. [Arduino pro mini](http://arduino.cc/en/Main/ArduinoBoardProMini): ~8USD;
* a [PCF8583 I2C real time clock](http://www.nxp.com/documents/data_sheet/PCF8583.pdf): free sample from NXP;
* resitors, crystal, "OR-ing" diodes, buzzer, &c: scrap + [Jaycar](http://www.jaycar.co.nz/), less than 10NZD;
* a power board with lots of space inside: free.

For this first milestone, I only want to be able to schedule events like turn heater, lighting or irrigation circuits on & off independantly.
Later I would like to take advantage of the remaining inputs of the embeded `MCU` to read temperature and humidity 
sensors like [this one](http://www.jaycar.co.nz/productView.asp?ID=XC4246) for example.

<!-- more -->

## Bitlash
The power board [8-bit MCU](http://www.atmel.com/Images/doc8161.pdf) will run an enhanced version of [Bitlash](http://bitlash.net/).
Enhancements include:

* a`time` user function allowing to set/get the date/time from a real time clock;
* a Unix like `cron` service relying on the Bitlash task scheduler: the cron service will check every seconds if cronjobs must be run;
* a few user functions used to interact with the cron service: `addcronjob`, `delcronjob`, `erasecrontab` & `lscrontab`;
* a basic watchdog mechanism which will warn me whem something's wrong (too many consecutive faulty cronjobs): the watchdog will play the Super Mario tune through the piezo buzzer every minutes in such case.

The code is available on github [here](https://github.com/ssrb/ChilliPowerBoard/blob/master/powerboard.ino).

In the following examples, I'm connected to my power board using `screen`:
{% highlight bash %}
screen /dev/ChilliPowerBoard 57600
{% endhighlight %}

### Switching power plugs

The `MCU` digital ouputs used to control the power plug relays go from `d2` (plug 1) to `d5` (plug 4).
For example, in order to switch the third plug on, one would type in the power board shell:
{% highlight bash %}
> d4=1
{% endhighlight %}

In order to toggle the second one:
{% highlight bash %}
> d3=!d3
{% endhighlight %}

This one could be the building block of an Annoy-a-tron:
{% highlight bash %}
> d2=random(2)
{% endhighlight %}

### Crontab

The crontab relies on Bitlash task scheduler. It runs once per second:
{% highlight bash %}
> ps
0: cronsvc
{% endhighlight %}

Toggle the plug #1 every 10 seconds:
{% highlight bash %}
> addcronjob(-10,-1,-1,-1,-1,-1,-1,"d2=!d2",0)
> lscrontab
JOB|   s    m    h  dow  dom    M    Y | CMD
 0 |*/10    *    *    *    *    *    * | d2=!d2
{% endhighlight %}

Delete cronjob 0:
{% highlight bash %}
> lscrontab
JOB|   s    m    h  dow  dom    M    Y | CMD
 0 |*/10    *    *    *    *    *    * | d2=!d2
> delcronjob(0)
> lscrontab
JOB|   s    m    h  dow  dom    M    Y | CMD
{% endhighlight %}

From monday to friday, Switch on plug #2 from 8am to 11:30am:
{% highlight bash %}
> addcronjob(0,0,8,0b01111100,-1,-1,-1,"d3=1",0)
> addcronjob(0,30,11,0b01111100,-1,-1,-1,"d3=0",1)
> lscrontab
JOB|   s    m    h  dow  dom    M    Y | CMD
 0 |   0    0    8   7C    *    *    * | d3=1
 1 |   0   30   11   7C    *    *    * | d3=0
{% endhighlight %}

Erase the crontab:
{% highlight bash %}
> erasecrontab
{% endhighlight %}

### EEPROM

All the cronjobs are stored in the MCU's' `EEPROM` so that the number of job is limited.
Since Bitlash is also using the EEPROM to store functions, we need to tell it at compile time how much of
it can be used (see the "Reserving EEPROM for Other Applications" section in the [Bitlash User's Guide](http://bitlash.net/bitlash-users-guide.pdf)).

Dump of the EEPROM:
{% highlight bash %}
> peep

E000:  cron svc$   cro nd;   $sta rtcr  ond$  cro  nd;  run   cron svc,  1000 ; $p  lug$ $...
E040:  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....
E080:  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....
E0C0:  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....
E100:  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....
E140:  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....
E180:  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....
E1C0:  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....
E200:  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....
E240:  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....
E280:  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....
E2C0:  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....
E300:  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....
E340:  .... ....  .... ....  .... ....  .... ....  .... ....  .... ....  .... ...�  .... ..d2
E380:  =!d2 $$$$  ���� ����  ���� ����  ���� ����  ���� ����  ���� ����  ���� ����  ���� ����
E3C0:  ���� ����  ���� ����  ���� ����  ���� ����  ���� ����  ���� ����  ���� ����  ���� ���.
{% endhighlight %}

## Bluetooth setup (Pimp my HC06)

First thing I'm going to do, is setup the BT module.

According to the product [datasheet]({{ site.url }}/{{ page.assets }}/HC-06_datasheet_201104_revised.pdf) (chapter 9 "AT command set"), one can change the name, the baudrate as well as the pin of the module sending
`AT` commands on its serial interface.

In order to do so, I will use an Arduino nano connected to my computer and loaded with a serial pipe program.
Then, I will just send the correct `AT` commands from a serial terminal.
I'm using the "serial monitor" built into the Arduino IDE but one could use gtkterm or even screen.

The pipe program:
{% highlight c linenos %}
#include <SoftwareSerial.h>

SoftwareSerial btSerialAdapter(2, 3); // RX, TX

void setup() { 
  Serial.begin(9600);
  while (!Serial);
  btSerialAdapter.begin(9600);
} 

void loop() { 
  if(Serial.available()) {
    btSerialAdapter.write(Serial.read());
  }
  if(btSerialAdapter.available()) {
    Serial.write(btSerialAdapter.read());
  }
} 
{% endhighlight %}

and the circuit:

![JY-MCU]({{ site.url }}/{{ page.assets }}/JY-MCU_linvor.png)

We're ready to send the `AT` commands:

<div class="galleria">
  <a href="{{ site.url }}/{{ page.assets }}/AT+VERSION.png">
   <img src="{{ site.url }}/{{ page.assets }}/AT+VERSION.png" 
   data-title="AT+VERSION" 
   data-description="This command displays the product version: linvor1.5"
   data-big="{{ site.url }}/{{ page.assets }}/AT+VERSION.png"/>
 </a>
 <a href="{{ site.url }}/{{ page.assets }}/AT+NAME.png">
   <img src="{{ site.url }}/{{ page.assets }}/AT+NAME.png" 
   data-title="AT+NAME"   	
   data-description="Set the name of the product to ChilliPowerBoard. Default was linvor. The BT manager reflects that change. It will be handy to distinguish among multiple power boards."
   data-big="{{ site.url }}/{{ page.assets }}/AT+NAME.png"/>
 </a>
 <a href="{{ site.url }}/{{ page.assets }}/AT+BAUD.png">
   <img src="{{ site.url }}/{{ page.assets }}/AT+BAUD.png" 
   data-title="AT+BAUD" 
   data-description="Set the UART to 57600Bd. Default was 9600Bd."
   data-big="{{ site.url }}/{{ page.assets }}/AT+BAUD.png"/>
 </a>
</div>

### Linux setup

Let's go through the few steps needed to easily talk to the power board from Linux.

* First thing is to retrieve the device address:
{% highlight bash linenos %}
$ hcitool scan
Scanning ...
  00:12:10:23:02:31 ChilliPowerBoard
{% endhighlight %}

* Then we can add some stuff in`/etc/bluetooth/rfcomm.conf` to connect auto-magically at startup:
{% highlight bash linenos %}
#
# RFCOMM configuration file.
#

rfcomm0 {
        # Automatically bind the device at startup
        bind yes;

        # Bluetooth address of the device
        device 00:12:10:23:02:31;

        # RFCOMM channel for the connection
        channel 1;

        # Description of the connection
        comment "Chilli Power Board";
}
{% endhighlight %}

* Then we setup permission on the `/dev/rfcomm0` file and tell modem-manager not to probe the power board.
In `/etc/udev/rules.d/100-ChilliPowerBoard.rules` (create this one):
{% highlight bash %}
KERNEL=="rfcomm0", MODE="0666", ENV{ID_MM_DEVICE_IGNORE}="1"
{% endhighlight %}

* Add ourselves to the `dialup` group (or whatever group is used for rfcomm on your distrib). In `/etc/group`:
{% highlight bash %}
dialout:x:20:sb
{% endhighlight %}

* Optionally create a user friendly symlink:
{% highlight bash %}
$ sudo ln -s /dev/rfcomm0 /dev/ChilliPowerBoard
{% endhighlight %}

* Manualy connect to the power board ?
{% highlight bash linenos %}
$ sudo rfcomm connect 0
Connected /dev/rfcomm0 to 00:12:10:23:02:31 on channel 1
Press CTRL-C for hangup
{% endhighlight %}

* Disconnect ...
{% highlight bash %}
$ sudo rfcomm release 0
{% endhighlight %}

* Talk to the board: 
{% highlight bash %}
$ screen /dev/ChilliPowerBoard 57600
{% endhighlight %}

* Optionally create a user friendly alias:
{% highlight bash %}
$ alias ChilliPowerBoard="screen /dev/ChilliPowerBoard 57600"
{% endhighlight %}

That's all ! The BT module is now ready and I will keep it on a corner of my desk until it's time to put everything together.
The next stage is to hack a little RTC circuit.

## Real Time Clock
Since I'm willing to schedule events at a given date and time on my power board, it needs to track time pretty accurately.
### Basic circuit
I had a spare `PCF8583` I2C [real time clock](http://en.wikipedia.org/wiki/Real-time_clock) I decided to use:

![PCF8583P]({{ site.url }}/{{ page.assets }}/PCF8583P.png)

All it needs is a 32.768Khz crystal (`\(2^{15}\)` cycles per seconds).
The basic circuit looks like this:

![RTC circuit on breadboard]({{ site.url }}/{{ page.assets }}/real_time_clock_basic_bboard.png)

According to [another datasheet](http://www.nxp.com/documents/user_manual/UM10301.pdf), page 10 of 52, the theoretical capacitance required at `OSCI` is 18pF. The corresponding capacitor should be connected to `\(V_{DD}\)` (page 11 of 52). I added one to the circuit in the end.
I choosed not to use a trimpot capacitor `\(C_{trim}\)` (5 to 25pF as per datasheet) since I can easily compensate for the expected 5 minutes error a year
by resetting the RTC via bluetooth, based on [NTP](http://en.wikipedia.org/wiki/Network_Time_Protocol) time for example (see page 24).

Also, page 12 of 37, it is said that
> If the alarm enable bit of the control and status register is reset (logic 0), a 1 Hz signal is observed on the interrupt pin `\(\overline{INT}\)`.

and page 13:
> In the clock mode, if the alarm enable is not activated (alarm enable bit of the control and 
> status register is logic 0), the interrupt output toggles at 1 Hz with a 50 % duty cycle (may 
> be used for calibration).

Here is what I observe (looks good to me):
 
![RTC 1Hz signal]({{ site.url }}/{{ page.assets }}/pcf8583_1hz_ref.bmp)

### Power backup

`RTC` circuits often come with a backup power. For example a 3V coin battery.
A simple way to implement this is a diode based "OR-ing" circuit.

Here is a little schematic and an attempt at explaining how it works (I hope you can read my handwriting):
![OR-ing]({{ site.url }}/{{ page.assets }}/oring_diodes.jpg)

The circuit upgraded with a backup power:
![RTC circuit with backup power on breadboard]({{ site.url }}/{{ page.assets }}/real_time_clock_backup_power_bboard.png)

I only had `1N4150` diodes available. Best practice for this circuit is to use [Schottky diode](http://en.wikipedia.org/wiki/Schottky_diode) in order to minimise the forward voltage drop.

I measured the current with and without primary power supplied and checked against the datasheet (page 21).
The measured current is higher than what's expected (4.3 vs ~3 micro amps @ 3V). I wonder why.

<div class="galleria">
  <img src="{{ site.url }}/{{ page.assets }}/real_time_clock_no_backup_current.jpg"
    data-title="With primary power supplied" 
    data-description="My guess is that I'm measuring the diode reverse current: 0.1 micro amp"
    data-big="{{ site.url }}/{{ page.assets }}/real_time_clock_no_backup_current.jpg"/>
  <img src="{{ site.url }}/{{ page.assets }}/real_time_clock_backup_current.jpg"
    data-title="Without primary power supplied" 
    data-description="The current supplied by the battery is 4.3 micro amps"
    data-big="{{ site.url }}/{{ page.assets }}/real_time_clock_backup_current.jpg"/>
  <img src="{{ site.url }}/{{ page.assets }}/PC8583PSupplyCurrentClockMode.png"
    data-title="Typical supply current in clock mode" 
    data-description="The measure is higher than what's expected on the datasheet: ~3 micro amps @ 3V"
    data-big="{{ site.url }}/{{ page.assets }}/PC8583PSupplyCurrentClockMode.png"/>
</div>

Since the circuit is simple, I decided to transfer it on a stripboard:

<div class="galleria">
  <img src="{{ site.url }}/{{ page.assets }}/real_time_clock_final1.jpg"/>
  <img src="{{ site.url }}/{{ page.assets }}/real_time_clock_final2.jpg"/>
  <img src="{{ site.url }}/{{ page.assets }}/real_time_clock_final3.jpg"/>
</div>

This wasn't a great idea. The clock is either running too fast or too slow. Just like [this guy](http://electronics.stackexchange.com/questions/79841/how-to-improve-i2c-rtc-accuracy).
Hopefuly the `UM10301` datasheet, section "14. PCB layout guidelines" will allow me
to create a proper PCB. But that will be another time.

Also it's worth noting that a [DX RTC module](http://dx.com/p/i2c-rtc-ds1307-24c32-real-time-clock-module-for-arduino-blue-149493)
costs about 3 yankee dollars.

That's all for the RTC hardware. Let's talk quickly about the software.

### Software
On the software side, a PCF8583 Arduino library is available on [github](https://github.com/edebill/PCF8583):


* it's incomplete because I found out that a dude called jiki974 added [support for the daily alarm](http://forum.arduino.cc/index.php/topic,33217.0.html) but these changes never made it to the repo;
* it's somehow buggy since it ignores the day of the week and doesn't really cope well with leap years.

Anyway, I submitted a [pull request](https://github.com/edebill/PCF8583/pull/1) to cope with these issues.

Interacting with the component is done via read and write operations from and to its `CMOS RAM`.
The map of the `RAM` is at available page 8 of the datasheet.

Ultimately this is what a `time` command using the PCF8583 library looks like in Bitlash:
{% highlight c linenos %}
numvar time() {
    if(getarg(0) == 0) {
      static const char* DayStrings[] = {
        "Sunday", "Monday", "Tuesday", 
        "Wednesday", "Thursday", "Friday", 
        "Saturday"};
      theRTClock.get_time();
      char t[50];
      sprintf_P(t, PSTR(" %04d/%02d/%02d %02d:%02d:%02d"),
        theRTClock.year, theRTClock.month, theRTClock.day,
        theRTClock.hour, theRTClock.minute, theRTClock.second);
      sp(DayStrings[theRTClock.get_day_of_week()]); sp(t); speol();
    } else if(getarg(0) == 6) {
       theRTClock.year= getarg(1);
       theRTClock.month = getarg(2);
       theRTClock.day = getarg(3);
       theRTClock.hour = getarg(4);
       theRTClock.minute = getarg(5);
       theRTClock.second = getarg(6);    
       theRTClock.set_time();
    }
   return 0;
}
{% endhighlight %}

## Flashing the MCU

I'm programming the MCU using a USBtinyISP clone, bypassing the bootloader. I find this method really reliable and 
it allows me to play with the fuse bits as well.

It looks like this:

![Programming the Mini using ISP]({{ site.url }}/{{ page.assets }}/mini_isp_prog.png)

On my Linux distribution, I setup credentials for the USBtinyISP by creating a `/etc/udev/rules.d/103-USBtinyISP.rules` file containing:
{% highlight bash %}
SUBSYSTEM=="usb", ATTR{idVendor}=="1781", ATTR{idProduct}=="0c9f", GROUP="plugdev", MODE="0666"
{% endhighlight %}

I must tell the Arduino IDE I'm using an ISP programmer: `File => Upload Using Programmer` or `Ctrl+Shift+U`, 
as well as select the proper programmer in `Tools => Programmer`. 
I could also use `avrdude` directly: 
{% highlight bash %}
$ avrdude -v -c usbtiny -p m328p -U flash:w:ChilliPowerBoard.hex
{% endhighlight %}

## Putting everything together

* Software: check !
* BT <=> serial pimping: check !
* DIY RTC: check ! (kinda)
* Flashed MCU: check !

We're good to go. Here is the plan:

![Fritzing breadboard]({{ site.url }}/{{ page.assets }}/fritzingchillipowerboard.png)

Its implementation looks like this:

![Everything together]({{ site.url }}/{{ page.assets }}/everything_together.png)

Now comes the real challenge: fit everything in the power board ...

<div class="galleria">
  <img src="{{ site.url }}/{{ page.assets }}/power_board_1.jpg"/>
  <img src="{{ site.url }}/{{ page.assets }}/power_board_3.jpg"/>
  <img src="{{ site.url }}/{{ page.assets }}/power_board_4.jpg"/>
  <img src="{{ site.url }}/{{ page.assets }}/power_board_5.jpg"/>
</div>

### Conclusion
It's a bit scarry and looks like a dangerous protoype. 
It works just fine.
In the near futur I plan to design a proper RTC circuit, hack an Android app to control
the board and ultimately read temperature and humidity sensors. That's all folks.

<script type="text/javascript" src="{{ site.url }}/rungalleria.js"></script>

