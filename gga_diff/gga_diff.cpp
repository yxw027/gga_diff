// gga_diff.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "../nemagga.h"

int main(int argc, char* argv[])
{
	gnssimu_lib::TNEMAGGAReader gga;
	if (argc >= 1)
		gga.ReadPOS(argv[1]);// "C:\\rtklib\\process\\qxwz\\82101250.pos");
#if 0
							 //if (argc < 4)
							 //		return;
							 //gnssimu_lib::gga_diff(argv[1], argv[2], argv[3]);
#else
							 //gnssimu_lib::gga_diff("D:\\RTK20170713\\Projects\\Prj_RTKOffLine\\rtkpost.sol", "D:\\RTK20170713\\0524_fix\\rtkpost.sol", "aa.txt");
	gnssimu_lib::gga_diff("D:\\share\\20170614\\BX306_05.nmea", "D:\\share\\20170614\\PP\\rtkpost.sol", "D:\\share\\20170614\\aa.csv");
#endif
	return 0;
}
