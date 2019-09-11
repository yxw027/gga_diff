// gga_diff.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "../nemagga.h"

using namespace gnssimu_lib;

void	att2C_nb(double* att, double C_nb[3][3])
{
	// attitude: roll, pitch and heading
	// attitude => C_nb
	double R = att[0], P = att[1], H = att[2];
	C_nb[0][0] = cos(H) * cos(P);
	C_nb[1][0] = sin(H) * cos(P);
	C_nb[2][0] = -sin(P);
	C_nb[0][1] = -sin(H) * cos(R) + cos(H) * sin(P) * sin(R);
	C_nb[1][1] = cos(H) * cos(R) + sin(H) * sin(P) * sin(R);
	C_nb[2][1] = cos(P) * sin(R);
	C_nb[0][2] = sin(H) * sin(R) + cos(H) * sin(P) * cos(R);
	C_nb[1][2] = -cos(H) * sin(R) + sin(H) * sin(P) * cos(R);
	C_nb[2][2] = cos(P) * cos(R);
	return;
}

void rotate_vector(double C[3][3], double* vec, int isTranspose)
{
	double temp[3] = { 0.0 };
	if (isTranspose == 0)
	{
		temp[0] = C[0][0] * vec[0] + C[0][1] * vec[1] + C[0][2] * vec[2];
		temp[1] = C[1][0] * vec[0] + C[1][1] * vec[1] + C[1][2] * vec[2];
		temp[2] = C[2][0] * vec[0] + C[2][1] * vec[1] + C[2][2] * vec[2];
	}
	else
	{
		temp[0] = C[0][0] * vec[0] + C[1][0] * vec[1] + C[2][0] * vec[2];
		temp[1] = C[0][1] * vec[0] + C[1][1] * vec[1] + C[2][1] * vec[2];
		temp[2] = C[0][2] * vec[0] + C[1][2] * vec[1] + C[2][2] * vec[2];
	}
	vec[0] = temp[0];
	vec[1] = temp[1];
	vec[2] = temp[2];
	return;
}


typedef struct
{
	int wn;
	double ws;
	double rpy[3]; /* Roll, Pitch and Yaw/Heading */
	double vxyz[3];
	double vENU[3];
	double wxyz[3];
	double fxyz[3];
	double xyz[3];
	double rmsRPY[3];
	double rmsENU[3];
	double rmsVenu[3];
	double blh[3];
	int ns;
	double hdop;
}solu1_t;

typedef struct
{
	int wn;
	double ws;
	double blh[3];
	double hdop;
	double pdop;
	int ns;
	double rmsENU[3];
	int type;
}solu2_t;

#define LAO_EST

unsigned long solu_diff(const char* fname1, const char* fname2, int opt)
{
	/* compare the soltion from NovAtel SPAN solution and NMEA GGA */
#if 0
	/* NMEA GGA */
	Week, GPSTime, Roll, Pitch, Heading, VX - ECEF, VY - ECEF, VZ - ECEF, VEast, VNorth, VUp, AngRateX, AngRateY, AngRateZ, AccBdyX, AccBdyY, AccBdyZ, X - ECEF, Y - ECEF, Z - ECEF, RollSD, PitchSD, HdngSD, SDEast, SDNorth, SDHeight, SD - VE, SD - VN, SD - VH, Latitude, Longitude, H - Ell, NS, HDOP
	(weeks, (sec), (deg), (deg), (deg), (m / s), (m / s), (m / s), (m / s), (m / s), (m / s), (deg / s), (deg / s), (deg / s), (m / s ^ 2), (m / s ^ 2), (m / s ^ 2), (m), (m), (m), (deg), (deg), (deg), (m), (m), (m), (m / s), (m / s), (m / s), (deg), (deg), (m), , (dop)
		2069, 423470, 0.391943239, 2.563304313, 60.98594284, 0, 0, 0, 0, 0, 0, 0.0041, -0.0044, 0.0023, 0.023, -0.038, -0.023, -2687728.289, -4281324.497, 3876272.875, 0.007077488, 0.007051909, 0.010994121, 0.011, 0.011, 0.021, 0.003, 0.003, 0.003, 37.66736958, -122.1197648, -19.819, 7, 1.4
#endif
	FILE * fSOL[3] = { NULL };
	char buffer[1024] = { 0 }, outfname[255] = { 0 };
	fSOL[0] = fopen(fname1, "r");
	fSOL[1] = fopen(fname2, "r");
	fSOL[1] = fopen(fname2, "r");

	memcpy(buffer, fname1, strlen(fname1));

	char* result1 = strrchr(buffer, '.');
	if (result1 != NULL) result1[0] = '\0';

	sprintf(outfname, "%s_diff.txt", buffer); 
	if (opt==0)
		fSOL[2] = fopen(outfname, "w");
	else
		fSOL[2] = fopen(outfname, "a");

	if (fSOL[0] == NULL || fSOL[1] == NULL)
	{
		if (fSOL[0] != NULL) fclose(fSOL[0]);
		if (fSOL[1] != NULL) fclose(fSOL[1]);
		return 0;
	}
	/* skip two lines */
	fgets(buffer, sizeof(buffer), fSOL[1]);
	fgets(buffer, sizeof(buffer), fSOL[1]);
	solu1_t sol[2] = { 0 };
	solu2_t sol2 = { 0 };
	double xyz_offset[] = { 0.1612,       -0.0036,        0.0859 };
	double lao[] = { 1.6066,        0.0371,       -1.0164 };
	unsigned long numofepoch = 0;
	while (!feof(fSOL[0]))
	{
		memset(&sol2, 0, sizeof(sol2));
		fgets(buffer, sizeof(buffer), fSOL[0]);
		if (sscanf(buffer, "%i,%lf,%lf,%lf,%lf,%lf,%lf,%i,%lf,%lf,%lf,%i"
			, &sol2.wn, &sol2.ws
			, sol2.blh + 0, sol2.blh + 1, sol2.blh + 2
			, &sol2.hdop, &sol2.pdop
			, &sol2.ns
			, sol2.rmsENU + 0, sol2.rmsENU + 1, sol2.rmsENU + 2
			, &sol2.type) < 12)
		{
			continue;
		}
		sol2.ws /= 1000.0;
		sol2.blh[0] *= PI / 180.0;
		sol2.blh[1] *= PI / 180.0;
		while (!feof(fSOL[1]))
		{
			if (sol[0].ws <= sol2.ws && sol[1].ws >= sol2.ws)
			{
				/* compare the solution */
				int ii = 0;
				break;
			}
			solu1_t sol1 = { 0 };
			fgets(buffer, sizeof(buffer), fSOL[1]);
			if (sscanf(buffer, "%i,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%i,%lf"
				, &sol1.wn, &sol1.ws
				, sol1.rpy + 0, sol1.rpy + 1, sol1.rpy + 2
				, sol1.vxyz + 0, sol1.vxyz + 1, sol1.vxyz + 2
				, sol1.vENU + 0, sol1.vENU + 1, sol1.vENU + 2
				, sol1.wxyz + 0, sol1.wxyz + 1, sol1.wxyz + 2
				, sol1.fxyz + 0, sol1.fxyz + 1, sol1.fxyz + 2
				, sol1.xyz + 0, sol1.xyz + 1, sol1.xyz + 2
				, sol1.rmsRPY + 0, sol1.rmsRPY + 1, sol1.rmsRPY + 2
				, sol1.rmsENU + 0, sol1.rmsENU + 1, sol1.rmsENU + 2
				, sol1.rmsVenu + 0, sol1.rmsVenu + 1, sol1.rmsVenu + 2
				, sol1.blh + 0, sol1.blh + 1, sol1.blh + 2
				, &sol1.ns, &sol1.hdop) < 34)
			{
				continue;
			}
			sol1.blh[0] *= PI / 180.0;
			sol1.blh[1] *= PI / 180.0;
			sol[0] = sol[1];
			sol[1] = sol1;
		}
		if (sol[0].ws == 0.0 && sol[1].ws >= sol2.ws)
		{
			printf("no overlap start at %10.3f\n", sol2.ws);
			break;
		}
		if (sol[0].ws <= sol2.ws && sol[1].ws >= sol2.ws)
		{
			/* interpolate */
			double dt1 = sol2.ws - sol[0].ws;
			double dt2 = sol2.ws - sol[1].ws;
			double dt = 0.0;
			double xyz_ref[3] = { 0 };
			double vxyz_ref[3] = { 0 };
			double att_ref[3] = { 0 };
			if (fabs(dt1) > fabs(dt2))
			{
				for (int i = 0; i < 3; ++i)
				{
					xyz_ref[i] = sol[1].xyz[i];
					vxyz_ref[i] = sol[1].vxyz[i];
					att_ref[i] = sol[1].rpy[i] * PI / 180.0;
				}
				dt = dt2;
			}
			else
			{
				for (int i = 0; i < 3; ++i)
				{
					xyz_ref[i] = sol[0].xyz[i];
					vxyz_ref[i] = sol[0].vxyz[i];
					att_ref[i] = sol[0].rpy[i] * PI / 180.0;
				}
				dt = dt1;
			}
			double xyz_[3] = { 0.0 };
			for (int i = 0; i < 3; ++i)
				xyz_[i] = xyz_ref[i] + vxyz_ref[i] * dt;
			double C_en[3][3] = { 0 };
			double xyz[3] = { 0 };
			blh2xyz(sol2.blh, xyz);
			blh2C_en(sol[0].blh, C_en);
			double dxyz[3] = { xyz[0] - xyz_[0]+xyz_offset[0], xyz[1] - xyz_[1] + xyz_offset[1], xyz[2] - xyz_[2] + xyz_offset[2] };
			double dned[3] = { 0 };
			double C_nb[3][3] = { 0 };
			att2C_nb(att_ref, C_nb);
			dned[0] = C_en[0][0] * dxyz[0] + C_en[1][0] * dxyz[1] + C_en[2][0] * dxyz[2];
			dned[1] = C_en[0][1] * dxyz[0] + C_en[1][1] * dxyz[1] + C_en[2][1] * dxyz[2];
			dned[2] = C_en[0][2] * dxyz[0] + C_en[1][2] * dxyz[1] + C_en[2][2] * dxyz[2];

#ifndef LAO_EST
			double db[3] = { lao[0],lao[1], lao[2] };
			rotate_vector(C_nb, db, 0);
#else
			double db[3] = { dned[0], dned[1], dned[2] };
			rotate_vector(C_nb, db, 1);
#endif

			dned[0] -= db[0];
			dned[1] -= db[1];
			dned[2] -= db[2];

			double level_b[3] = { 0 };
			double diffH = sqrt(dned[0] * dned[0] + dned[1] * dned[1]);
			double diff3 = sqrt(dned[0] * dned[0] + dned[1] * dned[1]+ dned[2] * dned[2]);

			fprintf(fSOL[2], "%10.3f,%14.9f,%14.9f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%02i,%02i,%7.3f,%7.3f,%7.3f,%10.4f,%10.4f,%10.4f,%7.3f,%7.3f,%7.3f,%10.4f,%10.4f,%10.4f\n"
				, sol2.ws, sol2.blh[0] * 180.0 / PI, sol2.blh[1] * 180.0 / PI, sol2.blh[2]
				, dned[0], dned[1], dned[2], diffH, diff3, sol2.ns, sol2.type, sol2.rmsENU[0], sol2.rmsENU[1], sol2.rmsENU[2], sol[0].vENU[1], sol[0].vENU[0],-sol[0].vENU[2]
				, db[0], db[1], db[2], att_ref[0], att_ref[1], att_ref[2]
				);
			++numofepoch;
		}
		int jj = 0;
	}
	if (fSOL[0] != NULL) fclose(fSOL[0]);
	if (fSOL[1] != NULL) fclose(fSOL[1]);
	if (fSOL[2] != NULL) fclose(fSOL[2]);
	return numofepoch;
}

int main(int argc, char* argv[])
{
	//solu_diff("D:\\aceinna\\tesla\\tesla_data\\248\\248_openRTK\\f9_190905_181418_DC-SF_sat.csv", "D:\\aceinna\\tesla\\tesla_data\\248\\reference\\NMPL18520008S_2019-09-05_22-35-53.csv", 0);
	//solu_diff("D:\\aceinna\\tesla\\tesla_data\\248\\248_openRTK\\f9_190905_181418_DC-SF_sat.csv", "D:\\aceinna\\tesla\\tesla_data\\248\\reference\\NMPL18520008S_2019-09-05_21-26-45.csv", 0);

	solu_diff("D:\\aceinna\\tesla\\tesla_data\\248\\248_openRTK\\f9_190905_181418_DC-SF_sat.csv", "D:\\aceinna\\tesla\\tesla_data\\248\\reference\\NMPL18520008S_2019-09-05_21-26-45.csv", 0);

	//solu_diff("D:\\aceinna\\tesla\\tesla_data\\248_rtklib\\f9_190905_181418_DC-SF.csv", "D:\\aceinna\\tesla\\tesla_data\\248\\reference\\NMPL18520008S_2019-09-05_21-26-45.csv", 0);
#if 0
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
#endif
	return 0;
}
