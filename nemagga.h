//------------------------------------------------------------------------------
// GNSS/IMU integration
// Developed by Dr. Yudan Yi
// Created on 11/21/2012
//------------------------------------------------------------------------------
#ifndef _NEMA_GGA_H_
#define _NEMA_GGA_H_
//------------------------------------------------------------------------------
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
//------------------------------------------------------------------------------
//#include "geoconst.h"
//------------------------------------------------------------------------------
namespace gnssimu_lib
{
#pragma warning (disable:4996)
	//--------------------------------------------------------------------------
	const double PI = 3.14159265358979;
	//--------------------------------------------------------------------------
	const double VEL_LIGHT = 299792458.0;				// Vacuum speed of light
	//--------------------------------------------------------------------------
	const double We_WGS84 = 7.2921151467e-5;
	const double GM_WGS84 = 3.9860050e+14;
	const double ae_WGS84 = 6378137.0;
	const double finv_WGS84 = 298.257223563;
	const double finv_Pz90 = 298.257839303;
	const double grav_WGS84 = 9.7803267714e0;
	//--------------------------------------------------------------------------
	class TGPGGA
	{
	public:
		//----------------------------------------------------------------------
		double time, blh[3], xyz[3], N, HDOP;	
		int numOfSat, solType, lineIndex;
		//----------------------------------------------------------------------
	public:
		~TGPGGA() {}
		TGPGGA()
		{
			memset(this, '\0', sizeof(TGPGGA)); 
		}
		TGPGGA(const TGPGGA& src)
		{
			memcpy(this, &src, sizeof(TGPGGA));
		}
		TGPGGA& operator= (const TGPGGA& src)
		{
			if (this!=&src)
			{
				memcpy(this, &src, sizeof(TGPGGA));
			}
			return *this;
		}
		//----------------------------------------------------------------------
		void ReSet()
		{
			memset(this, '\0', sizeof(TGPGGA)); 
		}
		//----------------------------------------------------------------------
		bool 	operator< (const TGPGGA& src) const { return time< src.time-0.001; };
		bool 	operator==(const TGPGGA& src) const { return fabs(time-src.time)<0.001; };
		//----------------------------------------------------------------------
		bool ParseGGA(const char *buffer)
		{
			ReSet();

			std::string lineReadStr(buffer), curstr;

			std::string::size_type nLoc = lineReadStr.find("$G"), nPreLoc; if (nLoc == std::string::npos) return false;
			nLoc = lineReadStr.find("GGA"); if (nLoc == std::string::npos) return false;
			//$GPGGA,201615.60,3946.9431210,N,08404.9232523,W,5,07,1.19,255.0887,M,0.0,M,0.0,0000*6A
			nPreLoc = 0;
			nLoc = lineReadStr.find(',', nPreLoc); if (nLoc == std::string::npos) return false;
			// GPGGA
			nPreLoc = nLoc+1; 
			nLoc = lineReadStr.find(',', nPreLoc); if (nLoc == std::string::npos) return false; 
			// UTC time
			curstr = lineReadStr.substr(nPreLoc, nLoc-nPreLoc); if (curstr.length()<6) return false;
			double hh = atof(curstr.substr(0,2).c_str());
			double mm = atof(curstr.substr(2,2).c_str());
			double ss = atof(curstr.substr(4).c_str());
			time = hh*3600.0+mm*60.0+ss;
			// latitude
			nPreLoc = nLoc+1; 
			nLoc = lineReadStr.find(',', nPreLoc); if (nLoc == std::string::npos) return false;
			curstr = lineReadStr.substr(nPreLoc, nLoc-nPreLoc); if (curstr.length()<4) return false; 
			double dd = atof(curstr.substr(0,2).c_str());
			       mm = atof(curstr.substr(2).c_str());
			blh[0] = (dd+mm/60.0)*PI/180.0;
			// S/N for latitude
			nPreLoc = nLoc+1; 
			nLoc = lineReadStr.find(',', nPreLoc); if (nLoc == std::string::npos) return false;
			curstr = lineReadStr.substr(nPreLoc, nLoc-nPreLoc); if (curstr.length()<1) return false;
			if (curstr[0]=='S'||curstr[0]=='s') blh[0] =-blh[0];
			// longitude
			nPreLoc = nLoc+1; 
			nLoc = lineReadStr.find(',', nPreLoc); if (nLoc == std::string::npos) return false;
			curstr = lineReadStr.substr(nPreLoc, nLoc-nPreLoc); if (curstr.length()<5) return false;
			dd = atof(curstr.substr(0,3).c_str());
			mm = atof(curstr.substr(3).c_str());
			blh[1] = (dd+mm/60.0)*PI/180.0;
			// E/W for longitude
			nPreLoc = nLoc+1; 
			nLoc = lineReadStr.find(',', nPreLoc); if (nLoc == std::string::npos) return false;
			curstr = lineReadStr.substr(nPreLoc, nLoc-nPreLoc); if (curstr.length()<1) return false;
			if (curstr[0]=='W'||curstr[0]=='w') blh[1] =-blh[1];
			// solution type
			nPreLoc = nLoc+1; 
			nLoc = lineReadStr.find(',', nPreLoc); if (nLoc == std::string::npos) return false;
			curstr = lineReadStr.substr(nPreLoc, nLoc-nPreLoc); if (curstr.length()<1) return false;
			solType = atoi(curstr.c_str());
			// number of satellite
			nPreLoc = nLoc+1; 
			nLoc = lineReadStr.find(',', nPreLoc); if (nLoc == std::string::npos) return false;
			curstr = lineReadStr.substr(nPreLoc, nLoc-nPreLoc); if (curstr.length()<1) return false;
			numOfSat = atoi(curstr.c_str());
			// HDOP
			nPreLoc = nLoc+1; 
			nLoc = lineReadStr.find(',', nPreLoc); if (nLoc == std::string::npos) return false;
			curstr = lineReadStr.substr(nPreLoc, nLoc-nPreLoc); if (curstr.length()<1) return false;
			HDOP = atof(curstr.c_str());
			// altitude
			nPreLoc = nLoc+1; 
			nLoc = lineReadStr.find(',', nPreLoc); if (nLoc == std::string::npos) return false;
			curstr = lineReadStr.substr(nPreLoc, nLoc-nPreLoc); if (curstr.length()<1) return false;
			blh[2] = atof(curstr.c_str());
			// M/m
			nPreLoc = nLoc+1; 
			nLoc = lineReadStr.find(',', nPreLoc); if (nLoc == std::string::npos) return false;
			curstr = lineReadStr.substr(nPreLoc, nLoc-nPreLoc); if (curstr.length()<1) return false;
			if (curstr[0]!='M'&&curstr[0]!='m') return false; 
			// geo height N
			nPreLoc = nLoc+1; 
			nLoc = lineReadStr.find(',', nPreLoc); if (nLoc == std::string::npos) return false;
			curstr = lineReadStr.substr(nPreLoc, nLoc-nPreLoc); if (curstr.length()<1) return false;
			N = atof(curstr.c_str());
			// M/m
			nPreLoc = nLoc+1; 
			nLoc = lineReadStr.find(',', nPreLoc); if (nLoc == std::string::npos) return false;
			curstr = lineReadStr.substr(nPreLoc, nLoc-nPreLoc); if (curstr.length()<1) return false;
			if (curstr[0]!='M'&&curstr[0]!='m') return false; 
			blh[2] += N;
			if (blh[0]==0.0&&blh[1]==0.0&&blh[2]==0.0) return false;
			return true;
		}
		//----------------------------------------------------------------------
		bool ParsePOS(const char *buffer)
		{
			ReSet();

			char key[255];
			int nLOC = 0;
			int nStart = 0;
			// 2018/05/05 04:08:51.000  -2162694.6558   4391224.2421   4075447.4248   2   6   0.6162   1.0989   1.0774  -0.5894   0.9792  -0.4501   0.00    1.4
			nStart = 0;         nLOC = 4; strncpy(key, buffer + nStart, sizeof(char)*nLOC); key[nLOC] = '\0'; int year = atoi(key);
			nStart += nLOC + 1; nLOC = 2; strncpy(key, buffer + nStart, sizeof(char)*nLOC); key[nLOC] = '\0'; int mon = atoi(key);
			nStart += nLOC + 1; nLOC = 2; strncpy(key, buffer + nStart, sizeof(char)*nLOC); key[nLOC] = '\0'; int day = atoi(key);
			nStart += nLOC + 1; nLOC = 2; strncpy(key, buffer + nStart, sizeof(char)*nLOC); key[nLOC] = '\0'; int hh = atoi(key);
			nStart += nLOC + 1; nLOC = 2; strncpy(key, buffer + nStart, sizeof(char)*nLOC); key[nLOC] = '\0'; int mm = atoi(key);
			nStart += nLOC + 1; nLOC = 6; strncpy(key, buffer + nStart, sizeof(char)*nLOC); key[nLOC] = '\0'; double ss = atof(key);
			nStart += nLOC;
			int flag = 0, ns = 0;
			double covP[6] = { 0.0 };
			double age = 0.0;
			double ratio = 0.0;
			sscanf(buffer + nStart, "%lf %lf %lf %i %i %lf %lf %lf %lf %lf %lf %lf %lf", xyz + 0, xyz + 1, xyz + 2, &flag, &ns
				, covP + 0, covP + 1, covP + 2, covP + 3, covP + 4, covP + 5
				, &age, &ratio);

			time = hh * 3600 + mm * 60 + ss;
			numOfSat = ns;
			if (flag == 1)
				solType = 4;
			else if (flag == 2)
				solType = 5;
			else
				solType = 1;

			//double time, blh[3], xyz[3], N, HDOP;
			//int numOfSat, solType, lineIndex;

			return true;
		}
		//----------------------------------------------------------------------
	};
	//--------------------------------------------------------------------------
	class TNEMAGGAReader
	{
	protected:
	public:
		std::vector<TGPGGA> vGGA_;
		~TNEMAGGAReader() {}
		TNEMAGGAReader()
		{
		}
		TNEMAGGAReader(const TNEMAGGAReader& src)
		{
			vGGA_ = src.vGGA_;
		}
		TNEMAGGAReader& operator= (const TNEMAGGAReader& src)
		{
			if (this!=&src)
			{
				vGGA_ = src.vGGA_;
			}
			return *this;
		}
		//----------------------------------------------------------------------
		void ReSet()
		{
			vGGA_.clear();
		}
		//----------------------------------------------------------------------

		void ReadGGA(const char *fileName)
		{
			FILE *fGGA = fopen(fileName, "r"); if (fGGA==NULL) { return; }
			char buffer[1024];
			TGPGGA gpGGA;
			int lineIndex = 0;
			while (!feof(fGGA))
			{
				fgets(buffer, sizeof(buffer), fGGA);

				++lineIndex;
				std::string lineReadStr(buffer);
				if (gpGGA.ParseGGA(buffer))
				{
					gpGGA.lineIndex = lineIndex;
					vGGA_.push_back(gpGGA);
				}
			}
			fclose(fGGA);
			int num_of_epoch = (int)vGGA_.size();
			int num_of_fix = 0;
			int num_of_flt = 0;
			double xyz1[3] = { 0.0 };
			double xyz2[3] = { 0.0 };
			double blh1[3] = { 0 };
			double C_en[3][3] = { 0 };
			double ned2[3] = { 0 };
			double term1[9] = { 0 };
			double term2[9] = { 0 };
			std::vector<TGPGGA>::iterator pGGA;
			for (pGGA = vGGA_.begin(); pGGA != vGGA_.end(); ++pGGA)
			{
				//blh2xyz(pGGA->blh, pGGA->xyz);
				//xyz1[0] += pGGA->xyz[0];
				//xyz1[1] += pGGA->xyz[1];
				//xyz1[2] += pGGA->xyz[2];
				if (pGGA->solType == 4)
				{
					//xyz2[0] += pGGA->xyz[0];
					//xyz2[1] += pGGA->xyz[1];
					//xyz2[2] += pGGA->xyz[2];
					++num_of_fix;
				}
				else if (pGGA->solType == 5)
					++num_of_flt;
			}
			/*
			if (num_of_epoch>0)
			{
				xyz1[0] /= num_of_epoch;
				xyz1[1] /= num_of_epoch;
				xyz1[2] /= num_of_epoch;
				xyz2blh(xyz1, blh1);
				blh2C_en(blh1, C_en);
				if (num_of_fix>0)
				{
					xyz2[0] /= num_of_fix;
					xyz2[1] /= num_of_fix;
					xyz2[2] /= num_of_fix;
					ned2[0] = C_en[0][0] * (xyz2[0] - xyz1[0]) + C_en[1][0] * (xyz2[1] - xyz1[1]) + C_en[2][0] * (xyz2[2] - xyz1[2]);
					ned2[1] = C_en[0][1] * (xyz2[0] - xyz1[0]) + C_en[1][1] * (xyz2[1] - xyz1[1]) + C_en[2][1] * (xyz2[2] - xyz1[2]);
					ned2[2] = C_en[0][2] * (xyz2[0] - xyz1[0]) + C_en[1][2] * (xyz2[1] - xyz1[1]) + C_en[2][2] * (xyz2[2] - xyz1[2]);
				}
				for (pGGA = vGGA_.begin(); pGGA != vGGA_.end(); ++pGGA)
				{
					double tempNED[3] = { 0 };
					tempNED[0] = C_en[0][0] * (pGGA->xyz[0] - xyz1[0]) + C_en[1][0] * (pGGA->xyz[1] - xyz1[1]) + C_en[2][0] * (pGGA->xyz[2] - xyz1[2]);
					tempNED[1] = C_en[0][1] * (pGGA->xyz[0] - xyz1[0]) + C_en[1][1] * (pGGA->xyz[1] - xyz1[1]) + C_en[2][1] * (pGGA->xyz[2] - xyz1[2]);
					tempNED[2] = C_en[0][2] * (pGGA->xyz[0] - xyz1[0]) + C_en[1][2] * (pGGA->xyz[1] - xyz1[1]) + C_en[2][2] * (pGGA->xyz[2] - xyz1[2]);
					term1[0] = tempNED[0] * tempNED[0];
					term1[1] = tempNED[0]*2.0;
					term1[2] = tempNED[0] * tempNED[0];

					if (pGGA->solType == 4)
					{
						xyz2[0] += pGGA->xyz[0];
						xyz2[1] += pGGA->xyz[1];
						xyz2[2] += pGGA->xyz[2];
						++num_of_fix;
					}
					else if (pGGA->solType == 5)
						++num_of_flt;
				}
			}
			*/
			int num_of_rtk = num_of_fix+num_of_flt;
			double fixrate1 = 0.0;
			double fixrate2 = 0.0;
			if (num_of_epoch > 0) fixrate1 = num_of_fix*100.0 / num_of_epoch;
			if (num_of_rtk > 0) fixrate2 = num_of_fix*100.0 / num_of_rtk;
			printf("%6i,%6i,%6i,%7.2f,%7.2f,%s\n", num_of_epoch, num_of_rtk, num_of_fix, fixrate1, fixrate2, fileName);
			return;
		}
		//----------------------------------------------------------------------
		void ReadPOS(const char *fileName)
		{
			FILE *fGGA = fopen(fileName, "r"); if (fGGA == NULL) { return; }
			char buffer[1024];
			TGPGGA gpGGA;
			int lineIndex = 0;
			double refxyz[3] = { 0.0 };
			while (!feof(fGGA))
			{
				fgets(buffer, sizeof(buffer), fGGA);
				++lineIndex;
				if (strlen(buffer) < 1) continue;
				if (buffer[0] == '%')
				{
					printf("%s", buffer);
					if (strstr(buffer, "ref pos") != NULL)
					{
						sscanf(buffer + 13, "%lf %lf %lf", refxyz + 0, refxyz + 1, refxyz + 2);
					}
					continue;
				}
				if (gpGGA.ParsePOS(buffer))
				{
					gpGGA.lineIndex = lineIndex;
					vGGA_.push_back(gpGGA);
				}
			}
			fclose(fGGA);
			int num_of_epoch = (int)vGGA_.size();
			int num_of_fix = 0;
			int num_of_flt = 0;
			double xyz1[3] = { 0.0 };
			double xyz2[3] = { 0.0 };
			double blh1[3] = { 0 };
			double C_en[3][3] = { 0 };
			double ned2[3] = { 0 };
			double term1[9] = { 0 };
			double term2[9] = { 0 };
			std::vector<TGPGGA>::iterator pGGA;
			for (pGGA = vGGA_.begin(); pGGA != vGGA_.end(); ++pGGA)
			{
				if (pGGA->solType == 4)
				{
					xyz2[0] += pGGA->xyz[0];
					xyz2[1] += pGGA->xyz[1];
					xyz2[2] += pGGA->xyz[2];
					++num_of_fix;
				}
			}
			double rmsXYZ[3] = { 0.0 };
			if (num_of_fix > 0)
			{
				xyz2[0] /= num_of_fix;
				xyz2[1] /= num_of_fix;
				xyz2[2] /= num_of_fix;
				for (pGGA = vGGA_.begin(); pGGA != vGGA_.end(); ++pGGA)
				{
					if (pGGA->solType == 4)
					{
						rmsXYZ[0] += (pGGA->xyz[0] - xyz2[0])*(pGGA->xyz[0] - xyz2[0]);
						rmsXYZ[1] += (pGGA->xyz[1] - xyz2[1])*(pGGA->xyz[1] - xyz2[1]);
						rmsXYZ[2] += (pGGA->xyz[2] - xyz2[2])*(pGGA->xyz[2] - xyz2[2]);
					}
				}
				rmsXYZ[0] = sqrt(rmsXYZ[0] / num_of_fix);
				rmsXYZ[1] = sqrt(rmsXYZ[1] / num_of_fix);
				rmsXYZ[2] = sqrt(rmsXYZ[2] / num_of_fix);
			}
			printf("%14.4f%14.4f%14.4f,%s\n", refxyz[0], refxyz[1], refxyz[2], fileName);
			printf("%14.4f%14.4f%14.4f%14.4f%14.4f%14.4f,%6i\n"
				, xyz2[0], xyz2[1], xyz2[2]
				, rmsXYZ[0], rmsXYZ[1], rmsXYZ[2], num_of_fix);
			printf("%14.4f%14.4f%14.4f\n"
				, xyz2[0] - refxyz[0], xyz2[1] - refxyz[1], xyz2[2] - refxyz[2]);
			return;
		}
		//----------------------------------------------------------------------

	};
	//--------------------------------------------------------------------------
	void	blh2xyz(const double *blh, double *xyz)
	{
		// lat, lon, ht => ecef xyz
		double a = ae_WGS84, finv = finv_WGS84;
		double f = 1.0 / finv, e2 = 2 * f - f*f;
		double lat = blh[0], lon = blh[1], ht = blh[2];
		double Rw = sqrt(1 - e2*sin(lat)*sin(lat));
		double Rn = a / Rw;
		xyz[0] = (Rn + ht)*cos(lat)*cos(lon);
		xyz[1] = (Rn + ht)*cos(lat)*sin(lon);
		xyz[2] = (Rn*(1 - e2) + ht)*sin(lat);
		return;
	}
	void	xyz2blh(const double *xyz, double *blh)
	{
		// ecef xyz => blh
		double a = ae_WGS84, finv = finv_WGS84;
		double f = 1.0 / finv, e2 = 2 * f - f*f;
		double x = xyz[0], y = xyz[1], z = xyz[2], lat, lon, ht;
		double R = sqrt(x*x + y*y + z*z);
		double ang = atan(fabs(z / sqrt(x*x + y*y)));
		if (z<0.0) ang = -ang;
		double lat1 = ang;
		double Rw = sqrt(1 - e2*sin(lat1)*sin(lat1));
		lat = atan(fabs(tan(ang)*(1 + (a*e2*sin(lat1)) / (z*Rw))));
		if (z<0.0) lat = -lat;
		while (fabs(lat - lat1)>1e-12)
		{
			lat1 = lat;
			Rw = sqrt(1 - e2*sin(lat1)*sin(lat1));
			lat = atan(fabs(tan(ang)*(1 + (a*e2*sin(lat1)) / (z*Rw))));
			if (z<0.0) lat = -lat;
		}
		if (lat>PI) lat = lat - 2.0*PI;
		if (fabs(x)<1e-12) { if (y >= 0.0) lon = PI / 2.0; else lon = 3.0*PI / 2.0; }
		else
		{
			lon = atan(fabs(y / x));
			if (x>0.0) { if (y >= 0.0) lon = lon; else lon = 2.0*PI - lon; }
			else { if (y >= 0.0) lon = PI - lon; else lon = PI + lon; }
		}
		Rw = sqrt(1 - e2*sin(lat)*sin(lat));
		double Rn = a / Rw;
		ht = R*cos(ang) / cos(lat) - Rn;
		if (lon>PI) lon = lon - 2.0*PI;
		blh[0] = lat;
		blh[1] = lon;
		blh[2] = ht;
		return;
	}
	void    blh2C_en(const double *blh, double C_en[3][3])
	{
		// blh => C_en
		double lat = blh[0], lon = blh[1];//, ht = blh[2];
		C_en[0][0] = -sin(lat)*cos(lon);
		C_en[1][0] = -sin(lat)*sin(lon);
		C_en[2][0] = cos(lat);
		C_en[0][1] = -sin(lon);
		C_en[1][1] = cos(lon);
		C_en[2][1] = 0.0;
		C_en[0][2] = -cos(lat)*cos(lon);
		C_en[1][2] = -cos(lat)*sin(lon);
		C_en[2][2] = -sin(lat);
		return;
	}
	void gga_diff(const char *fname1, const char *fname2, const char *oname)
	{
		TNEMAGGAReader gga1, gga2;
		gga1.ReadGGA(fname1);
		gga2.ReadGGA(fname2);
		std::vector<TGPGGA>::iterator pGGA1, pGGA2;
		double midblh[3] = { 0 };
		int numblh = 0;
		double midxyz[3];
		double C_en[3][3];
		double dxyz[3] = { 0 };
		double ned[3] = { 0 };
		double dned[3] = { 0 };
		for (pGGA1 = gga1.vGGA_.begin(); pGGA1 != gga1.vGGA_.end(); ++pGGA1)
		{
			blh2xyz(pGGA1->blh, pGGA1->xyz);
			midblh[0] += pGGA1->blh[0];
			midblh[1] += pGGA1->blh[1];
			midblh[2] += pGGA1->blh[2];
			++numblh;
		}
		if (numblh == 0) return;
		midblh[0] /= numblh;
		midblh[1] /= numblh;
		midblh[2] /= numblh;
		blh2xyz(midblh, midxyz);
		blh2C_en(midblh, C_en);
		for (pGGA1 = gga2.vGGA_.begin(); pGGA1 != gga2.vGGA_.end(); ++pGGA1)
			blh2xyz(pGGA1->blh, pGGA1->xyz);
		FILE *fout = fopen(oname, "w");
		int num_of_fix[2] = { 0 };
		int num_of_epoch = 0;
		for (pGGA1 = gga1.vGGA_.begin(); pGGA1 != gga1.vGGA_.end(); ++pGGA1)
		{
			pGGA2 = std::find(gga2.vGGA_.begin(), gga2.vGGA_.end(), *pGGA1);
			if (pGGA2==gga2.vGGA_.end())
			{
				continue;
			}
			dxyz[0] = pGGA2->xyz[0] - pGGA1->xyz[0];
			dxyz[1] = pGGA2->xyz[1] - pGGA1->xyz[1];
			dxyz[2] = pGGA2->xyz[2] - pGGA1->xyz[2];
			dned[0] = C_en[0][0] * dxyz[0] + C_en[1][0] * dxyz[1] + C_en[2][0] * dxyz[2];
			dned[1] = C_en[0][1] * dxyz[0] + C_en[1][1] * dxyz[1] + C_en[2][1] * dxyz[2];
			dned[2] = C_en[0][2] * dxyz[0] + C_en[1][2] * dxyz[1] + C_en[2][2] * dxyz[2];
			dxyz[0] = pGGA2->xyz[0] - midxyz[0];
			dxyz[1] = pGGA2->xyz[1] - midxyz[1];
			dxyz[2] = pGGA2->xyz[2] - midxyz[2];
			ned[0] = C_en[0][0] * dxyz[0] + C_en[1][0] * dxyz[1] + C_en[2][0] * dxyz[2];
			ned[1] = C_en[0][1] * dxyz[0] + C_en[1][1] * dxyz[1] + C_en[2][1] * dxyz[2];
			ned[2] = C_en[0][2] * dxyz[0] + C_en[1][2] * dxyz[1] + C_en[2][2] * dxyz[2];
			fprintf(fout, "%10.3f,%7.3f,%7.3f,%7.3f,%14.4f,%14.4f,%14.4f,%3i,%3i,%3i,%3i,%6i,%6i\n", pGGA1->time, dned[0], dned[1], dned[2], ned[0], ned[1], ned[2], pGGA1->numOfSat, pGGA2->numOfSat, pGGA1->solType, pGGA2->solType, pGGA1->lineIndex, pGGA2->lineIndex);
			if (pGGA1->solType == 4) ++num_of_fix[0];
			if (pGGA2->solType == 4) ++num_of_fix[1];
			++num_of_epoch;
		}
		double rate1 = 0.0;
		double rate2 = 0.0;
		if (num_of_epoch > 0) rate1 = num_of_fix[0] * 100.0 / num_of_epoch;
		if (num_of_epoch > 0) rate2 = num_of_fix[1] * 100.0 / num_of_epoch;
		printf("%6i,%6i,%6i,%7.2f,%7.2f\n", num_of_epoch, num_of_fix[0], num_of_fix[1], rate1, rate2);
		fclose(fout);
	}
	//--------------------------------------------------------------------------
#pragma warning (default:4996)
};
#endif