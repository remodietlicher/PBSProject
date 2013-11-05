#ifndef UTIL_H
#define UTIL_H

#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#if DO_PNG_OUT==1
	#include <zlib.h>
#endif

#ifndef WRITE_NORMALIZE
#define WRITE_NORMALIZE 0
#endif // WRITE_NORMALIZE


//-----------------------------------------------------------------------------
// helper functions

template < class T > inline T ABS( T a ) { 
	return (0 < a) ? a : -a ; 
}

template < class T > inline void SWAP_POINTERS( T &a, T &b ) { 
	T temp = a;
	a = b;
	b = temp;
}

template < class T > inline void CLAMP( T &a, T b=0., T c=1.) { 
	if(a<b) { a=b; return; }
	if(a>c) { a=c; return; }
}

template < class T > inline T MIN( T a, T b) {
	return (a < b) ? a : b;
}

template < class T > inline T MAX( T a, T b) {
	return (a > b) ? a : b;
}

template < class T > inline T MAX3( T a, T b, T c) {
	T max = (a > b) ? a : b;
	max = (max > c) ? max : c;
	return max;
}

template < class T > inline double MAX3V( T vec) {
	double max = (vec[0] > vec[1]) ? vec[0] : vec[1];
	max = (max > vec[2]) ? max : vec[2];
	return max;
}

template < class T > inline double MIN3V( T vec) {
	double min = (vec[0] < vec[1]) ? vec[0] : vec[1];
	min = (min < vec[2]) ? min : vec[2];
	return min;
}

//-----------------------------------------------------------------------------
// timing functions
#if _WIN32
#include<windows.h>
#else
#include <sys/time.h>
#endif

static long int getTime()
{
  long int ret = 0;
#ifdef _WIN32
  LARGE_INTEGER liTimerFrequency;
  QueryPerformanceFrequency(&liTimerFrequency);
  LARGE_INTEGER liLastTime;
  QueryPerformanceCounter(&liLastTime);
  ret = (INT)( ((double)liLastTime.QuadPart / liTimerFrequency.QuadPart)*1000 ); // - mFirstTime;
#else
  struct timeval tv;
  struct timezone tz;
  tz.tz_minuteswest = 0;
  tz.tz_dsttime = 0;
  gettimeofday(&tv,&tz);
  ret = (tv.tv_sec*1000)+(tv.tv_usec/1000); //-mFirstTime;
#endif
  return ret;
}
static std::string timeToString(long int usecs) { 
  int ss = (int)( ((double)usecs / 1000.0) );
  int ps = (int)( ((double)usecs - (double)ss*1000.0)/10.0 );
  std::ostringstream out;
  out << ss<<".";
  if(ps<10) out <<"0";
  out <<ps<<"s";
  //out << usecs<<"u "; // debug
  return out.str();
}


//-----------------------------------------------------------------------------
// png output
#ifndef DO_PNG_OUT
#define DO_PNG_OUT 0
#endif // DO_PNG_OUT
#if DO_PNG_OUT==1
#include <png.h>
#endif

//-----------------------------------------------------------------------------
//! write png image
static int writePng(const char *fileName, unsigned char **rowsp, int w, int h, bool normalize)
{
#if DO_PNG_OUT==1
	// defaults 
	const int colortype = PNG_COLOR_TYPE_RGBA;
	const int bitdepth = 8;
	png_structp png_ptr = NULL;
	png_infop info_ptr = NULL;
	png_bytep *rows = rowsp;

	FILE *fp = NULL;
	std::ostringstream outStr;
	outStr << "open '"<< std::string(fileName) << "' for writing";
	std::string doing = outStr.str() ;
	if (!(fp = fopen(fileName, "wb"))) goto fail;

	if(!png_ptr) {
		doing = "create png write struct";
		if (!(png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL))) goto fail;
	}
	if(!info_ptr) {
		doing = "create png info struct";
		if (!(info_ptr = png_create_info_struct(png_ptr))) goto fail;
	}

	if (setjmp(png_jmpbuf(png_ptr))) goto fail;
	doing = "init IO";
	png_init_io(png_ptr, fp);
	doing = "write header";
	png_set_IHDR(png_ptr, info_ptr, w, h, bitdepth, colortype, PNG_INTERLACE_NONE,
			PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
	doing = "write info";
	png_write_info(png_ptr, info_ptr);
	doing = "write image";
	png_write_image(png_ptr, rows);
	doing = "write end";
	png_write_end(png_ptr, NULL);
	doing = "write destroy structs";
	png_destroy_write_struct(&png_ptr, &info_ptr);

	fclose( fp );
	return 0;

fail:	
	std::cerr << "writePng: could not "<<doing<<" !\n";
	if(fp) fclose( fp );
	if(png_ptr || info_ptr) png_destroy_write_struct(&png_ptr, &info_ptr);
#endif // DO_PNG_OUT==1
	return -1;
}

//-----------------------------------------------------------------------------



using namespace std;

static void WritePPM(const char *FileName, unsigned char* Color, int Width, int Height)
{
	FILE* fp;
	fopen_s(&fp,FileName, "wb");
	if (fp==NULL) { printf("PPM ERROR (WritePPM) : unable to open %s!\n",FileName); return; }
	fprintf(fp, "P6\n%d %d\n255\n", Width, Height);
	fwrite(Color,1,Width*Height*3,fp);
	fclose(fp);
} 

static void dumpNumberedImage(int counter, std::string prefix, double* field, int xRes, int yRes)
{
#if DO_PNG_OUT==1
	char buffer[256];
	sprintf(buffer,"%04i", counter);
	std::string number = std::string(buffer);

	unsigned char pngbuf[xRes*yRes*4];
	unsigned char *rows[yRes];
	double *pfield = field;
	for (int j=0; j<yRes; j++) {
		for (int i=0; i<xRes; i++) {
			double val = *pfield;
			if(val>1.) val=1.;
			if(val<0.) val=0.;
			pngbuf[(j*xRes+i)*4+0] = (unsigned char)(val*255.); 
			pngbuf[(j*xRes+i)*4+1] = (unsigned char)(val*255.); 
			pngbuf[(j*xRes+i)*4+2] = (unsigned char)(val*255.); 
			pfield++;
			pngbuf[(j*xRes+i)*4+3] = 255;
		}
		rows[j] = &pngbuf[(yRes-j-1)*xRes*4];
	}
	std::string filenamePNG = prefix + number + std::string(".png");
	writePng(filenamePNG.c_str(), rows, xRes, yRes, false);
	printf("Writing %s\n", filenamePNG.c_str());
#endif
}


#endif
