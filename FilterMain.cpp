#include <stdio.h>
#include "cs1300bmp.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "Filter.h"
#include <immintrin.h>
#include <omp.h>

using namespace std;

#include "rdtsc.h"

//
// Forward declare the functions
//
Filter * readFilter(string filename);
double applyFilter(Filter *filter, cs1300bmp *input, cs1300bmp *output);

int
main(int argc, char **argv)
{

  if ( argc < 2) {
    fprintf(stderr,"Usage: %s filter inputfile1 inputfile2 .... \n", argv[0]);
  }

  //
  // Convert to C++ strings to simplify manipulation
  //
  string filtername = argv[1];

  //
  // remove any ".filter" in the filtername
  //
  string filterOutputName = filtername;
  string::size_type loc = filterOutputName.find(".filter");
  if (loc != string::npos) {
    //
    // Remove the ".filter" name, which should occur on all the provided filters
    //
    filterOutputName = filtername.substr(0, loc);
  }

  Filter *filter = readFilter(filtername);

  double sum = 0.0;
  int samples = 0;

  string inputFilename; //my change start
  string outputFilename;
  struct cs1300bmp *input;
  struct cs1300bmp *output;
  int ok;
  double sample; //my change end

  for (int inNum = 2; inNum < argc; inNum++) {
    inputFilename = argv[inNum];
    outputFilename = "filtered-" + filterOutputName + "-" + inputFilename;
    input = new struct cs1300bmp;
    output = new struct cs1300bmp;
    ok = cs1300bmp_readfile( (char *) inputFilename.c_str(), input);

    if ( ok ) {
      sample = applyFilter(filter, input, output);
      sum += sample;
      samples++;
      cs1300bmp_writefile((char *) outputFilename.c_str(), output);
    }
    delete input;
    delete output;
  }

  fprintf(stdout, "Average cycles per sample is %f\n", sum / samples);

}

class Filter *
readFilter(string filename)
{
  ifstream input(filename.c_str());

  if ( ! input.bad() ) {
    int size = 0;
    input >> size;
    Filter *filter = new Filter(size);
    int div;
    input >> div;
    filter -> setDivisor(div);

    int value; //my change

    for (int i=0; i < size; i++) {
      for (int j=0; j < size; j++) {
        input >> value;
        filter -> set(i,j,value);
      }
    }
    return filter;

  } else {
    cerr << "Bad input in readFilter:" << filename << endl;
    exit(-1);
  }
}


double
applyFilter(class Filter *filter, cs1300bmp *input, cs1300bmp *output)
{

  long long cycStart, cycStop;

  cycStart = rdtscll();

  output -> width = input -> width;
  output -> height = input -> height;

  int width = input -> width - 1; //start my change
  int height = input -> height - 1;
  double filterDivisor = 1 / (float) filter -> getDivisor();
  int* filterData = filter -> getdata();
  int rval;
  int gval; 
  int bval;
  pixel temp;
  pixel pix11;
  pixel pix12;
  pixel pix13;
  pixel pix21;
  pixel pix22;
  pixel pix23;
  pixel pix31;
  pixel pix32;
  pixel pix33;

omp_set_num_threads(15);
#pragma omp parallel for private(rval, gval, bval, temp, pix11, pix12, pix13, pix21, pix22, pix23, pix31, pix32, pix33)
  for(int row = 1; row < height; row = row + 1) { //swapped row and col for loop order. unroll these loops as well
    for(int col = 1; col < width; col = col + 1) {
        rval = 0;
        gval = 0;
        bval = 0;
        pix11 = input -> color[row-1][col-1];
        pix12 = input -> color[row-1][col];
        pix13 = input -> color[row-1][col+1];
        pix21 = input -> color[row][col-1];
        pix22 = input -> color[row][col];
        pix23 = input -> color[row][col+1];
        pix31 = input -> color[row+1][col-1];
        pix32 = input -> color[row+1][col];
        pix33 = input -> color[row+1][col+1];

        rval = pix11.red * filterData[0] + pix12.red * filterData[1] + pix13.red * filterData [2] + pix21.red * filterData[3] + pix22.red * filterData[4] + pix23.red * filterData[5] + pix31.red * filterData[6] + pix32.red * filterData[7]+ pix33.red * filterData[8]; 

        gval = pix11.green * filterData[0] + pix12.green * filterData[1] + pix13.green * filterData[2] + pix21.green * filterData[3] + pix22.green * filterData[4] + pix23.green * filterData[5] + pix31.green * filterData[6] + pix32.green * filterData[7] + pix33.green * filterData[8];

        bval = pix11.blue * filterData[0] + pix12.blue * filterData[1] + pix13.blue * filterData[2] + pix21.blue * filterData[3] + pix22.blue * filterData[4] + pix23.blue * filterData[5] + pix31.blue * filterData[6] + pix32.blue * filterData[7] + pix33.blue * filterData[8];
        
        rval = 	rval * filterDivisor;
        bval = bval * filterDivisor;
        gval = gval * filterDivisor;

        if ( rval  < 0 ) {
          rval = 0;
        }
        if ( rval  > 255 ) { 
          rval = 255;
        }

        if ( gval  < 0 ) {
          gval = 0;
        }
        if ( gval  > 255 ) { 
          gval = 255;
        }

        if ( bval  < 0 ) {
          bval = 0;
        }
        if ( bval  > 255 ) { 
          bval = 255;
        }
        temp.red = rval;
        temp.green = gval;
        temp.blue = bval;
        output -> color[row][col] = temp;
        
    }
  }

  cycStop = rdtscll();
  double diff = cycStop - cycStart;
  double diffPerPixel = diff / (output -> width * output -> height);

  fprintf(stderr, "Took %f cycles to process, or %f cycles per pixel\n", diff, diff / (output -> width * output -> height));
  return diffPerPixel;
}
