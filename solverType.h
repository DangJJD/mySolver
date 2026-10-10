#pragma once
#include <stdio.h>

double BiSec(char *problemType, size_t length, double args[], double boundsA, double boundsB);
double FalsePos(char *problemType, size_t length, double args[], double boundsA, double boundsB);
double Newton1(char *problemType, size_t length, double args[], double boundsA, double boundsB);
double Newton2(char *problemType, size_t length, double args[], double boundsA, double boundsB);
double BracketNewton1(char *problemType, size_t length, double args[], double boundsA, double boundsB);
double BracketNewton2(char *problemType, size_t length, double args[], double boundsA, double boundsB);
double findY(char *problemType, size_t length, double args[], double x);
