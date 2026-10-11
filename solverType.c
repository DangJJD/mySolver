#include "solverType.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double findY(char *problemType, size_t length, double args[], double x);

// NOTE: Cant use pi this way unless I calculate it for each func
// static const double pi = acos(-1); // found in math.h examples
static const double pi = M_PI;
static const int max_loops = 100;
static const double tolerance = 1e-8;

double BiSec(char *problemType, size_t length, double args[], double boundsA, double boundsB) {

	double boundsC = (boundsA + boundsB) / 2;
	double yA = findY(problemType, length, args, boundsA);
	double yB = findY(problemType, length, args, boundsB);
	double yC = findY(problemType, length, args, boundsC);
	
	if (yA == 0) return boundsA;
	if (yB == 0) return boundsB;
	if (yC == 0) return boundsC;

	for (int i = 0; i < max_loops; i++) {	
		if ( yA * yC < 0 ) {
			boundsB = boundsC;
		} else if ( yC * yB < 0 ) {
			boundsA = boundsC;
		} else { 
			fprintf(stderr, "Error: BiSec: No Root: A:%.17g B:%.17g", boundsA, boundsB);
			exit(1);
		}
		boundsC = (boundsA + boundsB) / 2;
		yA = findY(problemType, length, args, boundsA);
		yB = findY(problemType, length, args, boundsB);
		yC = findY(problemType, length, args, boundsC);
		if (fabs(yC) < tolerance) return boundsC;
	}
	fprintf(stderr, "Error: Hit Max Loops");
	exit(1);
}
double FalsePos(char *problemType, size_t length, double args[],
                double boundsA, double boundsB) {
    double yA = findY(problemType, length, args, boundsA);
    double yB = findY(problemType, length, args, boundsB);

    if (yA == 0) return boundsA;
    if (yB == 0) return boundsB;

    if (yA * yB > 0) {
        fprintf(stderr, "Error: No root in the interval.\n");
        exit(1);
    }

    for (int i = 0; i < max_loops; i++) {
        double xC = (boundsA * yB - boundsB * yA) / (yB - yA);
        double yC = findY(problemType, length, args, xC);

        if (fabs(yC) < tolerance) {
            return xC;
        }

        if (yA * yC < 0) {
            boundsB = xC;
            yB = yC;
        } else {
            boundsA = xC;
            yA = yC;
        }
    }

    fprintf(stderr, "Error: Hit Max Loops\n");
    exit(1);
}
double Newton1(char *problemType, size_t length, double args[], double boundsA, double boundsB) {
	(void)problemType; (void)length; (void)args; (void)boundsA; (void)boundsB;
	return 1;
}
double Newton2(char *problemType, size_t length, double args[], double boundsA, double boundsB) {
	(void)problemType; (void)length; (void)args; (void)boundsA; (void)boundsB;
	return 1;
}
double BracketNewton1(char *problemType, size_t length, double args[], double boundsA, double boundsB) {
	(void)problemType; (void)length; (void)args; (void)boundsA; (void)boundsB;
	return 1;
}
double BracketNewton2(char *problemType, size_t length, double args[], double boundsA, double boundsB) {
	(void)problemType; (void)length; (void)args; (void)boundsA; (void)boundsB;
	return 1;
}


double horner(size_t length, double arr[], double x) {

	double s = arr[0];
	for (size_t i = 1; i < length ; i++) {
		s = s * x + arr[i];
	}

	return s;
}

double findY(char *problemType, size_t length, double args[], double x){
	double y = 0.0;
	double c = args[0];
	if (strcmp("POLY", problemType) == 0) {
		y = horner(length, args, x); // NOTE: x: [-inf, inf]

	} else if (strcmp("COS", problemType) == 0) {
		if (0 <= x && x <= pi ) {
			y = cos(x) - c;
		} else {
			fprintf(stderr, "X is Out of Bounds ERROR\n");
			fprintf(stderr, "COS: X: %.17g\n", x);
			exit(1);
		}

	} else if (strcmp("SIN", problemType) == 0) {
		if (-1 * pi / 2 <= x && x <= pi / 2 ) {
			y = sin(x) - c;
		} else {
			fprintf(stderr, "X is Out of Bounds ERROR\n");
			fprintf(stderr, "SIN: X: %.17g\n", x);
			exit(1);
		}

	} else if (strcmp("TAN", problemType) == 0) {
		if (-1 * pi / 2 < x && x < pi / 2 ) {
			y = tan(x) - c;
		} else {
			fprintf(stderr, "X is Out of Bounds ERROR\n");
			fprintf(stderr, "SIN: X: %.17g\n", x);
			exit(1);
		}

	} else if (strcmp("EXP", problemType) == 0) {
		y = exp(x); // NOTE: x: [-inf, inf]

	} else if (strcmp("LOG", problemType) == 0) {
		if ( 0 < x ) {
			y = log(x) - c;
		} else {
			fprintf(stderr, "X is Out of Bounds ERROR\n");
			fprintf(stderr, "Log: X: %.17g\n", x);
			exit(1);
		}
	}

	return y;
}
