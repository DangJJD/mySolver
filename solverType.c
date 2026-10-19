#include "solverType.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double findY(char *problemType, size_t length, double args[], double x);
double findDerivative(char *problemType, size_t length, double args[], double x);
double findAnalyticalDerivative(char *problemType, size_t length, double args[], double x);


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
    double x = (boundsA + boundsB) / 2.0;

	printf("Initial Guess: %g\n", x);
        double A = findY(problemType, length, args, x);
        double dA = findAnalyticalDerivative(problemType, length, args, x);
	printf("Initial Y: %g\n", A);
	printf("Initial dY: %g\n", dA);
    for (int i = 0; i < max_loops; i++) {
        double y = findY(problemType, length, args, x);

        if (fabs(y) < tolerance) {
            return x;
        }

        double dy = findAnalyticalDerivative(problemType, length, args, x);

        if (fabs(dy) < 1e-12) {
            fprintf(stderr, "Error: Newton1 derivative is too small.\n");
            exit(1);
        }

        double nextX = x - ( y / dy );

        if (!isfinite(nextX)) {
            fprintf(stderr, "Error: Newton1 produced a non-finite value.\n");
            exit(1);
        }

        x = nextX;
    }

    printf("Last X value: %g\n", x);
    fprintf(stderr, "Error: Newton1 hit max loops.\n");
    exit(1);
}

double Newton2(char *problemType, size_t length, double args[],
               double boundsA, double boundsB) {
    double x = (boundsA + boundsB) / 2.0;

    for (int i = 0; i < max_loops; i++) {
        double y = findY(problemType, length, args, x);

        if (fabs(y) < tolerance) {
            return x;
        }

        double dy = findDerivative(problemType, length, args, x);

        if (fabs(dy) < 1e-12) {
            fprintf(stderr, "Error: Newton2 derivative is too small.\n");
            exit(1);
        }

        double nextX = x - y / dy;

        if (!isfinite(nextX)) {
            fprintf(stderr, "Error: Newton2 produced an invalid value.\n");
            exit(1);
        }

        x = nextX;
    }

    fprintf(stderr, "Error: Newton2 hit max loops.\n");
    exit(1);
}


double BracketNewton1(char *problemType, size_t length, double args[],
                      double boundsA, double boundsB) {
    double yA = findY(problemType, length, args, boundsA);
    double yB = findY(problemType, length, args, boundsB);

    if (yA == 0) return boundsA;
    if (yB == 0) return boundsB;

    if (yA * yB > 0) {
        fprintf(stderr, "Error: No root in the interval.\n");
        exit(1);
    }

    double x = (boundsA + boundsB) / 2.0;

    for (int i = 0; i < max_loops; i++) {
        double y = findY(problemType, length, args, x);

        if (fabs(y) < tolerance) return x;

        if (yA * y < 0) {
            boundsB = x;
            yB = y;
        } else {
            boundsA = x;
            yA = y;
        }

        double dy = findAnalyticalDerivative(
            problemType, length, args, x);

        double nextX = (boundsA + boundsB) / 2.0;

        if (fabs(dy) > 1e-12) {
            double candidate = x - y / dy;

            if (isfinite(candidate) &&
                candidate > boundsA && candidate < boundsB) {
                double candidateY = findY(
                    problemType, length, args, candidate);

                if (fabs(candidateY) < fabs(y)) {
                    nextX = candidate;
                }
            }
        }

        x = nextX;
    }

    fprintf(stderr, "Error: BracketNewton1 hit max loops.\n");
    exit(1);
}


double BracketNewton2(char *problemType, size_t length, double args[],
                      double boundsA, double boundsB) {
    double yA = findY(problemType, length, args, boundsA);
    double yB = findY(problemType, length, args, boundsB);

    if (yA == 0) return boundsA;
    if (yB == 0) return boundsB;

    if (yA * yB > 0) {
        fprintf(stderr, "Error: No root in the interval.\n");
        exit(1);
    }

    double x = (boundsA + boundsB) / 2.0;

    for (int i = 0; i < max_loops; i++) {
        double y = findY(problemType, length, args, x);

        if (fabs(y) < tolerance) return x;

        if (yA * y < 0) {
            boundsB = x;
            yB = y;
        } else {
            boundsA = x;
            yA = y;
        }

        double dy = findDerivative(problemType, length, args, x);
        double nextX = (boundsA + boundsB) / 2.0;

        if (fabs(dy) > 1e-12) {
            double candidate = x - y / dy;

            if (isfinite(candidate) &&
                candidate > boundsA && candidate < boundsB) {
                double candidateY = findY(
                    problemType, length, args, candidate);

                if (fabs(candidateY) < fabs(y)) {
                    nextX = candidate;
                }
            }
        }

        x = nextX;
    }

    fprintf(stderr, "Error: BracketNewton2 hit max loops.\n");
    exit(1);
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
    y = exp(x) - c;
  
} else if (strcmp("LOG", problemType) == 0) {
    if (0 < x) {
        y = log(x) - c;
    } else {
        fprintf(stderr, "X is Out of Bounds ERROR\n");
        fprintf(stderr, "Log: X: %.17g\n", x);
        exit(1);
    }
}


	
		
	return y;
}

double findDerivative(char *problemType, size_t length, double args[], double x) {
    double h = 1e-6 * (1.0 + fabs(x));
    double yPlus = findY(problemType, length, args, x + h);
    double yMinus = findY(problemType, length, args, x - h);

    return (yPlus - yMinus) / (2.0 * h);
}


double findAnalyticalDerivative(char *problemType, size_t length, double args[], double x) {
    if (strcmp(problemType, "POLY") == 0) {
        double derivative = 0.0;

        for (size_t i = 0; i < length - 1; i++) {
            derivative = derivative * x + args[i] * (double)(length - 1 - i);
        }

        return derivative;
    } else if (strcmp(problemType, "COS") == 0) {
        return -sin(x);
    } else if (strcmp(problemType, "SIN") == 0) {
        return cos(x);
    } else if (strcmp(problemType, "TAN") == 0) {
        double c = cos(x);
        return 1.0 / (c * c);
    } else if (strcmp(problemType, "EXP") == 0) {
        return exp(x);
    } else if (strcmp(problemType, "LOG") == 0) {
        return 1.0 / x;
    }

    fprintf(stderr, "Error: Unknown problem type in derivative.\n");
    exit(1);
}

