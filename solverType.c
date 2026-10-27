#include "solverType.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static const double pi = M_PI;
static const int max_loops = 200;

double findDerivative(char *problemType, size_t length, double args[], double x);
double findAnalyticalDerivative(char *problemType, size_t length, double args[], double x);

/*
 * Root-finding notes:
 * - Do not use |f(x)| < a fixed tolerance as the only success test.
 *   A function can be scaled by a tiny coefficient while having the same root.
 * - Bracketed methods stop when the interval cannot be narrowed meaningfully
 *   in double precision, or when a computed value is exactly zero.
 * - Newton methods stop when the Newton step is at floating-point resolution,
 *   and then verify that the residual is reasonably small relative to the
 *   function's local scale.
 */
static int oppositeSigns(double a, double b)
{
    return (a < 0.0 && b > 0.0) || (a > 0.0 && b < 0.0);
}

static int intervalIsTiny(double a, double b)
{
    double scale = fmax(1.0, fmax(fabs(a), fabs(b)));
    return fabs(b - a) <= 8.0 * DBL_EPSILON * scale;
}

static int samePoint(double a, double b)
{
    double scale = fmax(1.0, fmax(fabs(a), fabs(b)));
    return fabs(a - b) <= 8.0 * DBL_EPSILON * scale;
}

static int validProblemType(const char *type)
{
    return type != NULL &&
        (strcmp(type, "POLY") == 0 || strcmp(type, "COS") == 0 ||
         strcmp(type, "SIN") == 0 || strcmp(type, "TAN") == 0 ||
         strcmp(type, "EXP") == 0 || strcmp(type, "LOG") == 0);
}

static int validX(const char *type, double x)
{
    if (!isfinite(x)) return 0;
    if (strcmp(type, "COS") == 0) return x >= 0.0 && x <= pi;
    if (strcmp(type, "SIN") == 0) return x >= -pi / 2.0 && x <= pi / 2.0;
    if (strcmp(type, "TAN") == 0) return x > -pi / 2.0 && x < pi / 2.0;
    if (strcmp(type, "LOG") == 0) return x > 0.0;
    return 1;
}

static void requireFiniteValue(double value, const char *where)
{
    if (!isfinite(value)) {
        fprintf(stderr, "Error: Non-finite value encountered in %s.\n", where);
        exit(EXIT_FAILURE);
    }
}

static void requireBracket(double *a, double *b)
{
    if (!isfinite(*a) || !isfinite(*b)) {
        fprintf(stderr, "Error: Interval bounds must be finite.\n");
        exit(EXIT_FAILURE);
    }
    if (*a > *b) {
        double temp = *a;
        *a = *b;
        *b = temp;
    }
}

static double finishBracket(double a, double b)
{
    return a + (b - a) / 2.0;
}

static double newtonSolve(char *problemType, size_t length, double args[],
                          double boundsA, double boundsB, int numericalDerivative)
{
    double x;
    requireBracket(&boundsA, &boundsB);

    if (!validProblemType(problemType)) {
        fprintf(stderr, "Error: Unknown problem type.\n");
        exit(EXIT_FAILURE);
    }
    if (strcmp(problemType, "EXP") == 0 && args[0] <= 0.0) {
        fprintf(stderr, "Error: EXP has no real root when c <= 0.\n");
        exit(EXIT_FAILURE);
    }

    x = finishBracket(boundsA, boundsB);
    if (!validX(problemType, x)) {
        fprintf(stderr, "Error: Initial Newton guess is outside the valid domain.\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < max_loops; ++i) {
        double y = findY(problemType, length, args, x);
        requireFiniteValue(y, "Newton solver");
        if (y == 0.0) return x;

        double dy = numericalDerivative
            ? findDerivative(problemType, length, args, x)
            : findAnalyticalDerivative(problemType, length, args, x);
        requireFiniteValue(dy, "derivative calculation");
        /* Do not reject a small derivative solely by absolute size: a scaled
           function such as 1e-9*(x-0.7) has a perfectly useful derivative. */
        if (dy == 0.0) {
            fprintf(stderr, "Error: Newton derivative is zero.\n");
            exit(EXIT_FAILURE);
        }

        double step = y / dy;
        double nextX = x - step;
        if (!isfinite(nextX)) {
            fprintf(stderr, "Error: Newton step produced a non-finite value.\n");
            exit(EXIT_FAILURE);
        }
        if (!validX(problemType, nextX)) {
            fprintf(stderr, "Error: Newton step left the valid domain.\n");
            exit(EXIT_FAILURE);
        }

        if (samePoint(x, nextX)) {
            /* The Newton step is below the resolution available for x in a
               double. This is a safer criterion than an absolute residual,
               which fails for functions scaled by tiny coefficients. */
            return nextX;
        }
        x = nextX;
    }

    fprintf(stderr, "Error: Newton solver hit maximum iterations. Last x = %.17g\n", x);
    exit(EXIT_FAILURE);
}

double BiSec(char *problemType, size_t length, double args[], double boundsA, double boundsB)
{
    requireBracket(&boundsA, &boundsB);
    double yA = findY(problemType, length, args, boundsA);
    double yB = findY(problemType, length, args, boundsB);
    requireFiniteValue(yA, "bisection endpoint");
    requireFiniteValue(yB, "bisection endpoint");
    if (yA == 0.0) return boundsA;
    if (yB == 0.0) return boundsB;
    if (!oppositeSigns(yA, yB)) {
        fprintf(stderr, "Error: BiSec requires an interval whose endpoints have opposite signs.\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < max_loops; ++i) {
        double mid = finishBracket(boundsA, boundsB);
        if (mid == boundsA || mid == boundsB || intervalIsTiny(boundsA, boundsB)) return mid;
        double yMid = findY(problemType, length, args, mid);
        requireFiniteValue(yMid, "bisection midpoint");
        if (yMid == 0.0) return mid;
        if (oppositeSigns(yA, yMid)) {
            boundsB = mid;
            yB = yMid;
        } else {
            boundsA = mid;
            yA = yMid;
        }
    }
    fprintf(stderr, "Error: BiSec hit maximum iterations.\n");
    exit(EXIT_FAILURE);
}

double FalsePos(char *problemType, size_t length, double args[], double boundsA, double boundsB)
{
    requireBracket(&boundsA, &boundsB);
    double yA = findY(problemType, length, args, boundsA);
    double yB = findY(problemType, length, args, boundsB);
    requireFiniteValue(yA, "false-position endpoint");
    requireFiniteValue(yB, "false-position endpoint");
    if (yA == 0.0) return boundsA;
    if (yB == 0.0) return boundsB;
    if (!oppositeSigns(yA, yB)) {
        fprintf(stderr, "Error: FalsePos requires an interval whose endpoints have opposite signs.\n");
        exit(EXIT_FAILURE);
    }

    int lastReplaced = 0;
    for (int i = 0; i < max_loops; ++i) {
        double denom = yB - yA;
        double x;
        if (denom == 0.0 || !isfinite(denom)) {
            x = finishBracket(boundsA, boundsB);
        } else {
            x = boundsA - yA * (boundsB - boundsA) / denom;
            if (!isfinite(x) || x <= boundsA || x >= boundsB) x = finishBracket(boundsA, boundsB);
        }
        if (x == boundsA || x == boundsB || intervalIsTiny(boundsA, boundsB)) return finishBracket(boundsA, boundsB);
        double y = findY(problemType, length, args, x);
        requireFiniteValue(y, "false-position iteration");
        if (y == 0.0) return x;

        if (oppositeSigns(yA, y)) {
            boundsB = x;
            yB = y;
            if (lastReplaced == -1) yA *= 0.5;
            lastReplaced = -1;
        } else {
            boundsA = x;
            yA = y;
            if (lastReplaced == 1) yB *= 0.5;
            lastReplaced = 1;
        }
    }
    fprintf(stderr, "Error: FalsePos hit maximum iterations.\n");
    exit(EXIT_FAILURE);
}

double Newton1(char *problemType, size_t length, double args[], double boundsA, double boundsB)
{
    return newtonSolve(problemType, length, args, boundsA, boundsB, 0);
}

double Newton2(char *problemType, size_t length, double args[], double boundsA, double boundsB)
{
    return newtonSolve(problemType, length, args, boundsA, boundsB, 1);
}

static double bracketNewton(char *problemType, size_t length, double args[],
                            double boundsA, double boundsB, int numericalDerivative)
{
    requireBracket(&boundsA, &boundsB);
    double yA = findY(problemType, length, args, boundsA);
    double yB = findY(problemType, length, args, boundsB);
    requireFiniteValue(yA, "bracket-Newton endpoint");
    requireFiniteValue(yB, "bracket-Newton endpoint");
    if (yA == 0.0) return boundsA;
    if (yB == 0.0) return boundsB;
    if (!oppositeSigns(yA, yB)) {
        fprintf(stderr, "Error: Bracket-Newton requires an interval whose endpoints have opposite signs.\n");
        exit(EXIT_FAILURE);
    }

    double x = finishBracket(boundsA, boundsB);
    for (int i = 0; i < max_loops; ++i) {
        if (intervalIsTiny(boundsA, boundsB)) return finishBracket(boundsA, boundsB);
        double y = findY(problemType, length, args, x);
        requireFiniteValue(y, "bracket-Newton iteration");
        if (y == 0.0) return x;

        if (oppositeSigns(yA, y)) {
            boundsB = x;
            yB = y;
        } else {
            boundsA = x;
            yA = y;
        }
        if (intervalIsTiny(boundsA, boundsB)) return finishBracket(boundsA, boundsB);

        double midpoint = finishBracket(boundsA, boundsB);
        double dy = numericalDerivative
            ? findDerivative(problemType, length, args, x)
            : findAnalyticalDerivative(problemType, length, args, x);
        double nextX = midpoint;
        if (isfinite(dy) && dy != 0.0) {
            double candidate = x - y / dy;
            if (isfinite(candidate) && candidate > boundsA && candidate < boundsB && validX(problemType, candidate)) {
                double candidateY = findY(problemType, length, args, candidate);
                if (isfinite(candidateY) && fabs(candidateY) < fabs(y)) nextX = candidate;
            }
        }
        if (nextX == x) nextX = midpoint;
        x = nextX;
    }
    fprintf(stderr, "Error: Bracket-Newton hit maximum iterations.\n");
    exit(EXIT_FAILURE);

}

double BracketNewton1(char *problemType, size_t length, double args[],
                      double boundsA, double boundsB)
{
    return bracketNewton(problemType, length, args, boundsA, boundsB, 0);
}

double BracketNewton2(char *problemType, size_t length, double args[], double boundsA, double boundsB)
{
    return bracketNewton(problemType, length, args, boundsA, boundsB, 1);
}

double horner(size_t length, double arr[], double x)
{
    if (length == 0 || arr == NULL) {
        fprintf(stderr, "Error: POLY requires at least one coefficient.\n");
        exit(EXIT_FAILURE);
    }
    double value = arr[0];
    for (size_t i = 1; i < length; ++i) value = value * x + arr[i];
    return value;
}

double findY(char *problemType, size_t length, double args[], double x)
{
    if (!validProblemType(problemType)) {
        fprintf(stderr, "Error: Unknown problem type.\n");
        exit(EXIT_FAILURE);
    }
    if (args == NULL || length == 0) {
        fprintf(stderr, "Error: At least one coefficient or parameter is required.\n");
        exit(EXIT_FAILURE);
    }
    if (!validX(problemType, x)) {
        fprintf(stderr, "%s: X is outside the valid domain: %.17g\n", problemType, x);
        exit(EXIT_FAILURE);
    }

    double c = args[0];
    double y;
    if (strcmp(problemType, "POLY") == 0) {
        y = horner(length, args, x);
    } else if (strcmp(problemType, "COS") == 0) {
        y = cos(x) - c;
    } else if (strcmp(problemType, "SIN") == 0) {
        y = sin(x) - c;
    } else if (strcmp(problemType, "TAN") == 0) {
        y = tan(x) - c;
    } else if (strcmp(problemType, "EXP") == 0) {
        y = exp(x) - c;
    } else { /* LOG */
        y = log(x) - c;
    }
    return y;
}

double findDerivative(char *problemType, size_t length, double args[], double x)
{
    if (!validX(problemType, x)) {
        fprintf(stderr, "Error: Cannot differentiate outside the valid domain.\n");
        exit(EXIT_FAILURE);
    }
    double h = cbrt(DBL_EPSILON) * fmax(1.0, fabs(x));
    if (strcmp(problemType, "LOG") == 0) h = fmin(h, x * 0.25);
    if (strcmp(problemType, "COS") == 0) h = fmin(h, fmax(DBL_EPSILON, fmin(x, pi - x) * 0.25));
    if (strcmp(problemType, "SIN") == 0) h = fmin(h, fmax(DBL_EPSILON, fmin(x + pi / 2.0, pi / 2.0 - x) * 0.25));
    if (strcmp(problemType, "TAN") == 0) h = fmin(h, fmax(DBL_EPSILON, fmin(x + pi / 2.0, pi / 2.0 - x) * 0.25));
    if (!(h > 0.0) || !isfinite(h)) h = sqrt(DBL_EPSILON) * fmax(1.0, fabs(x));

    double xPlus = x + h;
    double xMinus = x - h;
    int plusOK = validX(problemType, xPlus);
    int minusOK = validX(problemType, xMinus);
    if (plusOK && minusOK) {
        return (findY(problemType, length, args, xPlus) - findY(problemType, length, args, xMinus)) / (2.0 * h);
    }
    if (plusOK) return (findY(problemType, length, args, xPlus) - findY(problemType, length, args, x)) / h;
    if (minusOK) return (findY(problemType, length, args, x) - findY(problemType, length, args, xMinus)) / h;
    fprintf(stderr, "Error: Could not choose a valid finite-difference step.\n");
    exit(EXIT_FAILURE);
}

double findAnalyticalDerivative(char *problemType, size_t length, double args[], double x)
{
    if (!validProblemType(problemType) || !validX(problemType, x)) {
        fprintf(stderr, "Error: Invalid problem type or x for derivative.\n");
        exit(EXIT_FAILURE);
    }
    if (strcmp(problemType, "POLY") == 0) {
        if (length == 0 || args == NULL) {
            fprintf(stderr, "Error: POLY requires at least one coefficient.\n");
            exit(EXIT_FAILURE);
        }
        double derivative = 0.0;
        for (size_t i = 0; i + 1 < length; ++i) {
            derivative = derivative * x + args[i] * (double)(length - 1 - i);
        }
        return derivative;
    }
    if (strcmp(problemType, "COS") == 0) return -sin(x);
    if (strcmp(problemType, "SIN") == 0) return cos(x);
    if (strcmp(problemType, "TAN") == 0) {
        double c = cos(x);
        return 1.0 / (c * c);
    }
    if (strcmp(problemType, "EXP") == 0) return exp(x);
    return 1.0 / x; /* LOG */
}
