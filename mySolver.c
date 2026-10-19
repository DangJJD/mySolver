#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "solverType.h"

void printUsage(void);

int main(int argc, char *argv[]) {
    int max = 20;

    /*
     * Minimum command:
     * ./mySolver COS BiSec 0 0 3.14
     * For POLY, at least one coefficient is required.
     */
    if (argc < 6) {
        fprintf(stderr, "Error: Not enough arguments given.\n");
        printUsage();
        return 1;
    }

    if (argc > max) {
        fprintf(stderr, "Error: Too many arguments. Max: %d\n", max);
        return 1;
    }

    char *problemType = argv[1];
    char *solverType = argv[2];

    /* Number of values between solver name and interval bounds */
    size_t length = (size_t)(argc - 5);

    double *args = malloc(length * sizeof(double));
    if (args == NULL) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        return 1;
    }

    /* Read the interval from the last two arguments */
    double bounds[2] = {
        atof(argv[argc - 2]),
        atof(argv[argc - 1])
    };

    /* Read coefficients or equation parameter */
    for (size_t i = 0; i < length; i++) {
        args[i] = atof(argv[i + 3]);
    }

    if (strcmp(problemType, "POLY") != 0 &&
        strcmp(problemType, "COS") != 0 &&
        strcmp(problemType, "SIN") != 0 &&
        strcmp(problemType, "TAN") != 0 &&
        strcmp(problemType, "EXP") != 0 &&
        strcmp(problemType, "LOG") != 0) {
        fprintf(stderr, "Error: Unknown problem type: %s\n", problemType);
        printUsage();
        free(args);
        return 1;
    }

    double output;

    if (strcmp(solverType, "BiSec") == 0) {
        output = BiSec(problemType, length, args, bounds[0], bounds[1]);
    } else if (strcmp(solverType, "FalsePos") == 0) {
        output = FalsePos(problemType, length, args, bounds[0], bounds[1]);
    } else if (strcmp(solverType, "Newton1") == 0) {
        output = Newton1(problemType, length, args, bounds[0], bounds[1]);
    } else if (strcmp(solverType, "Newton2") == 0) {
        output = Newton2(problemType, length, args, bounds[0], bounds[1]);
    } else if (strcmp(solverType, "Bracket-Newton1") == 0) {
        output = BracketNewton1(problemType, length, args, bounds[0], bounds[1]);
    } else if (strcmp(solverType, "Bracket-Newton2") == 0) {
        output = BracketNewton2(problemType, length, args, bounds[0], bounds[1]);
    } else {
        fprintf(stderr, "Error: Unknown solver type: %s\n", solverType);
        printUsage();
        free(args);
        return 1;
    }

    printf("Solver output: %.17g\n", output);

    free(args);
    return 0;
}

void printUsage(void) {
    printf("Usage:\n");
    printf("./mySolver <Problem Type> <Solver Type> <Arguments> <Range Start> <Range End>\n");
    printf("Problem types: POLY COS SIN TAN EXP LOG\n");
    printf("Solver types: BiSec FalsePos Newton1 Newton2 "
           "Bracket-Newton1 Bracket-Newton2\n");
}
