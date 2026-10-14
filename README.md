FCC-ee EFT analysis
====

Analysis for EFT in FCC-ee


Likelihood scan
====

Likelihood scan for one Wilson coefficient C, product of Poisson.
Input: 3 histograms, SM, C=+1, C=-1. Model, bin by bin:

    N(C) = SM + C * L + C^2 * Q

with L and Q extracted from C=+1 and C=-1. Data = SM (Asimov).

Create test input (input.root with histo_SM, histo_cP1, histo_cM1):

    root -l -b -q createTestInput.cxx

Run (default: input.root, histo_SM, histo_cP1, histo_cM1, C in [-1, 1], 201 points):

    root -l likelihood.cxx

With your own file and histograms:

    root -l 'likelihood.cxx("myfile.root", "histo_SM", "histo_cP1", "histo_cM1", -2, 2, 401)'

Output:

    best fit C, crossings at -2DlnL = 1 (68% CL) and 3.84 (95% CL), on screen
    canvas cc_scan:   -2 DlnL vs C
    canvas cc_inputs: the 3 input histograms

NB: if no crossing is found, enlarge the C range.
