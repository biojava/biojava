#! /usr/bin/env python
#
# Script to conservatively convert initial spaces to tabs
#
# The script first guesses the number of spaces used (e.g. 2,3,4, or 8).
# It then passes this to the unexpand utility for the actual conversion.
import sys,os
import optparse
import subprocess

def space_distribution(filename):
    """Reads the specified file and creates a histogram of the number of initial
    spaces in each line of the file.

    Returns a dictionary mapping n to the number of lines with n initial spaces
    """
    spaces = {}
    with open(filename,'r') as file:
        for line in file.readlines():
            s = 0
            for c in line:
                if c == ' ':
                    s += 1
                elif c == '*':
                    s -= 1
                    break
                else:
                    break
            spaces.setdefault(s,0)
            spaces[s] += 1
    return spaces

def max_compatible(distribution, threshold=.8, possible=[8,4,3,2]):
    """Guess the most likely number of spaces per tab for a file

    distribution: a dictionary giving the number of lines with a particular
        number of initial spaces
    threshold: Fraction of lines that must be correctly spaced
    possible: list of possible tabstops to consider

    Returns the highest tabstop such that at least threshold-fraction of the
        lines contain initial spaces that are a multiple of that tabstop.
        Returns None if no possiblities meet the threshold, or if the file
        did not contain sufficient lines to determine the tabstop
        unambiguously.
    """
    for poss in sorted(possible,reverse=True):
        compat = 0
        incompat = 0
        for spaces, count in distribution.items():
            # Only assess lines with spaces
            if spaces == 0:
                continue
            if spaces%poss == 0:
                compat += count
            else:
                incompat += count
        # unambiguous case
        if compat > 0 and incompat == 0:
            return poss
        # require sufficient lines with spaces
        if threshold < 1 and (compat+incompat)*(1-threshold) < 2:
            continue
        if compat >= (compat+incompat)*threshold:
            return poss

def unexpand(filename, spaces):
    """Calls the unexpand command to convert initial spaces into tabs
    """
    args = ["unexpand","--first-only","-t",str(spaces),filename]
    print " ".join(args)
    bits = subprocess.check_output(args)
    with open(filename,'wb') as out:
        out.write(bits)


if __name__ == "__main__":
    parser = optparse.OptionParser( usage="usage: python %prog [options] filenames..." )
    parser.add_option("-t","--threshold",help="fraction of lines with compatible spacing",
        dest="threshold",default=0.8, type='float',action='store')
    parser.add_option("-x","--unexpand",help="Unexpand spaces in the input files",
        dest="unexpand",default=False, action='store_true')
    (options,args) = parser.parse_args()

    if len(args) < 1:
        parser.print_usage()
        parser.exit("Error. No files")


    for f in args:
        dist = space_distribution(f)
        # skip tabbed files
        if dist.keys() == [0]:
            continue
        c = max_compatible(dist,options.threshold)
        if c is not None:
            if options.unexpand:
                unexpand(f,c)
            else:
                print "Match\t%s\t%d\t%s"%(f,c,dist)
        else:
            print "Unclear\t%s\t?\t%s"%(f,dist)

