


package org.biojava.utils.automata;

public interface PatternListener
{
    public void notifyHit(String patternID, int start, int end, boolean greedy);
}

