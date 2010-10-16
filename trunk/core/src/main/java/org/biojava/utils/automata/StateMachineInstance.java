


package org.biojava.utils.automata;


public interface StateMachineInstance
{
    /**
     * invoke transition from current state indicated by Symbol represented by symbol index.
     * @param symIdx alphabet index of the symbol encountered.
     * @return true if the symbol is valid and this state machine
     *         should continue to receive input.
     */
    public boolean transit(int symIdx);
    public StateMachineFactory parent();
    public int getStart();
}

