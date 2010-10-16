

package org.biojava.utils.automata;

/**
 * Class that produces StateMachineInstance objects.
 *
 * @author David Huen
 * @since 1.4
 */
public interface StateMachineFactory
{
    /**
     * Return a StateMachineInstance if the Symbol represented
     * by the symbol index is valid as the initial symbol of
     * the pattern. 
     * The returned StateMachineInstance will have its statepointer
     * updated to show receipt of the specified symbol. 
     * This method should not be used outside
     * of the package as it does no alphabet checks at all.
     * It should be package-private except I cannot define
     * an interface with such methods.
     *
     * @param symIdx alphabet index value for specified symbol.
     * @return an instance of StateMachineInstance if symbol
     * is valid otherwise null.
     */
    public StateMachineInstance startInstance(int symIdx, int start);
    public void setListener(PatternListener listener);
}

