



package org.biojava.utils.automata;

import org.biojava.bio.symbol.IllegalAlphabetException;

/**
 * Interface for classes that make StateMachineFactory instances
 * that in turn make StateMachineInstances that actually do the parsing.
 * Each instance of StateMachineToolkit is created to handle one
 * FiniteAlphabet only.
 *
 * @author David Huen
 * @since 1.4.
 */
interface StateMachineToolkit
{
    /**
     * Returns a StateMachineFactory that produces instances
     * that implement a StateMachine equivalent to the
     * specified FiniteAutomaton.
     * @param factoryName label for this factory.
     * @param fa The FiniteAutomaton that describes the StateMachine.
     */
    public StateMachineFactory getFactory(String factoryName, FiniteAutomaton fa)
        throws IllegalAlphabetException;
}

