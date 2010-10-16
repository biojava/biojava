

package org.biojava.utils.automata;

import java.util.Set;

import org.biojava.bio.symbol.Symbol;

public interface NfaBuilder
{
    public FiniteAutomaton getAutomaton();
    public FiniteAutomaton.Node getStart();
    public FiniteAutomaton.Node getEnd();
    public FiniteAutomaton.Node addNode(boolean isTerminal);
    public FiniteAutomaton.Transition addTransition(FiniteAutomaton.Node start, FiniteAutomaton.Node end, Symbol sym);
    public FiniteAutomaton.Transition addEpsilonTransition(FiniteAutomaton.Node start, FiniteAutomaton.Node end);
    public FiniteAutomaton.Transition addLambdaTransition(FiniteAutomaton.Node start, FiniteAutomaton.Node end);
    public FiniteAutomaton.NodeSet getNodes();
    public Set getTransitions();
    public FiniteAutomaton.NodeSet createNodeSet();
    public String toString();
}


