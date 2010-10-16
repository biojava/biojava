
package org.biojava.utils.automata;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.symbol.Symbol;

public class DfaBuilder
{
    private class DfaNode
    {
        FiniteAutomaton.NodeSet closureSet;
        FiniteAutomaton.Node node;

        private DfaNode(FiniteAutomaton.NodeSet closureSet)
        {
            //System.out.println("in DfaNode constructor");
            this.closureSet = closureSet;
            //System.out.println("initialising this.closureSet");
            node = dfa.addNode(isTerminal());
            //System.out.println("DfaNode created.");
        }

        private DfaNode(FiniteAutomaton.Node node, FiniteAutomaton.NodeSet closureSet)
        {
            this.closureSet = closureSet;
            this.node = node;
        }

        /**
         * checks whether this closure set has a terminal
         * Node in it.
         */
        private boolean isTerminal()
        {
            for (Iterator closI = closureSet.iterator(); closI.hasNext(); ) {
                FiniteAutomaton.Node currNode = (FiniteAutomaton.Node) closI.next();
                //System.out.println("isTerminal() checking: " + currNode);
                if (currNode.isTerminal())
                    return true;
            }
            return false;
        }

        public String toString()
        {
            StringBuffer output = new StringBuffer();

            output.append(node.toString());
            output.append("\n");
            output.append(closureSet.toString());
            output.append("\n");

            return output.toString();
        }
    }

    private Nfa nfa;
    private FiniteAutomaton dfa;
    private boolean converted = false;
    private Map closureSets = new HashMap();

    DfaBuilder(Nfa nfa)
        throws AutomatonException
    {
        this.nfa = nfa;

        dfa = new FiniteAutomaton("dfa_" + nfa.getName(), nfa.getAlphabet());

        // initialise DFA and the closureSets Map with the initial mapping.
        FiniteAutomaton.NodeSet initialSet = nfa.createNodeSet();
        initialSet.addNode(nfa.getStart());
        closureSets.put(initialSet, new DfaNode(dfa.getStart(), initialSet));
    }

    public FiniteAutomaton getDFA()
        throws AutomatonException
    {
        //System.out.println("getDFA called.");
        if (!converted) {
            constructSubsets();
            converted = true;
        }
        return dfa;
    }

    public void constructSubsets()
        throws AutomatonException
    {
        //System.out.println("constructSubsets() called.");
        // the initial NodeSet needs to contain the start state
        // and the lambda closure of the start state.
        FiniteAutomaton.NodeSet initialSet = nfa.createNodeSet();
        initialSet.addNodeSet(nfa.getLambdaClosure(nfa.getStart()));
        //System.out.println("adding initial node to dfa.");
        initialSet.addNode(nfa.getStart());
        
        //System.out.println("starting constructSubsets(...)");
        constructSubsets(getDfaNode(initialSet));
    }

    /**
     * Given a DFA node representing a Set of NFA nodes,
     * construct other DFA nodes that transit from it.
     */
    private void constructSubsets(DfaNode dfaNode)
        throws AutomatonException
    {
        //System.out.println("constructSubsets:\n" + dfaNode.toString());
        // retrieve the NFA nodes
        FiniteAutomaton.NodeSet closureSet = dfaNode.closureSet;

        // what Symbols have transitions from the source closure set?
        Set symbolSet = new HashSet();
        for (Iterator closI = closureSet.iterator(); closI.hasNext(); ) {
            org.biojava.utils.automata.Nfa.Node node = (org.biojava.utils.automata.Nfa.Node) closI.next();

            symbolSet.addAll(nfa.getSymbols(node));
        }

        // if there are lambda transitions from the source node set,
        // record this fact then remove LAMBDA from the symbol set.
        boolean isLambdaAffected = symbolSet.contains(Nfa.LAMBDA);
        if (isLambdaAffected) symbolSet.remove(Nfa.LAMBDA);

        // for each of the NFA nodes and each Symbol, compute
        // the next closure Sets and their associated DFA node.
        //System.out.println("constructSubsets going thru symbols for closure. " + symbolSet);
        for (Iterator symI = symbolSet.iterator(); symI.hasNext(); ) {
            Symbol currSymbol = (Symbol) symI.next();

            FiniteAutomaton.NodeSet closureForSymbol = nfa.createNodeSet();

            //System.out.println("constructSubsets going thru Nodes for Symbol.");
            for (Iterator closI = closureSet.iterator(); closI.hasNext(); ) {
                org.biojava.utils.automata.Nfa.Node node = (org.biojava.utils.automata.Nfa.Node) closI.next();

                // add closure set for current symbol for this NFA node to closureSet
                FiniteAutomaton.NodeSet currNodeSet = nfa.getClosure(node, currSymbol);
                //System.out.println("closure set for " + currSymbol.getName() + " from NFA node " + node.getID() + " is " + currNodeSet.toString());
                closureForSymbol.addNodeSet(currNodeSet);
                // add lambda closure to closureSet too.
                closureForSymbol.addNodeSet(nfa.getLambdaClosure(node));
            }

            if (!closureForSymbol.isEmpty()) {
                // check whether this NodeSet has had a DfaNode assigned to it.
                if (closureSets.containsKey(closureForSymbol)) {
                    // the new Transition ends with a known DfaNode.

                    // add transition from dfa Node to a existing Node for this closure set.
                    DfaNode dfaDestNode = getDfaNode(closureForSymbol);
                    dfa.addTransition(dfaNode.node, dfaDestNode.node, currSymbol);
                }
                else {
                    // this is a novel closure set.

                    // create a DfaNode for this closure set.
                    DfaNode dfaDestNode = getDfaNode(closureForSymbol);

                    // add transition from dfa Node to this new Node.
                    dfa.addTransition(dfaNode.node, dfaDestNode.node, currSymbol);

                    // call self recursively with the new Node.
                    constructSubsets(dfaDestNode);
                }
            }
        }
    }

    /**
     * get the DFA Node associated with the closure Set of NFA nodes.
     * If it does not exist, create it.
     */
    private DfaNode getDfaNode(FiniteAutomaton.NodeSet nfaNodes)
    {
        //System.out.println("getDfaNode called with " + nfaNodes.toString());
        DfaNode dfaNode = (DfaNode) closureSets.get(nfaNodes);

        //System.out.println("dfaNode is " + dfaNode);
        if (dfaNode == null) {
            //System.out.println("putting new DfaNode into closureSets.");
            dfaNode = new DfaNode(nfaNodes);
            closureSets.put(nfaNodes, dfaNode);
        }

        return dfaNode;
    }
}
