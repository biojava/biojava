


package org.biojava.utils.automata;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import org.biojava.bio.BioError;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;

/**
 * Class for modelling non-deterministic finite automata.
 * <p>
 * This implementation has epsilon and lambda transitions.
 * Both transitions are silent but the former is intended
 * to be optimised away while the latter must be retained
 * during optimisation.  This is necessary to implement
 * limited closure for the REs that one may want to build
 * with this NFA.
 *
 * @author David Huen
 * @since 1.4
 */
public class Nfa
    extends FiniteAutomaton
    implements NfaBuilder
{
    // Used to indicate a silent transition that can be munged and optimised away.
    static Symbol EPSILON = null;
    // Used to indicate a silent transition that must be preserved during munging.
    static Symbol LAMBDA = AlphabetManager.createSymbol("lambda");

    public Nfa(String name, FiniteAlphabet alfa)
    {
        super(name, alfa);
    }

    protected int alphaIndex(Symbol sym)
        throws IllegalSymbolException
    {
        if (sym == LAMBDA) return 998;
        else
            return super.alphaIndex(sym);
    }

    public boolean containsNode(Node node)
    {
        return nodes.contains(node);
    }

//    public void addNode() { super.addNode(); }

    /**
     * Add a silent optimisable transition to instance.
     */
    public Transition addEpsilonTransition(Node start, Node end)
    {
        return addTransition(start, end, EPSILON);
    }

    /**
     * Add a silent persistent transition to instance.
     */
    public Transition addLambdaTransition(Node start, Node end)
    {
        return addTransition(start, end, LAMBDA);
    }

    /**
     * merge all nodes linked by emission-less transitions.
     */
    void doEpsilonClosure()
    {
        // when accumulating closure set, ensure that it
        // start and end are in there, they become the
        // replacement.
        Set closure = new HashSet();

        boolean foundEpsilonTransitions;
        do {
            closure.clear();
            foundEpsilonTransitions = false;
            for (Iterator transI = transitions.iterator(); transI.hasNext(); ) {
                Transition currTransition = (Transition) transI.next();
    
                if (currTransition.sym == EPSILON) {
                    foundEpsilonTransitions = true;
                    if (closure.isEmpty()) {
                        // start a new closure
                        closure.add(currTransition.source);
                        closure.add(currTransition.dest);
                    }
                    else {
                        // if this transition is linked with any of those I
                        // have previously encountered, coalesce them.
                        if ((closure.contains(currTransition.source)) ||
                            (closure.contains(currTransition.dest))) {
                                closure.add(currTransition.source);
                                closure.add(currTransition.dest);
                        }
                    }
                }
            }

            // found a closure
            if (foundEpsilonTransitions) {
                // specify the Node that will act for rest
                // in closure set.
                boolean containsStart = closure.contains(start);
                boolean containsEnd = closure.contains(end);
                if (containsStart && containsEnd) {
                    throw new BioError("The epsilon transitions span entire model, you fool!");
                }

                Node vicar = null;
                if (containsStart) 
                    vicar = start;
                else if (containsEnd)
                    vicar = end;
                else
                    // silly way to have to retrieve an entry from a set....
                    vicar = (Node) closure.iterator().next();
                
                replaceNode(closure, vicar);
            }
        }
        while (foundEpsilonTransitions);
    }

    /**
     * Retrieve all Nodes reachable from the specified node by
     * emissionless lambda transitions.
     */
    // FIXME: take into consideration cycles of lambda transitions!
    NodeSet getLambdaClosure(Node node)
        throws AutomatonException
    {
        return _getLambdaClosure(node, createNodeSet());
    }

    private NodeSet _getLambdaClosure(Node node, NodeSet visitedNodes)
        throws AutomatonException
    {
        // I've visited this one too!
        visitedNodes.addNode(node);

        NodeSet closureSet = createNodeSet();

        NodeSet thisClosure = getClosure(node, LAMBDA);
        closureSet.addNodeSet(thisClosure);

        for (Iterator closI = thisClosure.iterator(); closI.hasNext(); ) {
            Node currNode = (Node) closI.next();
            if (visitedNodes.contains(currNode)) continue;
            closureSet.addNodeSet(_getLambdaClosure(currNode, visitedNodes));
        }

        return closureSet;
    }

    /**
     * goes thru data structures replacing every instance
     * of old Node with new Node.  Duplicate entries that
     * arise from the process are removed too.
     * <p>
     * The Nfa version of this method also strips
     * epsilon self-transitions.
     */
    private void replaceNode(Set oldNodes, Node newNode)
    {
        //System.out.println("oldNodes: " + oldNodes);
        //System.out.println("newNode:  " + newNode);
        // prepare to replace entire contents of transitions
        Transition [] transitionArray = new Transition[transitions.size()];

        // loop thru' all transitions replacing the oldNodes
        int j = 0;
        for (Iterator tranI = transitions.iterator(); tranI.hasNext();) {
            Transition currTransition = (Transition) tranI.next();
            //System.out.println(currTransition);
            if (oldNodes.contains(currTransition.source)) {
                currTransition.source = newNode;
            }
            if (oldNodes.contains(currTransition.dest)) {
                currTransition.dest = newNode;
            }
            transitionArray[j++] = currTransition;
        }

        // put back in transitions. Set behaviour will remove duplicates.
        transitions.clear();
        //System.out.println("j is " + j);
        for (int i=0; i < j; i++) {
            // put back in all non-silly transitions: epsilon self-transitions are silly.
            Transition currTransition = transitionArray[i];

            if ((currTransition.sym != EPSILON) ||
                (currTransition.source != currTransition.dest))
                transitions.add(currTransition);
        }

        // now clean up the nodes
        //System.out.println("removing oldNodes");
        for (Iterator oldI = oldNodes.iterator(); oldI.hasNext(); ) {
            nodes.remove(oldI.next());
        }

        nodes.add(newNode);
    }

}

