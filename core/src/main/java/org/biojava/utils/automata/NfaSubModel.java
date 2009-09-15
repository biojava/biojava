



package org.biojava.utils.automata;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.symbol.Symbol;

/**
 * This class caches a reference to all objects that
 * it directs its delegate to make.  These references
 * make it relatively easy for it to duplicate
 * all objects made through this class.
 */
public class NfaSubModel
    implements NfaBuilder
{
    private NfaBuilder delegate;
    private Set nodes = new HashSet();
    private Set transitions = new HashSet();

    private FiniteAutomaton.Node start = null;
    private FiniteAutomaton.Node end = null;

    NfaSubModel(NfaBuilder delegate)
    {
        this.delegate = delegate;
        start = addNode(false);
        end = addNode(false);
    }

    public FiniteAutomaton getAutomaton()
    {
        return delegate.getAutomaton();
    }

    public FiniteAutomaton.Node getStart() { return start; }
    public FiniteAutomaton.Node getEnd() { return end; }

    void setStart(FiniteAutomaton.Node start)
    {
        nodes.add(start);
        this.start = start;
    }

    void setEnd(FiniteAutomaton.Node end)
    {
        nodes.add(end);
        this.end = end;
    }

    public FiniteAutomaton.Node addNode(boolean isTerminal)
    {
        FiniteAutomaton.Node newNode = delegate.addNode(isTerminal);
        nodes.add(newNode);
        return newNode;
    }

    public FiniteAutomaton.Transition addTransition(FiniteAutomaton.Node start, FiniteAutomaton.Node end, Symbol sym)
    {
        FiniteAutomaton.Transition newTransition = delegate.addTransition(start, end, sym);
        transitions.add(newTransition);
        return newTransition;
    }

    public FiniteAutomaton.Transition addEpsilonTransition(FiniteAutomaton.Node start, FiniteAutomaton.Node end)
    {
        FiniteAutomaton.Transition newTransition = delegate.addEpsilonTransition(start, end);
        transitions.add(newTransition);
        return newTransition;
    }

    public FiniteAutomaton.Transition addLambdaTransition(FiniteAutomaton.Node start, FiniteAutomaton.Node end)
    {
        FiniteAutomaton.Transition newTransition = delegate.addLambdaTransition(start, end);
        transitions.add(newTransition);
        return newTransition;
    }

    public FiniteAutomaton.NodeSet getNodes()
    {
        FiniteAutomaton.NodeSet nodeSet = delegate.createNodeSet();
        nodeSet.addAll(nodes);
        return nodeSet;
    }

    public Set getTransitions()
    {
        return transitions;
    }

    public FiniteAutomaton.NodeSet createNodeSet()
    {
        return delegate.createNodeSet();
    }

    /**
     * Makes a deep clone of this instance.
     */
    public NfaSubModel duplicate()
    {
        // create a new NfaSubModel
        NfaSubModel newSubModel = new NfaSubModel(delegate);

        // create a mapping between old and new Nodes
        Map old2New = new HashMap();

        for (Iterator nodesI = nodes.iterator(); nodesI.hasNext(); ) {
            FiniteAutomaton.Node oldNode = (FiniteAutomaton.Node) nodesI.next();
            old2New.put(oldNode, newSubModel.addNode(oldNode.isTerminal()));
        } 

        // set up the start and end
        newSubModel.setStart((FiniteAutomaton.Node) old2New.get(start));
        newSubModel.setEnd((FiniteAutomaton.Node) old2New.get(end));

        // for each transition, create a new one with remapped
        // nodes.
        for (Iterator transI = transitions.iterator(); transI.hasNext(); ) {
            FiniteAutomaton.Transition oldTrans = (FiniteAutomaton.Transition) transI.next();

            newSubModel.addTransition((FiniteAutomaton.Node) old2New.get(oldTrans.source), 
                (FiniteAutomaton.Node) old2New.get(oldTrans.dest), oldTrans.sym);
        }

        return newSubModel;
    }

    public void append(NfaSubModel submodel)
    {
        addEpsilonTransition(end, submodel.getStart());
        setEnd(submodel.getEnd());
    }

    public String toString() { return delegate.toString(); }
}

