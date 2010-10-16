


package org.biojava.utils.automata;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.ListIterator;

import org.biojava.bio.BioException;
import org.biojava.bio.symbol.AlphabetIndex;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;

public class PatternBlitz
{
    FiniteAlphabet alfa;

    private LinkedList patternStore;
    private LinkedList instanceStore;
    private PatternListener listener;
    private StateMachineToolkit factory;
    private AlphabetIndex alfaIdx;

    public PatternBlitz(FiniteAlphabet alfa, StateMachineToolkit factory)
    {
        this.alfa = alfa;
        alfaIdx = AlphabetManager.getAlphabetIndex(alfa);
        patternStore = new LinkedList();
        instanceStore = new LinkedList();
        this.factory = factory;
    }

    public void lock()
    {

    }

    public void setListener(PatternListener listener)
    {
        this.listener = listener;
    }

    /**
     * add the specified regex to the patterns
     * used for searching.
     */
    public void addPattern(String pattern)
    {
        try {
            // create DFA for pattern
            FiniteAutomaton dfa = PatternMaker.compilePattern(pattern, alfa);

            // convert to state machine
            StateMachineFactory instanceFactory = factory.getFactory(pattern, dfa);
            instanceFactory.setListener(listener);

            // save in pattern store
            patternStore.addLast(instanceFactory);
        }
        catch (BioException be) {
            throw new AssertionError(be);
        }
    }

    private void scanPatterns(Symbol sym, int start)
        throws IllegalSymbolException
    {
        // convert symbol to its index
        int symIdx = alfaIdx.indexForSymbol(sym);

        // go thru' instance store with symbol index
        for (ListIterator instanceI = instanceStore.listIterator();
            instanceI.hasNext(); ) {

            // remove all terminal entries
            StateMachineInstance instance = (StateMachineInstance) instanceI.next();
            if (!instance.transit(symIdx))
                instanceI.remove();
        }

        // now traverse the pattern store to initiate new instances
        for (Iterator patternI = patternStore.iterator();
            patternI.hasNext(); ) {

            StateMachineFactory factory = (StateMachineFactory) patternI.next();

            StateMachineInstance newInstance;
            if ((newInstance = factory.startInstance(symIdx, start)) != null) {
                instanceStore.addLast(newInstance);
            }
        }
    }

    public void search(SymbolList sl)
        throws IllegalAlphabetException
    {
        // check compatible alphabets
        if (sl.getAlphabet() != alfa)
            throw new IllegalAlphabetException("incompatible alphabets");

        try {
            for (int i=1; i < sl.length(); i++) {
                Symbol sym = sl.symbolAt(i);

                scanPatterns(sym, i);
            }
        }
        catch (IllegalSymbolException ise) {
            throw new AssertionError(ise);
        }
    }

}

