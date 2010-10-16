/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */

package org.biojava.utils.automata;

import java.util.Iterator;
import java.util.Set;

import org.biojava.bio.symbol.AlphabetIndex;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;


public class ArrayStateMachineToolkit
    implements StateMachineToolkit
{
    /**
     * A StateMachine implementation in which the
     * transitions are maintained in an array. While
     * every ArrayStateMachine is initialised with a
     * FiniteAlphabet, symbol data is provided as an
     * AlphabetIndex int value for effiiciency
     * reasons and callers are responsible for validating
     * all such data prior to passing it to an instance
     * of this class.  It is not recommended that this
     * class is used with Alphabets with a large
     * number of Symbols as it uses a dense State
     * matrix implementation.
     *
     * @author David Huen
     * @since 1.4
     */
    private static class ArrayStateMachine
        implements StateMachineFactory
    {
        private static byte ERROR_STATE = Byte.MIN_VALUE;

        class Instance
            implements StateMachineInstance
        {
            private int statePointer = 0;
            private int start;
            private int end;
    
            private Instance(int start, int end, int statePointer)
            {
                this.start = start;
                this.end = end;
                this.statePointer = statePointer;
            }

            // this method should only be used by Comparators and
            // the equals() method.
            public int getStart() { return start; }

            /**
             * invoke transition from current state indicated by Symbol represented by symbol index.
             * @param symIdx alphabet index of the symbol encountered.
             * @return true if the symbol is valid and this state machine
             *         should continue to receive input.
             */
            public boolean transit(int symIdx)
            {
                byte dest = transitions[entryPoints[statePointer] + symIdx];
                if (dest == ERROR_STATE) return false;
                end++;
                statePointer = Math.abs(dest);
    
                if (dest <0) {
                    // notify listener
                    listener.notifyHit(name, start, end, false);
    
                    return false;
                }
                else
                    return true;
            }

            /**
             * Two Instances are equal if they are children of the
             * same ArrayStateMachine instance and have the same
             * start value.
             */
            public boolean equals(Object o)
            {
                if (!(o instanceof Instance)) return false;

                StateMachineInstance other = (Instance) o;

                if (other.parent() != ArrayStateMachine.this)
                    return false;

                return (start == other.getStart());
            }

            public StateMachineFactory parent() { return  ArrayStateMachine.this; }

        }
    
        class GreedyInstance
            implements StateMachineInstance
        {
            private int statePointer;
            private int start;
            private int end;
            boolean gotTerminationState = false;
    
            private GreedyInstance(int start, int end, int statePointer)
            {
                this.start = start;
                this.statePointer = statePointer;
            }

            // this method should only be used by Comparators and
            // the equals() method.
            public int getStart() { return start; }
 
            /**
             * invoke transition from current state indicated by Symbol represented by symbol index.
             * @param symIdx alphabet index of the symbol encountered.
             * @return true if the symbol is valid and this state machine
             *         should continue to receive input.
             */
            public boolean transit(int symIdx)
            {
                byte dest = transitions[entryPoints[statePointer] + symIdx];
                if (dest == ERROR_STATE) {
                    if (gotTerminationState) {
                        listener.notifyHit(name, start, end, true);
                    }
    
                    return false;
                }
                else {
                    end++;
                    statePointer = Math.abs(dest);
    
                    if (dest <0) {
                        // got a valid termination state, save it
                        gotTerminationState = true;
                        return false;
                    }
                    else
                        return true;
                }
            }

            /**
             * Two GreedyInstances are equal if they are children of the
             * same ArrayStateMachine instance and have the same
             * start value.
             */
            public boolean equals(Object o)
            {
                if (!(o instanceof Instance)) return false;

                StateMachineInstance other = (Instance) o;

                if (other.parent() != ArrayStateMachine.this)
                    return false;

                return (start == other.getStart());
            }

            public StateMachineFactory parent() { return  ArrayStateMachine.this; }

        }
    
        private String name;
        private FiniteAlphabet alfa;
        private boolean greedy;

        private PatternListener listener = null;
    
        /**
         * maps a Node ID to the index in the
         * array where its transition data is stored.
         */
        private int [] entryPoints;
    
        /**
         * transition matrix.  Organised as:-
         * [source nodeID * symbol idx] -&gt; dest nodeID.
         */
        private byte [] transitions;
    
        /**
         * Creates an ArrayStateMachine from a
         * FiniteAutomaton instance. The FiniteAutomaton
         * must be in its most compact form (node IDs
         * in contiguous running order, no duplicate/surplus
         * nodes/transitions).
         */
        ArrayStateMachine(String name, FiniteAutomaton fa, boolean greedy)
        {
            this.name = name;
            alfa = fa.getAlphabet();
            int alfaSize = alfa.size();
            entryPoints = new int [fa.nodeCount()];
            transitions = new byte [alfaSize * fa.nodeCount()];
            this.greedy = greedy;

            convert(fa);
            fa = null;   // release FA for GC.
        }
    
        public void setListener(PatternListener listener)
        {
            this.listener = listener;
        }
    
        private void convert(FiniteAutomaton fa)
        {
            try {
                AlphabetIndex alfaIdx = AlphabetManager.getAlphabetIndex(alfa);
        
                int idx = 0;
                int alfaSize = alfa.size();
        
                // initialise pointer array
                for (int i = 0; i < entryPoints.length; i++) {
                    entryPoints[i] = idx;
                    idx += alfaSize;
                }
        
                // initialise transition matrix
                for (int i = 0; i < transitions.length; i++) {
                    transitions[i] = ERROR_STATE;
                }
        
                // go thru all transitions, filling the transition matrix.
                Set transitionSet = fa.getTransitions();
        
                for (Iterator transI = transitionSet.iterator(); transI.hasNext(); ) {
                    FiniteAutomaton.Transition currTransition = (FiniteAutomaton.Transition) transI.next();
        
                    int symIdx = alfaIdx.indexForSymbol(currTransition.getSymbol());
                    //System.out.println(currTransition.getSymbol() + " " + symIdx + " " + currTransition.getSource().getID() + " " + currTransition.getDest().getID());
                    transitions[entryPoints[Math.abs(currTransition.getSource().getID())] + symIdx]
                        = (byte) currTransition.getDest().getID();
                }
            }
            catch (IllegalSymbolException ise) {
                throw new AssertionError(ise);
            }
        } 
    
        /**
         * get a StateMachineInstance to parse the sequence with.
         * @param start current Sequence coordinate.
         * @param greedy should greedy regex semantics be used?
         */
        StateMachineInstance getInstance(int start)
        {
            if (greedy)
                return new GreedyInstance(start, start, 0);
            else
                return new Instance(start, start, 0);
        }
    
        /**
         * Return a StateMachineInstance if the Symbol represented
         * by the symbol index is valid as the initial symbol of
         * the pattern.  This method must remain package private
         * as it does no alphabet checks at all.
         *
         * @param symIdx alphabet index value for specified symbol.
         * @param greedy should greedy regex semantics be used?
         * @return an instance of StateMachineInstance if symbol
         * is valid otherwise null.
         */
        public StateMachineInstance startInstance(int symIdx, int start)
        {
            int nextState = transitions[symIdx];
            //System.out.println("startInstance called with " + symIdx + " " + start + " " + nextState);
            if (nextState != ERROR_STATE) {
                if (greedy) {
                    return new GreedyInstance(start, start + 1, nextState);
                }
                else {
                    //System.out.println("creating new Instance.");
                    return new Instance(start, start + 1, nextState);
                }
            }
            else
                return null;
        }
    }

    private FiniteAlphabet alfa;
    private boolean greedy;

    ArrayStateMachineToolkit(FiniteAlphabet alfa, boolean greedy)
    {
        this.alfa = alfa;
        this.greedy = greedy;
    }

    public StateMachineFactory getFactory(String factoryName, FiniteAutomaton fa)
        throws IllegalAlphabetException
    {
        if (alfa != fa.getAlphabet())
            throw new IllegalAlphabetException("The FiniteAutomaton is defined on an Alphabet incompatible with this ArrayStateMachineToolKit.");
        return new ArrayStateMachine(factoryName, fa, greedy);
    }
}

