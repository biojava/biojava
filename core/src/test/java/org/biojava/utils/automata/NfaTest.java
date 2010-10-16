
package org.biojava.utils.automata;

import junit.framework.TestCase;

import org.biojava.bio.seq.DNATools;

public class NfaTest
    extends TestCase
{
    private static Nfa nfa;

    private Nfa.Node nT1LB;
    private Nfa.Node nT1preC;
    private Nfa.Node nT1postC;
    private Nfa.Node nT1preT;
    private Nfa.Node nT1postT;
    private Nfa.Node nT1RB;
    private Nfa.Node nT2LB;
    private Nfa.Node nT2preC1;
    private Nfa.Node nT2postC1;
    private Nfa.Node nT2preT1;
    private Nfa.Node nT2postT1;
    private Nfa.Node nT2preT2;
    private Nfa.Node nT2postT2;
    private Nfa.Node nT2preC2;
    private Nfa.Node nT2postC2;
    private Nfa.Node nT2preC3;
    private Nfa.Node nT2postC3;
    private Nfa.Node nT2preT3;
    private Nfa.Node nT2postT3;
    private Nfa.Node nT2RB;

    /**
     * Create the Nfa corresponding
     * to (c|t)*(ctt|c+t).
     */
    protected void setUp()
    {
        nfa = new Nfa("test Nfa", DNATools.getDNA());
        // (c|t)*
        nT1LB = nfa.addNode(false);
        nT1preC = nfa.addNode(false);
        nT1postC = nfa.addNode(false);
        nT1preT = nfa.addNode(false);
        nT1postT = nfa.addNode(false);
        nT1RB = nfa.addNode(false);

        nfa.addEpsilonTransition(nfa.getStart(), nT1LB);
        nfa.addEpsilonTransition(nT1LB, nT1preC);
        nfa.addEpsilonTransition(nT1LB, nT1preT);
        nfa.addTransition(nT1preC, nT1postC, DNATools.c());
        nfa.addTransition(nT1preT, nT1postT, DNATools.t());
        nfa.addEpsilonTransition(nT1postC, nT1RB);
        nfa.addEpsilonTransition(nT1postT, nT1RB);
        nfa.addEpsilonTransition(nfa.getStart(), nT1RB);
        nfa.addEpsilonTransition(nT1RB, nT1LB);

        // (ctt|c+t)
        nT2LB = nfa.addNode(false);
        nT2preC1 = nfa.addNode(false);
        nT2postC1 = nfa.addNode(false);
        nT2preT1 = nfa.addNode(false);
        nT2postT1 = nfa.addNode(false);
        nT2preT2 = nfa.addNode(false);
        nT2postT2 = nfa.addNode(false);
        nT2preC2 = nfa.addNode(false);
        nT2postC2 = nfa.addNode(false);
        nT2preC3 = nfa.addNode(false);
        nT2postC3 = nfa.addNode(false);
        nT2preT3 = nfa.addNode(false);
        nT2postT3 = nfa.addNode(false);
        nT2RB = nfa.addNode(false);

        nfa.addEpsilonTransition(nT1RB, nT2LB);

        nfa.addEpsilonTransition(nT2LB, nT2preC1);
        nfa.addTransition(nT2preC1, nT2postC1, DNATools.c());
        nfa.addEpsilonTransition(nT2postC1, nT2preT1);
        nfa.addTransition(nT2preT1, nT2postT1, DNATools.t());
        nfa.addEpsilonTransition(nT2postT1, nT2preT2);
        nfa.addTransition(nT2preT2, nT2postT2, DNATools.t());
        nfa.addEpsilonTransition(nT2postT2, nT2RB);

        nfa.addEpsilonTransition(nT2LB, nT2preC2);
        nfa.addTransition(nT2preC2, nT2postC2, DNATools.c());
        nfa.addEpsilonTransition(nT2postC2, nT2preC3);
        nfa.addTransition(nT2preC3, nT2postC3, DNATools.c());
        nfa.addEpsilonTransition(nT2postC3, nT2preC3);
        nfa.addEpsilonTransition(nT2postC3, nT2preT3);
        nfa.addTransition(nT2preT3, nT2postT3, DNATools.t());
        nfa.addEpsilonTransition(nT2postT3, nT2RB);

        nfa.addEpsilonTransition(nT2RB, nfa.getEnd());
    }

    /**
     * Test Nfa class by creating the Nfa corresponding
     * to (c|t)*(ctt|c+t).
     */
    public void testNfa()
    {
//        System.out.println(nfa.toString());
        nfa.doEpsilonClosure();
//        System.out.println(nfa.toString());
    }
}

