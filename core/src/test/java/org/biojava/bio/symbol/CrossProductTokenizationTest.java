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

package org.biojava.bio.symbol;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import junit.framework.TestCase;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.io.SymbolTokenization;

/**
 * Test for tokenization of cross-product symbols.
 *
 * @author Thomas Down
 * @since 1.2
 */

public class CrossProductTokenizationTest extends TestCase {
    private FiniteAlphabet simple_alpha;
    private FiniteAlphabet cross_alpha;
    private SymbolList cross_atomicSymbols;
    private SymbolTokenization nameTokenization;

    public CrossProductTokenizationTest(String name) {
	super(name);
    }

    protected void setUp() throws Exception {
	simple_alpha = DNATools.getDNA();
	cross_alpha = (FiniteAlphabet) AlphabetManager.getCrossProductAlphabet(Collections.nCopies(2, simple_alpha));

	List atomics = new ArrayList();
	for (Iterator i = cross_alpha.iterator(); i.hasNext(); ) {
	    atomics.add(i.next());
	}
	cross_atomicSymbols = new SimpleSymbolList(cross_alpha, atomics);

	nameTokenization = cross_alpha.getTokenization("name");
    }

    public void testAtomicNames()
        throws Exception
    {
	String tokenization = nameTokenization.tokenizeSymbolList(cross_atomicSymbols);
	SymbolList reconstruction = new SimpleSymbolList(nameTokenization, tokenization);
	assertTrue(SymbolUtils.compareSymbolLists(cross_atomicSymbols, reconstruction));
    }

//      public void testAmbiguousNames()
//          throws Exception
//      {
//  	String tokenization = nameTokenization.tokenizeSymbolList(simple_allSymbols);
//  	System.out.println(tokenization);
//  	SymbolList reconstruction = AlphabetManager.parse(nameTokenization, tokenization);
//  	assertTrue(SymbolUtils.compareSymbolLists(simple_allSymbols, reconstruction));
//      }

}
