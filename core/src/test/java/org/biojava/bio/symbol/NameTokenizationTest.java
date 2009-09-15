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
import java.util.Iterator;
import java.util.List;

import junit.framework.TestCase;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.io.SymbolTokenization;

public class NameTokenizationTest extends TestCase {
    private FiniteAlphabet simple_alpha;
    private SymbolList simple_atomicSymbols;
    private SymbolList simple_allSymbols;
    private SymbolTokenization nameTokenization;

    public NameTokenizationTest(String name) {
	super(name);
    }

    protected void setUp() throws Exception {
	simple_alpha = DNATools.getDNA();

	List atomics = new ArrayList();
	for (Iterator i = simple_alpha.iterator(); i.hasNext(); ) {
	    atomics.add(i.next());
	}
	simple_atomicSymbols = new SimpleSymbolList(simple_alpha, atomics);

	List all = new ArrayList(AlphabetManager.getAllSymbols(simple_alpha));
	simple_allSymbols = new SimpleSymbolList(simple_alpha, all);

	nameTokenization = simple_alpha.getTokenization("name");
    }

    public void testAtomicNames()
        throws Exception
    {
	String tokenization = nameTokenization.tokenizeSymbolList(simple_atomicSymbols);
	SymbolList reconstruction = new SimpleSymbolList(nameTokenization, tokenization);
	assertTrue(SymbolUtils.compareSymbolLists(simple_atomicSymbols, reconstruction));
    }

    public void testAmbiguousNames()
        throws Exception
    {
	String tokenization = nameTokenization.tokenizeSymbolList(simple_allSymbols);
	SymbolList reconstruction = new SimpleSymbolList(nameTokenization, tokenization);
	assertTrue(SymbolUtils.compareSymbolLists(simple_allSymbols, reconstruction));
    }

}
