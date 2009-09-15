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

import java.util.Iterator;

import junit.framework.TestCase;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.RNATools;

public class TranslationTest extends TestCase {
    private SymbolList dnaForm;
    private SymbolList rnaForm;
    private SymbolList pepForm;

    public TranslationTest(String name) {
	super(name);
    }

    protected void setUp() throws Exception {
	dnaForm = DNATools.createDNA("tttcctgtc");
	rnaForm = RNATools.createRNA("uuuccuguc");
	pepForm = ProteinTools.createProtein("FPV");
    }


    public void testTranscribe()
        throws Exception
    {
	assertTrue(compareSymbolList(RNATools.transcribe(dnaForm), rnaForm));
    }

    public void testTranslate()
        throws Exception
    {
	assertTrue(compareSymbolList(RNATools.translate(rnaForm), pepForm));
    }

    private boolean compareSymbolList(SymbolList sl1, SymbolList sl2) {
	if (sl1.length() != sl2.length()) {
	    return false;
	}

	Iterator si1 = sl1.iterator();
	Iterator si2 = sl2.iterator();
	while (si1.hasNext()) {
	    if (! (si1.next() == si2.next())) {
		return false;
	    }
	}

	return true;
    }
}
