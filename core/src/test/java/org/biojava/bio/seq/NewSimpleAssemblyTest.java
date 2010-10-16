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

package org.biojava.bio.seq;

import java.util.Iterator;

import junit.framework.TestCase;

import org.biojava.bio.Annotation;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.bio.symbol.SymbolList;

/**
 * Tests for SimpleAssembly.  By dependancy, this also
 * tests ProjectedFeatureHolder and SimpleAssembly.
 *
 * @author Thomas Down
 * @since 1.2
 */

public class NewSimpleAssemblyTest extends TestCase
{
    protected Sequence fragment1;
    protected Sequence fragment2;
    protected Sequence fragment3;
    protected Sequence fragment4;
    protected Sequence fragment5;
    protected Sequence fragment6;
    protected Sequence assembly;

    public NewSimpleAssemblyTest(String name) {
	super(name);
    }

    protected void setUp() throws Exception {
	fragment1 = new SimpleSequence(DNATools.createDNA("aacgta"),
				       "fragment1",
				       "fragment1",
				       Annotation.EMPTY_ANNOTATION);
	fragment2 = new SimpleSequence(DNATools.createDNA("ttgatgc"),
				       "fragment2",
				       "fragment2",
				       Annotation.EMPTY_ANNOTATION);
        fragment3 = new SimpleSequence(DNATools.createDNA("aagctt"),
                                       "fragment3",
                                       "fragment3",
                                       Annotation.EMPTY_ANNOTATION);
        fragment4 = new SimpleSequence(DNATools.createDNA("ttagat"),
                                       "fragment4",
                                       "fragment4",
                                       Annotation.EMPTY_ANNOTATION);
        fragment5 = new SimpleSequence(DNATools.createDNA("ggatcc"),
                                       "fragment5",
                                       "fragment5",
                                       Annotation.EMPTY_ANNOTATION);
        fragment6 = new SimpleSequence(DNATools.createDNA("atctgg"),
                                       "fragment6",
                                       "fragment6",
                                       Annotation.EMPTY_ANNOTATION);

	StrandedFeature.Template sft = new StrandedFeature.Template();
	sft.type = "project_test";
	sft.source = "test";
	sft.annotation = Annotation.EMPTY_ANNOTATION;
	sft.strand = StrandedFeature.NEGATIVE;
	sft.location = new RangeLocation(1, 3);
	fragment2.createFeature(sft);

	assembly = new NewSimpleAssembly(31, "test", "test");
	ComponentFeature.Template templ = new ComponentFeature.Template();
	templ.type = "component";
	templ.source = "test";
	templ.annotation = Annotation.EMPTY_ANNOTATION;

	templ.strand = StrandedFeature.POSITIVE;
	templ.location = new RangeLocation(1, 4);
	templ.componentSequence = fragment1;
	templ.componentLocation = new RangeLocation(2, 5);
	assembly.createFeature(templ);

	templ.strand = StrandedFeature.NEGATIVE;
	templ.location = new RangeLocation(6, 12);
	templ.componentSequence = fragment2;
	templ.componentLocation = new RangeLocation(1, fragment2.length());
	assembly.createFeature(templ);

        templ.strand = StrandedFeature.POSITIVE;
        templ.location = new RangeLocation(14, 19);
        templ.componentSequence = fragment3;
        templ.componentLocation = new RangeLocation(1, fragment3.length());
        assembly.createFeature(templ);

        templ.strand = StrandedFeature.POSITIVE;
        templ.location = new RangeLocation(18, 23);
        templ.componentSequence = fragment4;
        templ.componentLocation = new RangeLocation(1, fragment4.length());
        assembly.createFeature(templ);

        templ.strand = StrandedFeature.POSITIVE;
        templ.location = new RangeLocation(26, 31);
        templ.componentSequence = fragment5;
        templ.componentLocation = new RangeLocation(1, fragment5.length());
        assembly.createFeature(templ);

        templ.strand = StrandedFeature.POSITIVE;
        templ.location = new RangeLocation(22, 27);
        templ.componentSequence = fragment6;
        templ.componentLocation = new RangeLocation(1, fragment6.length());
        assembly.createFeature(templ);
    }

    public void testAssembledSymbols()
	throws Exception
    {
	assertTrue(compareSymbolList(assembly,
				     DNATools.createDNA("acgtngcatcaanaagcttagatctggatcc")));
	assertTrue(compareSymbolList(assembly.subList(1,4),
				     DNATools.createDNA("acgt")));
	assertTrue(compareSymbolList(assembly.subList(1,5),
				     DNATools.createDNA("acgtn")));
	assertTrue(compareSymbolList(assembly.subList(5,5),
				     DNATools.createDNA("n")));
	assertTrue(compareSymbolList(assembly.subList(11, 12),
				     DNATools.createDNA("aa")));
    }

    public void testProjectedFeatures()
        throws Exception
    {
	FeatureHolder f = assembly.filter(new FeatureFilter.ByType("project_test"), true);
	assertEquals(f.countFeatures(), 1);
	Feature pf = (Feature) f.features().next();
	Location pfl = pf.getLocation();
	assertEquals(pfl.getMin(), 10);
	assertEquals(pfl.getMax(), 12);
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
