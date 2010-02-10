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

package org.biojava.bio.seq.impl;

import java.util.Iterator;
import java.util.SortedSet;
import java.util.TreeSet;

import junit.framework.TestCase;

import org.biojava.bio.Annotation;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.GappedSequence;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.RangeLocation;

/**
 * Test for GappedSequence.  By dependancy, this also tests
 * some functionality from GappedSymbolList.
 *
 * @author Thomas Down
 * @since 1.3
 */

public class GappedSequenceTest extends TestCase
{
    protected Sequence seq;
    protected GappedSequence gappedSeq;

    public GappedSequenceTest(String name) {
         super(name);
    }

    protected void setUp() throws Exception {
		seq = new SimpleSequence(DNATools.createDNA("aacgtaggttccatgc"),
					       "fragment1",
					       "fragment1",
					       Annotation.EMPTY_ANNOTATION);
		
		Feature.Template sft = new Feature.Template();
		sft.type = "normal";
		sft.source = "test";
		sft.annotation = Annotation.EMPTY_ANNOTATION;
		sft.location = new RangeLocation(8, 10);
		seq.createFeature(sft);
	
		sft.type = "split";
		sft.location = new RangeLocation(2, 10);
		seq.createFeature(sft);
	
		gappedSeq = new SimpleGappedSequence(seq);
		gappedSeq.addGapsInSource(5, 2);
    }

    /**
     * Ensure that a location wholly contained by a contiguous
     * block in the gapped view is simply translated to the
     * appropriate position
     */

    public void testNonSplit()
	throws Exception
    {
	FeatureHolder fh = gappedSeq.filter(new FeatureFilter.ByType("normal"), false);
	Feature f = (Feature) fh.features().next();
	Location fl = f.getLocation();
	assertTrue(fl.isContiguous());
	assertEquals(fl.getMin(), 10);
	assertEquals(fl.getMax(), 12);
    }

    /**
     * Ensure that a location straddling a gap is projected
     * to a non-contiguous location, with appropriate
     * boundaries.
     */

    public void testSplit()
	throws Exception
    {
	FeatureHolder fh = gappedSeq.filter(new FeatureFilter.ByType("split"), false);
	Feature f = (Feature) fh.features().next();
	Location fl = f.getLocation();
	assertTrue(!fl.isContiguous());
	assertEquals(fl.getMin(), 2);
	assertEquals(fl.getMax(), 12);

	Iterator bi = fl.blockIterator();
	SortedSet sb = new TreeSet(Location.naturalOrder);
	while (bi.hasNext()) {
	    sb.add(bi.next());
	}
	assertEquals(sb.size(), 2);

	bi = sb.iterator();

	Location block = (Location) bi.next();
	assertEquals(block.getMin(), 2);
	assertEquals(block.getMax(), 4);

	block = (Location) bi.next();
	assertEquals(block.getMin(), 7);
	assertEquals(block.getMax(), 12);
    }

    public void testRemoveRemoteFeature()
    		throws Exception
    	{
		seq = new SimpleSequence(DNATools.createDNA("aacgtaggttccatgc"),
			       "fragment1",
			       "fragment1",
			       Annotation.EMPTY_ANNOTATION);

		Feature.Template sft = new Feature.Template();
		sft.type = "normal";
		sft.source = "test";
		sft.annotation = Annotation.EMPTY_ANNOTATION;
		sft.location = new RangeLocation(8, 10);
		seq.createFeature(sft);
		
		sft.type = "split";
		sft.location = new RangeLocation(2, 10);
		seq.createFeature(sft);
		
		gappedSeq = new SimpleGappedSequence(seq);
		gappedSeq.addGapsInSource(5, 2);
		
		Feature split = (Feature) gappedSeq.filter(new FeatureFilter.ByType("split")).features().next();
		assertEquals(seq.countFeatures(), 2);
		gappedSeq.removeFeature(split);
		assertEquals(seq.countFeatures(), 1);
    	}
    
    public void testRemoveLocalFeature()
		throws Exception
	{
        gappedSeq = new SimpleGappedSequence(seq);
        gappedSeq.addGapsInSource(5, 2);
        
		Feature.Template sft = new Feature.Template();
		sft.type = "local";
		sft.source = "test";
		sft.annotation = Annotation.EMPTY_ANNOTATION;
		sft.location = new RangeLocation(8, 10);
		seq.createFeature(sft);
		
		Feature local = (Feature) gappedSeq.filter(new FeatureFilter.ByType("local")).features().next();
		assertEquals(gappedSeq.countFeatures(), 3);
		gappedSeq.removeFeature(local);
		assertEquals(gappedSeq.countFeatures(), 2);
	}
}
