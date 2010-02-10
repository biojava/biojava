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

package org.biojava.bio.seq.homol;

import junit.framework.TestCase;

import org.biojava.bio.Annotation;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.utils.ChangeListener;

public class SimpleSimilarityPairFeatureTest extends TestCase
{
    protected Sequence qSeq;
    protected Sequence sSeq;

    protected SimilarityPairFeature qf;
    protected SimilarityPairFeature sf;

    public SimpleSimilarityPairFeatureTest(String name)
    {
        super(name);
    }

    protected void setUp() throws Exception
    {
        qSeq = new SimpleSequence(DNATools.createDNA("aacgtaggttccatgc"),
                                 "fragment1",
                                 "fragment1",
                                 Annotation.EMPTY_ANNOTATION);
        sSeq = new SimpleSequence(DNATools.createDNA("ttaacgtttttttttt"),
                                  "fragment2",
                                  "fragment2",
                                  Annotation.EMPTY_ANNOTATION);

        // Query sequence feature
        SimilarityPairFeature.Template qt =
            new SimilarityPairFeature.Template();
        qt.type       = "similarity";
        qt.source     = "test";
        qt.location   = new RangeLocation(1, 5);
        qt.strand     = StrandedFeature.POSITIVE;
        qt.score      = 1.0D;
        qt.annotation = Annotation.EMPTY_ANNOTATION;

        // Subject sequence feature
        SimilarityPairFeature.Template st =
            new SimilarityPairFeature.Template();
        st.type       = "similarity";
        st.source     = "test";
        st.location   = new RangeLocation(3, 7);
        st.strand     = StrandedFeature.POSITIVE;
        st.score      = 1.0D;
        st.annotation = Annotation.EMPTY_ANNOTATION;

        qt.alignment = SimilarityPairFeature.EMPTY_PAIRWISE;
        st.alignment = SimilarityPairFeature.EMPTY_PAIRWISE;

        qf = (SimilarityPairFeature) qSeq.createFeature(qt);
        sf = (SimilarityPairFeature) sSeq.createFeature(st);

        sf.setSibling(qf);
        qf.setSibling(sf);

        qf.addChangeListener(ChangeListener.ALWAYS_VETO);
        sf.addChangeListener(ChangeListener.ALWAYS_VETO);
    }

    public void testGetType()
    {
        assertEquals("similarity", qf.getType());
        assertEquals("similarity", sf.getType());
    }

    public void testGetSource()
    {
        assertEquals("test", qf.getSource());
        assertEquals("test", sf.getSource());
    }

    public void testGetLocation()
    {
         assertEquals(1, qf.getLocation().getMin());
         assertEquals(5, qf.getLocation().getMax());
         assertEquals(3, sf.getLocation().getMin());
         assertEquals(7, sf.getLocation().getMax());
    }

    public void testGetStrand()
    {
        assertEquals(StrandedFeature.POSITIVE, qf.getStrand());
        assertEquals(StrandedFeature.POSITIVE, sf.getStrand());
    }

    public void testGetScore()
    {
        assertEquals(1.0D, qf.getScore(), 0.0D);
        assertEquals(1.0D, sf.getScore(), 0.0D);
    }

    public void testGetAnnotation()
    {
        assertEquals(Annotation.EMPTY_ANNOTATION, qf.getAnnotation());
        assertEquals(Annotation.EMPTY_ANNOTATION, sf.getAnnotation());
    }

    public void testGetAlignment()
    {
        assertEquals(SimilarityPairFeature.EMPTY_PAIRWISE, qf.getAlignment());
        assertEquals(SimilarityPairFeature.EMPTY_PAIRWISE, sf.getAlignment());
    }

    public void testGetSibling()
    {
        assertEquals(qf, sf.getSibling());
        assertEquals(sf, qf.getSibling());
    }
}
