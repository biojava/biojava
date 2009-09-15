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

package org.biojava.bio.seq.project;

import junit.framework.TestCase;

import org.biojava.bio.Annotation;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FilterUtils;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.seq.projection.ProjectedFeatureHolder;
import org.biojava.bio.seq.projection.ReparentContext;
import org.biojava.bio.seq.projection.TranslateFlipContext;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * Tests for ProjectedFeatureHolder
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @since 1.3
 */

public class ProjectedFeatureHolderTest extends TestCase
{
    public ProjectedFeatureHolderTest(String name) {
        super(name);
    }

  public void testReparentContext()
  throws Exception
  {
    Sequence seq = new SimpleSequence(
            DNATools.createDNA("gattaca"),
            "test",
            "test",
            Annotation.EMPTY_ANNOTATION
    );
    Sequence seq2 = new SimpleSequence(
            DNATools.createDNA("-------"),
            "test",
            "test",
            Annotation.EMPTY_ANNOTATION
    );
    Feature.Template template = new Feature.Template();
    template.type = "test";
    template.source = "foo";
    template.location = new RangeLocation(2, 4);
    template.annotation = Annotation.EMPTY_ANNOTATION;
    Feature origFeat = seq.createFeature(template);

    ProjectedFeatureHolder pfh = new ProjectedFeatureHolder(
            new ReparentContext(seq2, seq));

    assertEquals("Same size", pfh.countFeatures(), seq.countFeatures());

    Feature projFeat = (Feature) pfh.features().next();
    assertEquals("getSource()", origFeat.getSource(), projFeat.getSource());
    assertEquals("getType()", origFeat.getType(), projFeat.getType());
    assertEquals("getLocation()", origFeat.getLocation(), projFeat.getLocation());
    assertEquals("origFeat.getParent()", origFeat.getParent(), seq);
    assertEquals("projFeat.getParent()", projFeat.getParent(), seq2);
  }

  public void testReparentContext_create()
  throws Exception
  {
    Sequence seq = new SimpleSequence(
            DNATools.createDNA("gattaca"),
            "test",
            "test",
            Annotation.EMPTY_ANNOTATION
    );
    Sequence seq2 = new SimpleSequence(
            DNATools.createDNA("-------"),
            "test",
            "test",
            Annotation.EMPTY_ANNOTATION
    );
    Feature.Template template = new Feature.Template();
    template.type = "test";
    template.source = "foo";
    template.location = new RangeLocation(2, 4);
    template.annotation = Annotation.EMPTY_ANNOTATION;

    ProjectedFeatureHolder pfh = new ProjectedFeatureHolder(
            new ReparentContext(seq2, seq));

    Feature projFeat = pfh.createFeature(template);

    Feature origFeat = (Feature) seq.features().next();

    assertEquals("Same size", pfh.countFeatures(), seq.countFeatures());

    assertEquals("getSource()", origFeat.getSource(), projFeat.getSource());
    assertEquals("getType()", origFeat.getType(), projFeat.getType());
    assertEquals("getLocation()", origFeat.getLocation(), projFeat.getLocation());
    assertEquals("origFeat.getParent()", origFeat.getParent(), seq);
    assertEquals("projFeat.getParent()", projFeat.getParent(), seq2);
  }

  public void testTranslateFlipContext_translateOnly()
          throws Exception {
    Sequence seq = new SimpleSequence(
            DNATools.createDNA("gattaca"),
            "test",
            "test",
            Annotation.EMPTY_ANNOTATION
    );
    Sequence seq2 = new SimpleSequence(
            DNATools.createDNA("-------"),
            "test",
            "test",
            Annotation.EMPTY_ANNOTATION
    );
    StrandedFeature.Template template = new StrandedFeature.Template();
    template.type = "test";
    template.source = "foo";
    template.location = new RangeLocation(2, 4);
    template.annotation = Annotation.EMPTY_ANNOTATION;
    template.strand = StrandedFeature.NEGATIVE;

    StrandedFeature origFeat = (StrandedFeature) seq.createFeature(template);

    ProjectedFeatureHolder pfh = new ProjectedFeatureHolder(
            new TranslateFlipContext(seq2, seq, 3, false));

    assertEquals("Same size", pfh.countFeatures(), seq.countFeatures());

    StrandedFeature projFeat = (StrandedFeature) pfh.features().next();

    assertEquals("getSource()", origFeat.getSource(), projFeat.getSource());
    assertEquals("getType()", origFeat.getType(), projFeat.getType());
    assertEquals("getLocation()", origFeat.getLocation().translate(3), projFeat.getLocation());
    assertEquals("getStrand()", origFeat.getStrand(), projFeat.getStrand());
    assertEquals("projFeat.getParent()", projFeat.getParent(), seq2);
    assertEquals("pfh.filter(5,7).size() == 1",
                 1,
                 pfh.filter(FilterUtils.containedByLocation(projFeat.getLocation())).countFeatures());
    assertTrue("pfh.filter(5,7).contains(projFeat)",
               pfh.filter(FilterUtils.containedByLocation(projFeat.getLocation())).containsFeature(projFeat));
    assertEquals("projFeat.filter(POSITIVE).size() == 1",
                 1,
                 pfh.filter(FilterUtils.byStrand(origFeat.getStrand())).countFeatures());
  }


  public void testTranslateFlipContext_translateAndFlip()
          throws Exception {
    Sequence seq = new SimpleSequence(
            DNATools.createDNA("gattaca"),
            "test",
            "test",
            Annotation.EMPTY_ANNOTATION
    );
    Sequence seq2 = new SimpleSequence(
            DNATools.createDNA("-------"),
            "test",
            "test",
            Annotation.EMPTY_ANNOTATION
    );
    StrandedFeature.Template template = new StrandedFeature.Template();
    template.type = "test";
    template.source = "foo";
    template.location = new RangeLocation(2, 3);
    template.annotation = Annotation.EMPTY_ANNOTATION;
    template.strand = StrandedFeature.NEGATIVE;

    StrandedFeature origFeat = (StrandedFeature) seq.createFeature(template);

    ProjectedFeatureHolder pfh = new ProjectedFeatureHolder(
            new TranslateFlipContext(seq2, seq, seq.length() + 1, true));

    assertEquals("Same size", pfh.countFeatures(), seq.countFeatures());

    StrandedFeature projFeat = (StrandedFeature) pfh.features().next();

    assertEquals("getSource()", origFeat.getSource(), projFeat.getSource());
    assertEquals("getType()", origFeat.getType(), projFeat.getType());
    assertEquals("getLocation()",
                 new RangeLocation(5,6),
                 projFeat.getLocation());
    assertEquals("getStrand()", StrandedFeature.POSITIVE, projFeat.getStrand());
    assertEquals("projFeat.getParent()", projFeat.getParent(), seq2);
    assertEquals("pfh.filter(5,6).size() == 1",
                 1,
                 pfh.filter(FilterUtils.containedByLocation(projFeat.getLocation())).countFeatures());
    assertTrue("pfh.filter(5,6).contains(projFeat)",
               pfh.filter(FilterUtils.containedByLocation(projFeat.getLocation())).containsFeature(projFeat));
    assertEquals("projFeat.filter(NEGATIVE).size() == 0",
                 0,
                 pfh.filter(FilterUtils.byStrand(StrandedFeature.NEGATIVE)).countFeatures());
    assertEquals("projFeat.filter(POSITIVE).size() == 1",
                 1,
                 pfh.filter(FilterUtils.byStrand(StrandedFeature.POSITIVE)).countFeatures());
  }

  public void testTranslateFlipContext_translateAndFlip_create()
          throws Exception {
    Sequence seq = new SimpleSequence(
            DNATools.createDNA("gattaca"),
            "test",
            "test",
            Annotation.EMPTY_ANNOTATION
    );
    Sequence seq2 = new SimpleSequence(
            DNATools.createDNA("-------"),
            "test",
            "test",
            Annotation.EMPTY_ANNOTATION
    );
    StrandedFeature.Template template = new StrandedFeature.Template();
    template.type = "test";
    template.source = "foo";
    template.location = new RangeLocation(5, 6);
    template.annotation = Annotation.EMPTY_ANNOTATION;
    template.strand = StrandedFeature.POSITIVE;

    ProjectedFeatureHolder pfh = new ProjectedFeatureHolder(
            new TranslateFlipContext(seq2, seq, seq.length() + 1, true));

    StrandedFeature projFeat = (StrandedFeature) pfh.createFeature(template);

    Feature origFeat = (StrandedFeature) seq.features().next();

    assertEquals("Same size", pfh.countFeatures(), seq.countFeatures());

    assertEquals("getSource()", origFeat.getSource(), projFeat.getSource());
    assertEquals("getType()", origFeat.getType(), projFeat.getType());
    assertEquals("getLocation()",
                 new RangeLocation(5,6),
                 projFeat.getLocation());
    assertEquals("getStrand()", StrandedFeature.POSITIVE, projFeat.getStrand());
    assertEquals("projFeat.getParent()", projFeat.getParent(), seq2);
    assertEquals("pfh.filter(5,6).size() == 1",
                 1,
                 pfh.filter(FilterUtils.containedByLocation(projFeat.getLocation())).countFeatures());
    assertTrue("pfh.filter(5,6).contains(projFeat)",
               pfh.filter(FilterUtils.containedByLocation(projFeat.getLocation())).containsFeature(projFeat));
    assertEquals("projFeat.filter(NEGATIVE).size() == 0",
                 0,
                 pfh.filter(FilterUtils.byStrand(StrandedFeature.NEGATIVE)).countFeatures());
    assertEquals("projFeat.filter(POSITIVE).size() == 1",
                 1,
                 pfh.filter(FilterUtils.byStrand(StrandedFeature.POSITIVE)).countFeatures());
  }

    public void testFeatureChangeEvent()
        throws Exception
    {
        Sequence seq = new SimpleSequence(
                DNATools.createDNA("gattaca"),
                "test",
                "test",
                Annotation.EMPTY_ANNOTATION
        );
        Feature.Template template = new Feature.Template();
        template.type = "test";
        template.source = "foo";
        template.location = new RangeLocation(2, 4);
        template.annotation = Annotation.EMPTY_ANNOTATION;
        Feature seqFeature = seq.createFeature(template);
        
        ProjectedFeatureHolder pfh = new ProjectedFeatureHolder(
                new TranslateFlipContext(seq, seq, 7, true));
        Feature pfhFeature = (Feature) pfh.filter(new FeatureFilter.ByType("test")).features().next();
        
        pfhFeature.addChangeListener(ChangeListener.ALWAYS_VETO, ChangeType.UNKNOWN);
        boolean vetoed = false;
        try {
        	synchronized(seqFeature) {
        		seqFeature.setLocation(new RangeLocation(1, 3));
        	}
        } catch (ChangeVetoException cve) {
            vetoed = true;
        }
        assertTrue(vetoed);
        
        pfhFeature.removeChangeListener(ChangeListener.ALWAYS_VETO, ChangeType.UNKNOWN);
        seqFeature.setLocation(new RangeLocation(1, 3));
    }
}
