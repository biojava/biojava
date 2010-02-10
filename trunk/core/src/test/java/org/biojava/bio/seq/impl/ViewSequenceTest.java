package org.biojava.bio.seq.impl;

import junit.framework.TestCase;

import org.biojava.bio.Annotation;
import org.biojava.bio.SmallAnnotation;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceTools;
import org.biojava.bio.symbol.DummySymbolList;
import org.biojava.bio.symbol.LocationTools;
import org.biojava.ontology.OntoTools;

/**
 *
 *
 * @author Matthew Pocock
 */
public class ViewSequenceTest
        extends TestCase
{
  public void testEquivalentOriginals()
          throws Throwable
  {
    Sequence original = SequenceTools.createSequence(new DummySymbolList(DNATools.getDNA(), 1000),
                                                     "dummy",
                                                     "dummy",
                                                     Annotation.EMPTY_ANNOTATION);
    Feature.Template tplt = new Feature.Template();
    tplt.annotation = new SmallAnnotation();
    tplt.annotation.setProperty("geneID", "braca1");
    tplt.location = LocationTools.makeLocation(1,2);
    tplt.type = "synth";
    tplt.typeTerm = OntoTools.ANY;
    tplt.source = "my head";
    tplt.sourceTerm = OntoTools.ANY;

    Feature made = original.createFeature(tplt);
    ViewSequence view = new ViewSequence(original);
    Feature viewed = (Feature) view.features().next();
    Feature.Template copy = viewed.makeTemplate();

    System.out.println("template: " + tplt);
    System.out.println("orignal:  " + made.makeTemplate());
    System.out.println("copy:     " + copy);
    assertEquals("templates are equal", tplt, copy);
    assertEquals("made are equal", tplt, made.makeTemplate());
    assertEquals("both are equal", made.makeTemplate(), copy);
  }

  public void testEquivalentViews()
          throws Throwable
  {
    Sequence original = SequenceTools.createSequence(new DummySymbolList(DNATools.getDNA(), 1000),
                                                     "dummy",
                                                     "dummy",
                                                     Annotation.EMPTY_ANNOTATION);
    Feature.Template tplt = new Feature.Template();
    tplt.annotation = new SmallAnnotation();
    tplt.annotation.setProperty("geneID", "braca1");
    tplt.location = LocationTools.makeLocation(1, 2);
    tplt.type = "synth";
    tplt.typeTerm = OntoTools.ANY;
    tplt.source = "my head";
    tplt.sourceTerm = OntoTools.ANY;

    ViewSequence view = new ViewSequence(original);
    Feature made = view.createFeature(tplt);

    Feature viewed = (Feature) view.features().next();
    Feature.Template copy = viewed.makeTemplate();

    System.out.println("template: " + tplt);
    System.out.println("orignal:  " + made.makeTemplate());
    System.out.println("copy:     " + copy);
    assertEquals("templates are equal", tplt, copy);
    assertEquals("made are equal", tplt, made.makeTemplate());
    assertEquals("both are equal", made.makeTemplate(), copy);
  }
}