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

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.Frame;
import org.biojava.bio.seq.FramedFeature;
import org.biojava.bio.seq.RNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.IllegalAlphabetException;

/**
 * Title:        SimpleFramedFeature.<p>
 * Description:  A no frills implementation of FramedFeature.<p>
 * Copyright:    Copyright (c) 2001.<p>
 *
 * @author Mark Schreiber
 * @version 1.0
 */

public class SimpleFramedFeature extends SimpleStrandedFeature implements FramedFeature, Frame {
  private FramedFeature.ReadingFrame readingFrame;

  public SimpleFramedFeature(Sequence sourceSeq, FeatureHolder parent, FramedFeature.Template template)
    throws IllegalAlphabetException {
    super(sourceSeq,parent,template);
    this.readingFrame = template.readingFrame;
    if (sourceSeq.getAlphabet() == RNATools.getRNA() && template.strand == NEGATIVE) {
      throw new IllegalAlphabetException("Cannot create a FramedFeature on the negative strand of an RNA");
    }
    else if (sourceSeq.getAlphabet() != DNATools.getDNA()) {
      //throw new IllegalAlphabetException("Cannot create a FramedFeature on a sequence of type "+
      //                                    sourceSeq.getAlphabet().getName());
    }
  }
  public ReadingFrame getReadingFrame() {
    return readingFrame;
  }

  public int getFrame(){
    return readingFrame.getFrame();
  }

  public Feature.Template makeTemplate(){
    FramedFeature.Template ft = new FramedFeature.Template();
    fillTemplate(ft);
    return ft;
  }

  protected void fillTemplate(FramedFeature.Template ft){
    super.fillTemplate(ft);
	ft.readingFrame = getReadingFrame();
  }
  public String toString(){
   return super.toString() + " "+getStrand().getToken()+""
           + getReadingFrame().getFrame();
  }

}
