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

package org.biojava.bio.seq.io;

import org.biojava.bio.seq.Feature;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.Symbol;

/**
 * Adapter class for SeqIOListener that has empty methods.
 *
 * @author Matthew Pocock
 * @since 1.1
 * @see org.biojavax.bio.seq.io.RichSeqIOAdapter
 */

public class SeqIOAdapter implements SeqIOListener {
  public void startSequence()
  throws ParseException {}
  
  public void endSequence()
  throws ParseException {}
  
  public void setName(String name)
  throws ParseException {}
  
  public void setURI(String uri)
  throws ParseException {}
  
  public void addSymbols(Alphabet alpha, Symbol[] syms, int start, int length)
  throws IllegalAlphabetException {}
  
  public void addSequenceProperty(Object key, Object value)
  throws ParseException {}
  
  public void startFeature(Feature.Template templ)
  throws ParseException {}
  
  public void endFeature()
  throws ParseException {}
  
  public void addFeatureProperty(Object key, Object value)
  throws ParseException {}
}
