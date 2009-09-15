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

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.Symbol;

/**
 * Base-class for builders that pass filtered events onto another builder.
 *
 * @author Matthew Pocock
 * @since 1.1
 * @see org.biojavax.bio.seq.io.RichSeqIOAdapter
 */
public class SequenceBuilderFilter implements SequenceBuilder {
  /**
   * Delegate that will recieve all of the forwarded event.
   */
  private SequenceBuilder delegate;
  
  /**
   * Create a new SeqIOFilter that will forward events on to another listener.
   *
   * @param delegate  the SequenceBuilder to wrap
   */
  public SequenceBuilderFilter(SequenceBuilder delegate) {
    this.delegate = delegate;
  }
  
  /**
   * Retrieve the delegate that is wrapped.
   *
   * @return the SequenceBuilder delegate
   */
  public SequenceBuilder getDelegate() {
    return delegate;
  }
  
  public void startSequence() throws ParseException {
    delegate.startSequence();
  }
  
  public void endSequence() throws ParseException {
    delegate.endSequence();
  }
  
  public void setName(String name) throws ParseException {
    delegate.setName(name);
  }
  
  public void setURI(String uri) throws ParseException {
    delegate.setURI(uri);
  }
  
  public void addSymbols(Alphabet alpha, Symbol[] syms, int start, int length)
  throws IllegalAlphabetException {
    delegate.addSymbols(alpha, syms, start, length);
  }
  
  public void addSequenceProperty(Object key, Object value) throws ParseException {
    delegate.addSequenceProperty(key, value);
  }
  
  public void startFeature(Feature.Template templ) throws ParseException {
    delegate.startFeature(templ);
  }
  
  public void endFeature() throws ParseException {
    delegate.endFeature();
  }
  
  public void addFeatureProperty(Object key, Object value) throws ParseException {
    delegate.addFeatureProperty(key, value);
  }
  
  public Sequence makeSequence() throws BioException {
    return delegate.makeSequence();
  }
}
