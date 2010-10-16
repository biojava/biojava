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

import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.SmallAnnotation;
import org.biojava.bio.seq.impl.DummySequence;
import org.biojava.bio.seq.impl.RevCompSequence;
import org.biojava.bio.seq.impl.SimpleGappedSequence;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.seq.impl.SubSequence;
import org.biojava.bio.seq.impl.ViewSequence;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.DummySymbolList;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeVetoException;

/**
 * Methods for manipulating sequences.
 *
 * @author Matthew Pocock
 */
public final class SequenceTools {
  private SequenceTools() {
  }
    
  public static Sequence createSequence(
    SymbolList syms, String uri, String name, Annotation ann
  ) {
    return new SimpleSequence(syms, uri, name, ann);
  }

  /**
   * Extract a sub-sequence from a sequence.
   *
   * <p>
   * The sub-sequence will be indexed from 1 through to (end-start+1). An index
   * of i in the sub-sequence corresponds to (i+start-1) in the original.
   * All features from the original sequence will be projected down into this
   * co-ordinate system. All features overlapping the edges will be given fuzzy
   * locations.
   * </p>
   *
   * @param seq   the sequence to sub-sequence
   * @param start the first index to include in the sub-sequence
   * @param end   the last index to include in the sub-sequence
   * @return a view Sequence for this region
   * @throws IndexOutOfBoundsException if start or end are not in seq, or if
   *     end < start
   */
  public static Sequence subSequence(Sequence seq, int start, int end)
  throws IndexOutOfBoundsException {
    return new SubSequence(seq, start, end);
  }

  /**
   * Extract a sub-sequence from a sequence.
   *
   * <p>
   * The sub-sequence will be indexed from 1 through to (end-start+1). An index
   * of i in the sub-sequence corresponds to (i+start-1) in the original.
   * All features from the original sequence will be projected down into this
   * co-ordinate system. All features overlapping the edges will be given fuzzy
   * locations.
   * </p>
   *
   * @param seq   the sequence to sub-sequence
   * @param start the first index to include in the sub-sequence
   * @param end   the last index to include in the sub-sequence
   * @param name  a new name to give to this sub-sequence
   * @return a view Sequence for this region
   * @throws IndexOutOfBoundsException if start or end are not in seq, or if
   *     end < start
   */
  public static Sequence subSequence(Sequence seq, int start, int end, String name)
  throws IndexOutOfBoundsException {
    return new SubSequence(seq, start, end, name);
  }

  /**
   * Extract a sub-sequence from a sequence.
   *
   * <p>
   * The sub-sequence will be indexed from 1 through to (end-start+1). If the
   * strand is NEGATIVE, all features will be flipped in the same manner as
   * the reverseComplement method. If it is UNKNOWN or
   * POSITIVE, then this is identical to the other subSequence methods.
   * </p>
   *
   * @param seq   the sequence to sub-sequence
   * @param start the first index to include in the sub-sequence
   * @param end   the last index to include in the sub-sequence
   * @param name  a new name to give to this sub-sequence
   * @param strand a StrandedFeature.Strand indicating which strand the
   *    sub-sequence should be on
   * @return a view Sequence for this region
   * @throws IndexOutOfBoundsException if start or end are not in seq, or if
   *     end < start
   */
  public static Sequence subSequence(
    Sequence seq,
    int start,
    int end,
    String name,
    StrandedFeature.Strand strand
  ) throws IndexOutOfBoundsException, IllegalAlphabetException {
    Sequence s = subSequence(seq, start, end, name);
    if(strand == StrandedFeature.NEGATIVE) {
      s = reverseComplement(s);
    }
    return s;
  }

  /**
   * Reverse-complement a sequence, and flip all of its features.
   *
   * @param seq  the Sequence to reverse-complement
   * @return  the flipped Sequence
   * @throws IllegalAlphabetException  if the symbols in the sequence can not be
   *     complemented
   */
  public static Sequence reverseComplement(Sequence seq)
  throws IllegalAlphabetException {
    return new RevCompSequence(seq);
  }

  /**
   * Create a new sequence that has all of the data in the original, but allows
   * new features and top-level annotations to be added independantly. Use this
   * as a scratch-space.
   *
   * @param seq  the Sequence to view
   * @return a new ViewSequence
   */
  public static ViewSequence view(Sequence seq) {
    return new ViewSequence(seq);
  }

  /**
   * Create a new sequence that has all of the data in the original, but allows
   * new features and top-level annotations to be added independantly. Use this
   * as a scratch-space.
   *
   * @param seq  the Sequence to view
   * @param name a new name for the sequence
   * @return a new ViewSequence with the new name
   */
  public static ViewSequence view(Sequence seq, String name) {
    return new ViewSequence(seq, name);
  }

  /**
   * Creates a new Sequence with the data of the old but with a different
   * FeatureRealizer that will be applied to new Features.
   *
   * @param seq the Sequence to wrap
   * @param fr the new FeatureRealizer
   * @return the new ViewSequence
   */
  public static ViewSequence view(Sequence seq, FeatureRealizer fr){
    return new ViewSequence(seq, fr);
  }

  /**
   * Create a new gapped sequence for a sequence.
   *
   * <p>
   * The gapped sequence can be used to insert gaps. The features on the
   * underlying sequence will be projected onto the view taking the gaps into
   * account.
   * </p>
   *
   * @param seq
   * @return a GappedSequence view of seq
   */
  public static GappedSequence gappedView(Sequence seq) {
    return new SimpleGappedSequence(seq);
  }

  /**
   * Mask of a sequence.
   *
   * <P>
   * This will return a view of a sequence where everything outside loc is
   * dropped. This includes all symbols, which become gaps, and all features,
   * which behave in a similar manner to those produced by subSequence().
   * </p>
   *
   * @param seq  the Sequence to mask
   * @param loc  the region to retain
   * @return  a Sequence viewing just the retained portion of seq
   * @throws IndexOutOfBoundsException  if loc is not totaly within seq
   * @throws IllegalArgumentException  fixme: not sure where this comes from
   */
  public static Sequence maskSequence(Sequence seq, RangeLocation loc)
  throws IndexOutOfBoundsException, IllegalArgumentException {
    GappedSequence gSeq = gappedView(subSequence(
            seq,
            loc.getMin(),
            loc.getMax(),
            seq.getName() + ":" + loc.toString()));
    gSeq.addGapsInSource(1, loc.getMin());
    gSeq.addGapsInSource(seq.length(), gSeq.length() - gSeq.length());

    return gSeq;
  }

  /**
   * Create a new Sequence that has no annotation, no features and a zero-length
   * symbol list.
   *
   * Instantiate this if an API requres a sequence, but you can't be bothered
   * or are not able to provide full sequence information.
   * 
   * It is sometimes usefull to create a dummy sequence and then wrap this in
   * a view.
   *
   * @param uri  the URI to give the dummy sequence
   * @param name the name of the dummy sequence
   * @return a dummy Sequence
   */
  public static Sequence createDummy(String uri, String name) {
    return new DummySequence(uri, name);
  }

  /**
   * Create a new Sequence that contains a single symbol repeated over and over.
   *
   * @param alpha   the Alphabet this sequence is over
   * @param length  the length of the sequence
   * @param sym     the symbol returned by every call to symbolAt
   * @param uri     the URI of the sequence
   * @param name    the name of the sequence
   * @return  a new sequence of the right length
   * @throws IllegalSymbolException if sym is not in alpha
   *
   * @since 1.4
   */
  public static Sequence createDummy(
          Alphabet alpha, int length, Symbol sym,
          String uri, String name)
          throws IllegalSymbolException
  {
    return createSequence(new DummySymbolList(alpha, length, sym),
                          uri, name, new SmallAnnotation());
  }

  /**
   * Add features to a sequence that contain the same information as all
   * those in a feature holder.
   *
   * @param seq  the Sequence to add features to
   * @param fh  the features to add
   * @throws ChangeVetoException if the sequence could not be modified
   * @throws BioException if there was an error creating the features
   */
  public static void addAllFeatures(Sequence seq, FeatureHolder fh)
  throws
    ChangeVetoException,
    BioException
  {
    addFeatures(seq, fh);
  }

  private static void addFeatures(FeatureHolder toAddTo, FeatureHolder thingsToAdd)
  throws
    ChangeVetoException,
    BioException
  {
    for(Iterator i = thingsToAdd.features(); i.hasNext(); ) {
      Feature f2add = (Feature) i.next();
      Feature added = toAddTo.createFeature(f2add.makeTemplate());
      addFeatures(added, f2add);
    }
  }
}
