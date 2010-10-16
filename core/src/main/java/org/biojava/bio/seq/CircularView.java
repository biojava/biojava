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

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.impl.ViewSequence;
import org.biojava.bio.symbol.CircularLocation;
import org.biojava.bio.symbol.Edit;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeVetoException;

/**
 * <p>
 * A circular view onto another Sequence object.  The class allows for
 * reinterpretation of locations and indices onto the sequence to allow for
 * overlapping of the origin. The origin is assumed to be the first symbol.
 * Future versions may support changing the origin.
 * </p>
 * In biojavax {@link org.biojavax.bio.seq.RichSequence RichSequences} intrinsicly know about
 * circularity. No view is required. We strongly recommend using RichSequence
 * if possible.
 *
 *
 * @author Mark Schreiber
 * @version 1.2
 * @since 1.1
 * @see org.biojavax.bio.seq.RichSequence 
 */

public class CircularView extends ViewSequence{
  public CircularView(Sequence seq, FeatureRealizer fr){
    super(seq, fr);
  }

  public CircularView(Sequence seq){
    super(seq);
  }

  private int realValue(int val){
    val = (val % length());
    if(val < 1) val = length() + val;
    return val;
  }

  /**
   * <p>
   * Over rides ViewSequence. Allows any integer index, positive or negative
   * to return a symbol via the equation
   * <CODE>val = (val % length);</CODE>
   * Note that base zero is the base immediately before base 1 which is of course
   * the last base of the sequence.
   * </p>
   * @param index the index of the <code>Symbol </code>requested.
   * @return the <code>Symbol </code>specified by the <code>index.</code>
   */
  public Symbol symbolAt(int index){
    return super.symbolAt(realValue(index));
  }

  /**
   * <p>
   * Over rides ViewSequence. Allows any integer index, positive or negative
   * to return a symbol via the equation
   * <CODE>val = (val % length);</CODE>
   * </p>
   *
   * <p>
   * Will return a linear String which can, if nescessary, span the origin.
   * </p>
   * @param start the index of the fist base
   * @param end the index of the last base
   * @return a <code>String </code>representation of the tokenized <code>Symbol</code>s
   */
  public String subStr(int start, int end){

    start = realValue(start);
    end = realValue(end);

    if(start <= end){
      return super.subStr(start, end);
    }
    else{
      String toEnd = super.subStr(start,super.length());
      String fromStart = super.subStr(1,end);
      return toEnd + fromStart;
    }
  }

  /**
   * Over rides ViewSequence to allow the use of locations that have
   * coordinates outside of the sequence length (which are needed to
   * describe locations that overlap the origin of a circular sequence).
   *
   * @since 1.2
   * @throws BioException if a non circular location is added that exceeds the
   * 'boundaries' of the sequence.
   * @throws ChangeVetoException if the sequence is locked.
   * @param template the template of the feature to be created.
   * @return the feature created you can use the template of the returned feature
   * to create another of the same type.
   */
  public Feature createFeature(Feature.Template template)
        throws ChangeVetoException, BioException
    {
      Location loc = template.location;
      if(loc.getMax() > length() && (loc instanceof CircularLocation == false)){
        throw new BioException("Only CircularLocations may exceed sequence length");
      }
      Feature f = realizeFeature(this, template);
      ((SimpleFeatureHolder)getAddedFeatures()).addFeature(f);
      return f;
    }

  /**
   * <p>
   * Over rides ViewSequence. Allows any integer index, positive or negative
   * to return a symbol via the equation
   * <CODE>index = ((index -1) % length)+1</CODE>
   * </p>
   *
   * <p>
   * Will return a linear SymbolList which can ,if nescessary, span the origin.
   * </p>
   * @param start the first base of the sublist
   * @param end the last base of the sublist
   * @return a <code>SymbolList </code>containing the <code>Symbols</code> from
   * <code>start</code> to <code>end</code> inclusive
   */
  public SymbolList subList(int start, int end){

    start = realValue(start);
    end = realValue(end);

    if(start <= end){
      return super.subList(start, end);
    }
    else{
      SimpleSymbolList fromStart = new SimpleSymbolList(super.subList(1,end));
      SimpleSymbolList toEnd = new SimpleSymbolList(super.subList(start,length()));
      Edit edit = new Edit(toEnd.length() +1, 0, fromStart);
      try{
        toEnd.edit(edit);
      }catch(Exception e){
        throw new BioError("Couldn't construct subList, this shouldn't happen",e);
      }
      return toEnd;
    }
  }
}

