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
 * Created on 01-21-2010
 */
package org.biojava3.core.sequence.template;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.biojava3.core.exceptions.CompoundNotFoundError;
import org.biojava3.core.util.Equals;
import org.biojava3.core.util.Hashcoder;


/**
 *
 * @author Andy Yates
 *
 * @param <C> The compound this set will contain
 */
public abstract class AbstractCompoundSet<C extends Compound> implements CompoundSet<C> {

  private Map<CharSequence, C> charSeqToCompound = new HashMap<CharSequence, C>();
  private int maxCompoundCharSequenceLength = -1;
  private Boolean compoundStringLengthEqual = null;
  
  Map<C,Set<C>> equivalentsMap = new HashMap<C, Set<C>>();

  protected void addCompound(C compound, C lowerCasedCompound, Iterable<C> equivalents) {
    addCompound(compound);
    addCompound(lowerCasedCompound);

    addEquivalent(compound, lowerCasedCompound);
    addEquivalent(lowerCasedCompound, compound);

    for(C equivalent: equivalents) {
      addEquivalent(compound, equivalent);
      addEquivalent(equivalent, compound);
      addEquivalent(lowerCasedCompound, equivalent);
      addEquivalent(equivalent, lowerCasedCompound);
    }
  }

  protected void addCompound(C compound, C lowerCasedCompound, C... equivalents) {
    List<C> equiv = new ArrayList<C>(equivalents.length);
    equiv.addAll(Arrays.asList(equivalents));
    addCompound(compound, lowerCasedCompound, equiv);
  }

  protected void addEquivalent(C compound, C equivalent) {
	 Set<C> s = equivalentsMap.get(compound);
	 if ( s == null){
		 s = new HashSet<C>();
		 equivalentsMap.put(compound, s);
	 }
	  
    s.add( equivalent);
  }

  protected void addCompound(C compound) {
    charSeqToCompound.put(compound.toString(), compound);
    maxCompoundCharSequenceLength = -1;
    compoundStringLengthEqual = null;
  }

  public String getStringForCompound(C compound) {
    return compound.toString();
  }

  public C getCompoundForString(String string) {
    if(string == null) {
      throw new IllegalArgumentException("Given a null CharSequence to process");
    }

    if (string.length()==0) {
      return null;
    }

    if (string.length() > getMaxSingleCompoundStringLength()) {
      throw new IllegalArgumentException("CharSequence supplied is too long.");
    }

    return charSeqToCompound.get(string);
  }

  public int getMaxSingleCompoundStringLength() {
    if(maxCompoundCharSequenceLength == -1) {
      for(C compound: charSeqToCompound.values()) {
        int size = getStringForCompound(compound).length();
        if(size > maxCompoundCharSequenceLength) {
          maxCompoundCharSequenceLength = size;
        }
      }
    }
    return maxCompoundCharSequenceLength;
  }

    @Override
    public boolean isCompoundStringLengthEqual() {
        if(compoundStringLengthEqual == null) {
            int lastSize = -1;
            compoundStringLengthEqual = Boolean.TRUE;
            for(CharSequence c: charSeqToCompound.keySet()) {
                int size = c.length();
                if(lastSize != -1) {
                    lastSize = size;
                    continue;
                }
                if(lastSize != size) {
                    compoundStringLengthEqual = Boolean.FALSE;
                    break;
                }
            }
        }
        return compoundStringLengthEqual;
    }

  public boolean hasCompound(C compound) {
    C retrievedCompound = getCompoundForString(compound.toString());
    return (retrievedCompound == null) ? false : true;
  }

  public boolean compoundsEquivalent(C compoundOne, C compoundTwo) {
    assertCompound(compoundOne);
    assertCompound(compoundTwo);
    return equivalentsMap.get(compoundOne).contains(compoundTwo);
  }

  public Set<C> getEquivalentCompounds(C compound) {
    return equivalentsMap.get(compound);
  }

  public boolean compoundsEqual(C compoundOne, C compoundTwo) {
    assertCompound(compoundOne);
    assertCompound(compoundTwo);
    return compoundOne.equalsIgnoreCase(compoundTwo);
  }

  public void verifySequence(Sequence<C> sequence) throws CompoundNotFoundError {
    for(C compound: sequence) {
      assertCompound(compound);
    }
  }

  public List<C> getAllCompounds() {
    return new ArrayList<C>(charSeqToCompound.values());
  }

  private void assertCompound(C compound) {
    boolean okay = hasCompound(compound);
    if(! okay) {
      throw new CompoundNotFoundError("The CompoundSet "+
          getClass().getSimpleName()+" knows nothing about the compound "+
          compound);
    }
  }

    @Override
    public boolean isComplementable() {
        return false;
    }

    @Override
    public int hashCode() {
        int s = Hashcoder.SEED;
        s = Hashcoder.hash(s, charSeqToCompound);
        s = Hashcoder.hash(s, equivalentsMap);
        return s;
    }

    @Override
    @SuppressWarnings("unchecked")
    public boolean equals(Object o) {
        if(Equals.classEqual(this, o)) {
            AbstractCompoundSet<C> that = (AbstractCompoundSet<C>)o;
            return  Equals.equal(charSeqToCompound, that.charSeqToCompound) &&
                    Equals.equal(equivalentsMap, that.equivalentsMap);
        }
        return false;
    }


}
