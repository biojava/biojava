package org.biojava3.core.sequence.template;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.biojava3.core.sequence.compound.NucleotideCompound;

/**
 * Provides a set of static methods to be used as static imports when needed
 * across multiple Sequence implementations but inheritance gets in the way.
 *
 * @author ayates
 */
public class SequenceMixin {

  /**
   * For the given vargs of compounds this method counts the number of
   * times those compounds appear in the given sequence
   *
   * @param sequence The {@link Sequence} to perform the count on
   * @param compounds The compounds to look for
   * @param <C> The type of compound we are looking for
   * @return The number of times the given compounds appear in this Sequence
   */
  public static <C extends Compound> int countCompounds(
      Sequence<C> sequence, C... compounds) {
    int count = 0;
    Set<C> compoundSet = new HashSet<C>(Arrays.asList(compounds));
    for (C currentCompound : sequence) {
      if(compoundSet.contains(currentCompound)) {
        count++;
      }
    }
    return count;
  }

  /**
   * Returns the count of GC in the given sequence
   *
   * @param sequence The {@link NucleotideCompound} {@link Sequence} to perform
   * the GC analysis on
   * @return The number of GC compounds in the sequence
   */
  public static int countGC(Sequence<NucleotideCompound> sequence) {
    CompoundSet<NucleotideCompound> cs = sequence.getCompoundSet();
    NucleotideCompound G = cs.getCompoundForString("G");
    NucleotideCompound C = cs.getCompoundForString("C");
    NucleotideCompound g = cs.getCompoundForString("g");
    NucleotideCompound c = cs.getCompoundForString("c");
    return countCompounds(sequence, G, C, g, c);
  }

  public static <C extends Compound> StringBuilder toStringBuilder(Sequence<C> sequence) {
    StringBuilder sb = new StringBuilder(sequence.getLength());
    for(C compound: sequence) {
      sb.append(compound.toString());
    }
    return sb;
  }

  public static <C extends Compound> List<C> toList(Sequence<C> sequence) {
    List<C> list = new ArrayList<C>(sequence.getLength());
    for(C compound: sequence) {
      list.add(compound);
    }
    return list;
  }

  public static <C extends Compound> int indexOf(Sequence<C> sequence, C compound) {
    int index = 1;
    for(C currentCompound: sequence) {
      if(currentCompound.equals(compound)) {
        return index;
      }
      index++;
    }
    return 0;
  }

  public static <C extends Compound> int lastIndexOf(Sequence<C> sequence, C compound) {
    for(int index = sequence.getLength(); index >= 1; index--) {
      C currentCompound = sequence.getCompoundAt(index);
      if(currentCompound.equals(compound)) {
        return index;
      }
    }
    return 0;
  }
}
