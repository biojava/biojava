package org.biojava3.core.sequence.template;

import java.util.Arrays;
import java.util.HashSet;
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
   * @param <T> The type of compound we are looking for
   * @return The number of times the given compounds appear in this Sequence
   */
  public static <T extends Compound> int countCompounds(
      Sequence<T> sequence, T... compounds) {
    int count = 0;
    Set<T> compoundSet = new HashSet<T>(Arrays.asList(compounds));
    for (T currentCompound : sequence) {
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
}
