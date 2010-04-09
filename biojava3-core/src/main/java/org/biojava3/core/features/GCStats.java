/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.core.features;

import java.util.Collection;

import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.template.Sequence;
import org.biojava3.core.sequence.template.SequenceMixin;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class GCStats {

  /**
   * Calculates the total percentage of GC found in the given sequences.
   *
   * @param sequences List of Sequences to do GC count on
   * @return The percentage of GC in the Sequences
   */
  public static double getGCStats(Collection<Sequence<NucleotideCompound>> sequences) {
    double gcCount = 0;
    double total = 0;
    for (Sequence<NucleotideCompound> sequence : sequences) {
      gcCount += SequenceMixin.countGC(sequence);
      total += sequence.getLength();
    }
    return (gcCount / total) * 100.0;
  }

  /**
   * Calculates the total percentage of GC found in the given sequences.
   *
   * @param sequences List of Strings to do GC count on
   * @return The percentage of GC in the strings
   */
  public static double getGCStatsString(Collection<String> sequences) {
    double gcCount = 0;
    double total = 0;
    for (String sequence : sequences) {
      char[] dna = sequence.toCharArray();
      for (char d : dna) {
        if (d == 'G' || d == 'C' || d == 'g' || d == 'c') {
          gcCount++;
        }
        total++;
      }
    }
    return (gcCount / total) * 100.0;
  }

  /**
   * Calculates the percentage of GC in the given sequence
   *
   * @param sequence The {@link NucleotideCompound} {@link Sequence} to perform
   * the GC analysis on
   * @return The percentage of GC in the sequence
   */
  public static double getGCStats(Sequence<NucleotideCompound> sequence) {
    double count = SequenceMixin.countGC(sequence);
    double size = sequence.getLength();
    return (count / size) * 100.0;
  }
}
