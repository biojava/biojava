/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.genome.parsers.gff;

import java.util.Collection;
import org.biojava3.core.sequence.DNASequence;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class GCStats {

    public static double getGCStats(Collection<DNASequence> sequences) {
        double gcCount = 0;
        double total = 0;

        for (DNASequence sequence : sequences) {
            char[] dna = sequence.toString().toCharArray();
            for (char d : dna) {
                if (d == 'G' || d == 'C' || d == 'g' || d == 'c') {
                    gcCount++;
                }
                total++;
            }
        }

        return (gcCount / total) * 100.0;
    }

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
}
