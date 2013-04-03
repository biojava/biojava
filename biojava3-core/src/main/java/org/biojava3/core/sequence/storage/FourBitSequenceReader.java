/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.core.sequence.storage;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.Sequence;

/**
 *
 * Four bit encoding of the bit formats. This can support up to 16 compounds
 * from a compound set. To allow us to support the redundant set of Nucleotide
 * compounds this class will use case-insensitive encoding. The values assigned
 * to these compounds is also done at runtime; if you want a predictable
 * ordering then override and use your own encodings. However all
 * encodings are calculated using lexographical ordering of the compounds
 * so if a CompoundSet does not change then this encoding should not cauuse
 * a problem.
 *
 * @author ayates
 */
public class FourBitSequenceReader<C extends Compound> extends BitSequenceReader<C> {

    public FourBitSequenceReader(Sequence<C> sequence) {
        super(new FourBitArrayWorker<C>(sequence), sequence.getAccession());
    }

    public FourBitSequenceReader(String sequence, CompoundSet<C> compoundSet) {
        this(sequence, compoundSet, new AccessionID("Unknown"));
    }

    public FourBitSequenceReader(String sequence, CompoundSet<C> compoundSet, AccessionID accession) {
        super(new FourBitArrayWorker<C>(sequence, compoundSet), accession);
    }

    public FourBitSequenceReader(FourBitArrayWorker<C> worker) {
        super(worker, new AccessionID("unknown"));
    }

    public FourBitSequenceReader(FourBitArrayWorker<C> worker, AccessionID accession) {
        super(worker, accession);
    }

    /**
     * A four bit per compound implementation of the bit array worker code. This
     * version can handle upto 16 compounds but this does mean that its ability
     * to compress a normal sequence is halved (compared to the 1/4 performance
     * seen with the 2bit workers).
     *
     * @param <C> Must extend NucleotideCompound
     */
    public static class FourBitArrayWorker<C extends Compound> extends BitArrayWorker<C> {

        public FourBitArrayWorker(CompoundSet<C> compoundSet, int length) {
            super(compoundSet, length);
        }

        public FourBitArrayWorker(CompoundSet<C> compoundSet, int[] sequence) {
            super(compoundSet, sequence);
        }

        public FourBitArrayWorker(Sequence<C> sequence) {
            super(sequence);
        }

        public FourBitArrayWorker(String sequence, CompoundSet<C> compoundSet) {
            super(sequence, compoundSet);
        }
        /**
         * Masking value used for extracting the right most 2 bits from a byte
         */
        private final static byte MASK = (byte) ((int) Math.pow(2, 0) | (int) Math.pow(2, 1) | (int) Math.pow(2, 2) | (int) Math.pow(2, 3));

        
        protected byte bitMask() {
            return MASK;
        }

        
        protected int compoundsPerDatatype() {
            return 8;
        }

        /**
         * Returns a Map which encodes the contents of CompoundSet. This
         * version is case-insensitive i.e. C and c both encode for the same
         * position. We sort lexigraphically so if the compound set has
         * not changed then neither will this.
         */
        
        protected Map<C, Integer> generateCompoundsToIndex() {
            final CompoundSet<C> cs = getCompoundSet();
            Map<C, Integer> map = new HashMap<C, Integer>();
            int index = 0;
            for (C currentCompound : sortedCompounds(cs)) {
                C upperCasedCompound = getOptionalUpperCasedCompound(currentCompound, cs);

                //if it has the uppercased compound then set this
                //compounds' value to that one
                if (map.containsKey(upperCasedCompound)) {
                    map.put(currentCompound, map.get(upperCasedCompound));
                } else {
                    map.put(currentCompound, index++);
                }
            }

            return map;
        }

        private C getOptionalUpperCasedCompound(C currentCompound, CompoundSet<C> cs) {
            C upperCasedCompound = null;
            String upperCasedString = cs.getStringForCompound(currentCompound).toUpperCase();
            if (cs.getCompoundForString(upperCasedString) == null) {
                upperCasedCompound = currentCompound;
            } else {
                upperCasedCompound = cs.getCompoundForString(upperCasedString);
            }
            return upperCasedCompound;
        }

        private List<C> sortedCompounds(final CompoundSet<C> cs) {
            List<C> compounds = new ArrayList<C>(cs.getAllCompounds());
            Collections.sort(compounds, new Comparator<C>() {

                
                public int compare(C o1, C o2) {
                    String s1 = cs.getStringForCompound(o1);
                    String s2 = cs.getStringForCompound(o2);
                    return s1.compareTo(s2);
                }
            });
            return compounds;
        }

        /**
         * Returns a List which reverse encodes the Compound, Integer map
         */
        
        protected List<C> generateIndexToCompounds() {
            CompoundSet<C> cs = getCompoundSet();
            Map<C, Integer> lookup = getCompoundsToIndexLookup();
            Map<Integer, C> tempMap = new HashMap<Integer, C>();
            //First get the reverse lookup working
            for (C compound : lookup.keySet()) {
                C upperCasedCompound = getOptionalUpperCasedCompound(compound, cs);
                Integer pos = lookup.get(upperCasedCompound);
                tempMap.put(pos, upperCasedCompound);
            }

            //Then populate the results by going back through the sorted integer keys
            List<C> compounds = new ArrayList<C>();
            List<Integer> keys = new ArrayList<Integer>(tempMap.keySet());
            Collections.sort(keys);
            for (Integer key : keys) {
                compounds.add(tempMap.get(key));
            }

            return compounds;
        }
    }
}
