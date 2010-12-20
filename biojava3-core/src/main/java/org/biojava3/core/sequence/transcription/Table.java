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
package org.biojava3.core.sequence.transcription;

import java.util.List;

import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.util.Equals;
import org.biojava3.core.util.Hashcoder;

/**
 * Provides a way of separating us from the specific {@link IUPACTable} even
 * though this is the only implementing class for the interface.
 *
 * @author ayates
 */
public interface Table {

    List<Codon> getCodons(CompoundSet<NucleotideCompound> nucelotides,
            CompoundSet<AminoAcidCompound> aminoAcids);

    CompoundSet<Codon> getCodonCompoundSet(
            final CompoundSet<NucleotideCompound> rnaCompounds,
            final CompoundSet<AminoAcidCompound> aminoAcidCompounds);

    /**
     * Returns true if the given compound could have been a start amino acid;
     * this does not assert if the codon that actually coded for the amino
     * acid was a start codon. This is as accurate a call as we can make with an
     * {@link AminoAcidCompound}.
     */
    boolean isStart(AminoAcidCompound compound);

    /**
     * Instance of a Codon which is 3 {@link NucleotideCompound}s, its
     * corresponding {@link AminoAcidCompound} and if it is a start or stop codon.
     * The object implements hashCode & equals but according to the nucleotide
     * compounds & not to the designation of it being a start, stop & amino
     * acid compound
     *
     * @author ayates
     *
     */
    public static class Codon implements Compound {

        private final CaseInsensitiveTriplet triplet;
        private final boolean start;
        private final boolean stop;
        private final AminoAcidCompound aminoAcid;
        private final String stringified;

        public Codon(CaseInsensitiveTriplet triplet, AminoAcidCompound aminoAcid, boolean start,
                boolean stop) {
            this.triplet = triplet;
            this.start = start;
            this.stop = stop;
            this.aminoAcid = aminoAcid;
            this.stringified = triplet.toString();
        }

        public Codon(CaseInsensitiveTriplet triplet) {
            this(triplet, null, false, false);
        }

        public NucleotideCompound getOne() {
            return triplet.getOne();
        }

        public NucleotideCompound getTwo() {
            return triplet.getTwo();
        }

        public NucleotideCompound getThree() {
            return triplet.getThree();
        }

        public boolean isStart() {
            return start;
        }

        public boolean isStop() {
            return stop;
        }

        public AminoAcidCompound getAminoAcid() {
            return aminoAcid;
        }

        public CaseInsensitiveTriplet getTriplet() {
            return triplet;
        }

        @Override
        public boolean equals(Object obj) {
            boolean equals = false;
            if(Equals.classEqual(this, obj)) {
                Codon casted = (Codon) obj;
                equals =   Equals.equal(getTriplet(), casted.getTriplet()) &&
                            Equals.equal(isStart(), casted.isStart()) &&
                            Equals.equal(isStop(), casted.isStop()) &&
                            Equals.equal(getAminoAcid(), casted.getAminoAcid());
            }
            return equals;
        }

        @Override
        public int hashCode() {
            int result = Hashcoder.SEED;
            result = Hashcoder.hash(result, getTriplet());
            result = Hashcoder.hash(result, isStop());
            result = Hashcoder.hash(result, isStart());
            result = Hashcoder.hash(result, getAminoAcid());
            return result;
        }

        @Override
        public String toString() {
            return stringified;
        }

        @Override
        public boolean equalsIgnoreCase(Compound compound) {
            return toString().equalsIgnoreCase(compound.toString());
        }

        @Override
        public String getDescription() {
            throw new UnsupportedOperationException("Not supported");
        }

        @Override
        public String getLongName() {
            throw new UnsupportedOperationException("Not supported");
        }

        @Override
        public Float getMolecularWeight() {
            throw new UnsupportedOperationException("Not supported");
        }

        @Override
        public String getShortName() {
            return stringified;
        }

        @Override
        public void setDescription(String description) {
            throw new UnsupportedOperationException("Not supported");
        }

        @Override
        public void setLongName(String longName) {
            throw new UnsupportedOperationException("Not supported");
        }

        @Override
        public void setMolecularWeight(Float molecularWeight) {
            throw new UnsupportedOperationException("Not supported");
        }

        @Override
        public void setShortName(String shortName) {
            throw new UnsupportedOperationException("Not supported");
        }
    }

    /**
     * Class used to hold three nucleotides together and allow for equality
     * to be assessed in a case insensitive manner.
     */
    public static class CaseInsensitiveTriplet {

        private final NucleotideCompound one;
        private final NucleotideCompound two;
        private final NucleotideCompound three;

        private transient boolean hashSet = false;
        private transient int hash;
        private transient boolean stringSet = false;
        private transient String stringify;

        public CaseInsensitiveTriplet(NucleotideCompound one,
                NucleotideCompound two, NucleotideCompound three) {
            this.one = one;
            this.two = two;
            this.three = three;

        }

        public NucleotideCompound getOne() {
            return one;
        }

        public NucleotideCompound getTwo() {
            return two;
        }

        public NucleotideCompound getThree() {
            return three;
        }

        @Override
        public boolean equals(Object obj) {
            boolean equals = false;
            if(Equals.classEqual(this, obj)) {
                CaseInsensitiveTriplet casted = (CaseInsensitiveTriplet) obj;
                return toString().equals(casted.toString());
            }
            return equals;
        }

        @Override
        public int hashCode() {
            if(!hashSet) {
                hash = toString().hashCode();
                hashSet = true;
            }
            return hash;
        }

        @Override
        public String toString() {
            if(!stringSet) {
                stringify = getOne().getUpperedBase() +
                    getTwo().getUpperedBase() +
                    getThree().getUpperedBase();
            }
            return stringify;
        }

        /**
         * Attempts to provide an int version of this codon which multiplies
         * each position by
         */
        public int intValue() {
            return (16 * compoundToInt(getOne())) +
                    (4 * compoundToInt(getTwo())) +
                    (compoundToInt(getThree()));
        }

        public int compoundToInt(NucleotideCompound c) {
            char b = c.getUpperedBase().charAt(0);
            return (int)b;
//            int v = -1;
//            if('A' == b) {
//                v = 1;
//            }
//            else if('C' == b) {
//                v = 2;
//            }
//            else if('G' == b) {
//                v = 3;
//            }
//            else if('T' == b || 'U' == b) {
//                v = 4;
//            }
//            return v;
        }
    }
}
