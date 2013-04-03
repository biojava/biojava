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

import java.util.EnumMap;
import java.util.Map;

import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.compound.DNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.compound.RNACompoundSet;
import org.biojava3.core.sequence.io.IUPACParser;
import org.biojava3.core.sequence.io.ProteinSequenceCreator;
import org.biojava3.core.sequence.io.RNASequenceCreator;
import org.biojava3.core.sequence.io.IUPACParser.IUPACTable;
import org.biojava3.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.Sequence;
import org.biojava3.core.sequence.transcription.Table.Codon;

/**
 * Used as a way of encapsulating the data structures required to parse
 * DNA to a Protein sequence.
 *
 * In order to build one look at {@ TranscriptionEngine.Builder} which provides
 * intelligent defaults & allows you to build an engine which is nearly
 * the same as the default one but with a few changes. All of the engine
 * is customisable.
 *
 * By default the code will attempt to:
 *
 * <ul>
 * <li>Trim Stops</li>
 * <li>Convert initiating codons to M</li>
 * <li>Allow for the fuzzy translation of Codons i.e. if it contains an
 * N that produces a {@link Sequence}&lt;{@link{AminoAcidCompound}&gt;
 * with an X at that position
 * </ul>
 *
 * @author ayates
 */
public class TranscriptionEngine {

    private static final class IOD {

        public static final TranscriptionEngine INSTANCE = new TranscriptionEngine.Builder().build();
    }

    /**
     * Default instance to use when Transcribing from DNA -&gt; RNA -&gt;
     * Protein. If you require anything that is not a default setting then look
     * at {@ TranscriptionEngine.Builder} for customisation options.
     */
    public static TranscriptionEngine getDefault() {
        return IOD.INSTANCE;
    }
    private final Table table;
    private final RNAToAminoAcidTranslator rnaAminoAcidTranslator;
    private final DNAToRNATranslator dnaRnaTranslator;
    private final SequenceCreatorInterface<AminoAcidCompound> proteinSequenceCreator;
    private final SequenceCreatorInterface<NucleotideCompound> rnaSequenceCreator;
    private final CompoundSet<NucleotideCompound> dnaCompounds;
    private final CompoundSet<NucleotideCompound> rnaCompounds;
    private final CompoundSet<AminoAcidCompound> aminoAcidCompounds;

    private TranscriptionEngine(
            Table table,
            RNAToAminoAcidTranslator rnaAminoAcidTranslator,
            DNAToRNATranslator dnaRnaTranslator,
            SequenceCreatorInterface<AminoAcidCompound> proteinSequenceCreator,
            SequenceCreatorInterface<NucleotideCompound> rnaSequenceCreator,
            CompoundSet<NucleotideCompound> dnaCompounds,
            CompoundSet<NucleotideCompound> rnaCompounds,
            CompoundSet<AminoAcidCompound> aminoAcidCompounds) {
        this.table = table;
        this.rnaAminoAcidTranslator = rnaAminoAcidTranslator;
        this.dnaRnaTranslator = dnaRnaTranslator;
        this.proteinSequenceCreator = proteinSequenceCreator;
        this.rnaSequenceCreator = rnaSequenceCreator;
        this.dnaCompounds = dnaCompounds;
        this.rnaCompounds = rnaCompounds;
        this.aminoAcidCompounds = aminoAcidCompounds;
    }

    /**
     * Quick method to let you go from a CDS to a Peptide quickly. It assumes
     * you are translating only in the first frame
     *
     * @param dna The CDS to translate
     * @return The Protein Sequence
     */
    public Sequence<AminoAcidCompound> translate(Sequence<NucleotideCompound> dna) {
        Map<Frame, Sequence<AminoAcidCompound>> trans =
                multipleFrameTranslation(dna, Frame.ONE);
        return trans.get(Frame.ONE);
    }

    /**
     * A way of translating DNA in a number of frames
     *
     * @param dna The CDS to translate
     * @param frames The Frames to translate in
     * @return All generated protein sequences in the given frames. Can have
     * null entries
     */
    public Map<Frame, Sequence<AminoAcidCompound>> multipleFrameTranslation(
            Sequence<NucleotideCompound> dna, Frame... frames) {
        Map<Frame, Sequence<AminoAcidCompound>> results =
                new EnumMap<Frame, Sequence<AminoAcidCompound>>(Frame.class);
        for (Frame frame : frames) {
            Sequence<NucleotideCompound> rna =
                    getDnaRnaTranslator().createSequence(dna, frame);
            Sequence<AminoAcidCompound> peptide =
                    getRnaAminoAcidTranslator().createSequence(rna);
            results.put(frame, peptide);
        }
        return results;
    }

    public Table getTable() {
        return table;
    }

    public RNAToAminoAcidTranslator getRnaAminoAcidTranslator() {
        return rnaAminoAcidTranslator;
    }

    public DNAToRNATranslator getDnaRnaTranslator() {
        return dnaRnaTranslator;
    }

    public SequenceCreatorInterface<AminoAcidCompound> getProteinSequenceCreator() {
        return proteinSequenceCreator;
    }

    public SequenceCreatorInterface<NucleotideCompound> getRnaSequenceCreator() {
        return rnaSequenceCreator;
    }

    public CompoundSet<NucleotideCompound> getDnaCompounds() {
        return dnaCompounds;
    }

    public CompoundSet<NucleotideCompound> getRnaCompounds() {
        return rnaCompounds;
    }

    public CompoundSet<AminoAcidCompound> getAminoAcidCompounds() {
        return aminoAcidCompounds;
    }

    /**
     * This class is the way to create a {@link TranslationEngine}.
     */
    public static class Builder {

        private Table table;
        private RNAToAminoAcidTranslator rnaAminoAcidTranslator;
        private DNAToRNATranslator dnaRnaTranslator;
        private SequenceCreatorInterface<AminoAcidCompound> proteinSequenceCreator;
        private SequenceCreatorInterface<NucleotideCompound> rnaSequenceCreator;
        private CompoundSet<NucleotideCompound> dnaCompounds;
        private CompoundSet<NucleotideCompound> rnaCompounds;
        private CompoundSet<AminoAcidCompound> aminoAcidCompounds;
        private boolean initMet = true;
        private boolean trimStop = true;
        private boolean translateNCodons = true;
        private boolean decorateRna = false;

        /**
         * The method to finish any calls to the builder with which returns
         * a transcription engine. The engine is designed to provide everything
         * required for transcription to those classes which will do the
         * transcription.
         */
        public TranscriptionEngine build() {
            return new TranscriptionEngine(
                    getTable(),
                    getRnaAminoAcidTranslator(),
                    getDnaRnaTranslator(),
                    getProteinCreator(),
                    getRnaCreator(),
                    getDnaCompounds(),
                    getRnaCompounds(),
                    getAminoAcidCompounds());
        }

        //---- START OF BUILDER METHODS
        /**
         * Uses the static instance of {@link IUPACParser} to find instances of
         * {@link IUPACTable}s by ID.
         */
        public Builder table(Integer id) {
            table = IUPACParser.getInstance().getTable(id);
            return this;
        }

        /**
         * Uses the static instance of {@link IUPACParser} to find instances of
         * {@link IUPACTable}s by its String name
         */
        public Builder table(String name) {
            table = IUPACParser.getInstance().getTable(name);
            return this;
        }

        public Builder table(Table table) {
            this.table = table;
            return this;
        }

        public Builder dnaCompounds(CompoundSet<NucleotideCompound> compounds) {
            this.dnaCompounds = compounds;
            return this;
        }

        public Builder rnaCompounds(CompoundSet<NucleotideCompound> compounds) {
            this.rnaCompounds = compounds;
            return this;
        }

        public Builder aminoAcidsCompounds(CompoundSet<AminoAcidCompound> compounds) {
            this.aminoAcidCompounds = compounds;
            return this;
        }

        public Builder dnaRnaTranslator(DNAToRNATranslator translator) {
            this.dnaRnaTranslator = translator;
            return this;
        }

        public Builder rnaAminoAcidTranslator(RNAToAminoAcidTranslator translator) {
            this.rnaAminoAcidTranslator = translator;
            return this;
        }

        public Builder proteinCreator(SequenceCreatorInterface<AminoAcidCompound> creator) {
            this.proteinSequenceCreator = creator;
            return this;
        }

        public Builder rnaCreator(SequenceCreatorInterface<NucleotideCompound> creator) {
            this.rnaSequenceCreator = creator;
            return this;
        }

        public Builder initMet(boolean initMet) {
            this.initMet = initMet;
            return this;
        }

        public Builder trimStop(boolean trimStop) {
            this.trimStop = trimStop;
            return this;
        }

        public Builder translateNCodons(boolean translateNCodons) {
            this.translateNCodons = translateNCodons;
            return this;
        }

        /**
         * Performs an optimisation where RNASequences are not translated into
         * their own objects but are views onto the base DNA sequence.
         */
        public Builder decorateRna(boolean decorateRna) {
            this.decorateRna = decorateRna;
            return this;
        }

        //------ INTERNAL BUILDERS with defaults if exists
        private CompoundSet<NucleotideCompound> getDnaCompounds() {
            if (dnaCompounds != null) {
                return dnaCompounds;
            }
            return DNACompoundSet.getDNACompoundSet();
        }

        private CompoundSet<NucleotideCompound> getRnaCompounds() {
            if (rnaCompounds != null) {
                return rnaCompounds;
            }
            return RNACompoundSet.getRNACompoundSet();
        }

        private CompoundSet<AminoAcidCompound> getAminoAcidCompounds() {
            if (aminoAcidCompounds != null) {
                return aminoAcidCompounds;
            }
            return AminoAcidCompoundSet.getAminoAcidCompoundSet();
        }

        private DNAToRNATranslator getDnaRnaTranslator() {
            if (dnaRnaTranslator != null) {
                return dnaRnaTranslator;
            }
            return new DNAToRNATranslator(new RNASequenceCreator(getRnaCompounds()),
                    getDnaCompounds(), getRnaCompounds(), isDecorateRna());
        }

        private RNAToAminoAcidTranslator getRnaAminoAcidTranslator() {
            if (rnaAminoAcidTranslator != null) {
                return rnaAminoAcidTranslator;
            }
            return new RNAToAminoAcidTranslator(
                    getProteinCreator(), getRnaCompounds(), getCodons(),
                    getAminoAcidCompounds(), getTable(), isTrimStop(), isInitMet(), isTranslateNCodons());
        }

        private CompoundSet<Codon> getCodons() {
            return getTable().getCodonCompoundSet(getRnaCompounds(), getAminoAcidCompounds());
        }

        private SequenceCreatorInterface<AminoAcidCompound> getProteinCreator() {
            if (proteinSequenceCreator != null) {
                return proteinSequenceCreator;
            }
            return new ProteinSequenceCreator(getAminoAcidCompounds());
        }

        private SequenceCreatorInterface<NucleotideCompound> getRnaCreator() {
            if (rnaSequenceCreator != null) {
                return rnaSequenceCreator;
            }
            return new RNASequenceCreator(getRnaCompounds());
        }

        private Table getTable() {
            if (table != null) {
                return table;
            }
            table(1); //Will set table to default IUPAC codee
            return table;
        }

        private boolean isInitMet() {
            return initMet;
        }

        private boolean isTrimStop() {
            return trimStop;
        }
        private boolean isTranslateNCodons() {
            return translateNCodons;
        }
        private boolean isDecorateRna() {
            return decorateRna;
        }
    }
}
