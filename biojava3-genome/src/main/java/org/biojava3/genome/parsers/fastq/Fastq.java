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
package org.biojava3.genome.parsers.fastq;

import java.util.ArrayList;
import java.util.List;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.features.FeatureInterface;
import org.biojava3.core.sequence.features.QualityFeature;
import org.biojava3.core.sequence.template.AbstractSequence;

/**
 * FASTQ formatted sequence.
 *
 * @since 3.0.3
 */
public final class Fastq {

    /**
     * Description of this FASTQ formatted sequence.
     */
    private final String description;
    /**
     * Sequence for this FASTQ formatted sequence.
     */
    private final String sequence;
    /**
     * Quality scores for this FASTQ formatted sequence.
     */
    private final String quality;
    /**
     * FASTQ sequence format variant for this FASTQ formatted sequence.
     */
    private final FastqVariant variant;

    /**
     * Create a new FASTQ formatted sequence from the specified description,
     * sequence, quality scores, and sequence format variant.
     *
     * @param description description of this FASTQ formatted sequence, must not
     * be null
     * @param sequence sequence for this FASTQ formatted sequence, must not be
     * null
     * @param quality quality scores for this FASTQ formatted sequence, must not
     * be null
     * @param variant FASTQ sequence format variant for this FASTQ formatted
     * sequence, must not be null
     */
    Fastq(final String description, final String sequence, final String quality,
            final FastqVariant variant) {
        if (description == null) {
            throw new IllegalArgumentException("description must not be null");
        }
        if (sequence == null) {
            throw new IllegalArgumentException("sequence must not be null");
        }
        if (quality == null) {
            throw new IllegalArgumentException("quality must not be null");
        }
        if (variant == null) {
            throw new IllegalArgumentException("variant must not be null");
        }
        this.description = description;
        this.sequence = sequence;
        this.quality = quality;
        this.variant = variant;
    }

    Fastq(final DNASequence sequence, final FastqVariant variant) {
        if (sequence == null) {
            throw new IllegalArgumentException("sequence must not be null");
        }

        this.description = sequence.getOriginalHeader();
        this.sequence = sequence.getSequenceAsString();
        this.variant = variant;
        
        List<FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound>> features = sequence.getFeaturesByType("quality");
        
        FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound> qualfeature = null;
        for (FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound> feature : features) {
            if (feature instanceof QualityFeature) {
                //TODO better implementation for multiple quality features per DNASequence
                qualfeature = feature;
                break;
            }
        }
        if (qualfeature != null) {
            this.quality = encodeFeature(qualfeature, variant);
        } else {
            throw new IllegalArgumentException("sequence must have a quality feature");
        }
    }

    /**
     * Return the description of this FASTQ formatted sequence. The description
     * will not be null.
     *
     * @return the description of this FASTQ formatted sequence
     */
    public String getDescription() {
        return description;
    }

    /**
     * Return the sequence for this FASTQ formatted sequence. The sequence will
     * not be null.
     *
     * @return the sequence for this FASTQ formatted sequence
     */
    public String getSequence() {
        return sequence;
    }

    /**
     * Return the quality scores for this FASTQ formatted sequence. The quality
     * scores will not be null.
     *
     * @return the quality scores for this FASTQ formatted sequence
     */
    public String getQuality() {
        return quality;
    }

    /**
     * Return the FASTQ sequence format variant for this FASTQ formatted
     * sequence. The FASTQ sequence format variant will not be null.
     *
     * @return the FASTQ sequence format variant for this FASTQ formatted
     * sequence
     */
    public FastqVariant getVariant() {
        return variant;
    }

    /**
     * Return a biojava DNASequence with the quality values added as Feature.
     * The encoded quality values are automatically transformed into Phred
     * values.
     *
     * @return a biojava DNASequence with the quality values as feature
     */
    public DNASequence getDNASequence() {
        DNASequence seq = new DNASequence(sequence);
        seq.setOriginalHeader(description);
        QualityFeature feat = new QualityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound>("quality", "sequencing");
        feat.setQualities(getPhredQualities());
        seq.addFeature(1, seq.getLength(), feat);
        return seq;
    }

    private List<Number> getPhredQualities() {
        List<Number> qualities = new ArrayList<Number>();
        if (this.variant.equals(FastqVariant.FASTQ_SANGER)) {
            // Phred+33,  raw reads typically (0, 40)
            for (int i = 0; i < quality.length(); i++) {
                int charval = quality.charAt(i);
                qualities.add(charval - 33);
            }
        } else if (this.variant.equals(FastqVariant.FASTQ_SOLEXA)) {
            // Solexa+64, raw reads typically (-5, 40)
            for (int i = 0; i < quality.length(); i++) {
                int charval = quality.charAt(i);
                int phredvalue = getPhredFromSolexa(charval - 64);
                qualities.add(phredvalue);
            }
        } else if (this.variant.equals(FastqVariant.FASTQ_ILLUMINA)) {
            // Phred+64,  raw reads typically (0, 40) or (3,40)
            for (int i = 0; i < quality.length(); i++) {
                int charval = quality.charAt(i);
                qualities.add(charval - 64);
            }

        } else if (this.variant.equals(FastqVariant.FASTQ_NEW_ILLUMINA)) {
            // Phred+33,  raw reads typically (0, 41)
            for (int i = 0; i < quality.length(); i++) {
                int charval = quality.charAt(i);
                qualities.add(charval - 33);
            }
        } else {
            throw new IllegalArgumentException("unsupported variant");
        }

        return qualities;
    }

    private int getPhredFromError(double error) {
        return (int) Math.round(-10 * Math.log10(error));
    }

    private int getSolexaFromError(double error) {
        return (int) Math.round(-10 * Math.log10(error / (1 - error)));
    }

    private int getPhredFromSolexa(int solexa) {
        return (int) Math.round(10 * Math.log10(Math.pow(10, (solexa / 10.0)) + 1));
    }

    private int getSolexaFromPhred(int phred) {
        return (int) Math.round(10 * Math.log10(Math.pow(10, (phred / 10.0)) - 1));
    }

    
    private String encodeFeature(FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound> feature, FastqVariant variant) {
        if (feature instanceof QualityFeature) {
            if (feature.getLocations().getLength() != this.sequence.length()) {
                //TODO handle quality features with different locations more gracefully
                throw new IllegalArgumentException("DNASequence quality feature must cover the whole sequence");
            }
            QualityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> feat = (QualityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound>) feature;
            StringBuilder sb = new StringBuilder();
            List<Number> qualities = feat.getQualities();
            for (Number number : qualities) {
                sb.append(encodeQuality(number.intValue()));
            }
            return sb.toString();
        }
        return null;
    }

    private char encodeQuality(int quality) {
        if (this.variant.equals(FastqVariant.FASTQ_SANGER)) {
            char encoding = (char) (quality + 33);
            return encoding;
        } else if (this.variant.equals(FastqVariant.FASTQ_SOLEXA)) {
            int solexa = getSolexaFromPhred(quality);
            if (solexa < -5) {
                solexa = -5;
            }
            char encoding = (char) (solexa + 64);
            return encoding;
        } else if (this.variant.equals(FastqVariant.FASTQ_ILLUMINA)) {
            char encoding = (char) (quality + 64);
            return encoding;
        } else if (this.variant.equals(FastqVariant.FASTQ_NEW_ILLUMINA)) {
            char encoding = (char) (quality + 33);
            return encoding;
        } else {
            throw new IllegalArgumentException("unsupported variant");
        }
    }
}