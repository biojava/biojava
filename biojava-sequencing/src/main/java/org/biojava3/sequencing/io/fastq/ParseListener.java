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
package org.biojava3.sequencing.io.fastq;

import java.io.IOException;

/**
 * Low-level event based parser callback.
 *
 * @since 3.0.3
 */
public interface ParseListener
{
    /**
     * Notify this parse listener of a description line.
     *
     * @param description description line
     * @throws IOException if an I/O error occurs
     */
    void description(String description) throws IOException;

    /**
     * Notify this parse listener of a sequence line.
     *
     * <p>
     * Note that the sequence in FASTQ format may contain end-of-line characters,
     * so both this method and <code>appendSequence(String)</code> may be called per FASTQ
     * formatted sequence.
     * </p>
     *
     * @param sequence sequence line
     * @throws IOException if an I/O error occurs
     */
    void sequence(String sequence) throws IOException;

    /**
     * Notify this parse listener of an additional sequence line.
     *
     * <p>
     * Note that the sequence in FASTQ format may contain end-of-line characters,
     * so this method may be called more than once per FASTQ formatted sequence.
     * </p>
     *
     * @param sequence additional sequence line
     * @throws IOException if an I/O error occurs
     */
    void appendSequence(String sequence) throws IOException;

    /**
     * Notify this parse listener of a repeat description line.
     *
     * @param repeatDescription repeat description line
     * @throws IOException if an I/O error occurs
     */
    void repeatDescription(String repeatDescription) throws IOException;

    /**
     * Notify this listener of a quality line.
     *
     * <p>
     * Note that the quality scores in FASTQ format may contain end-of-line characters,
     * so both this method and <code>appendQuality(String)</code> may be called per FASTQ
     * formatted sequence.
     * </p>
     *
     * @param quality quality line
     * @throws IOException if an I/O error occurs
     */
    void quality(String quality) throws IOException;

    /**
     * Notify this listener of a quality line.
     *
     * <p>
     * Note that the quality scores in FASTQ format may contain end-of-line characters,
     * so this method may be called more than once per FASTQ formatted sequence.
     * </p>
     *
     * @param quality additional quality line
     * @throws IOException if an I/O error occurs
     */
    void appendQuality(String quality) throws IOException;

    /**
     * Notify this listener the FASTQ formatted sequence is complete.
     *
     * @throws IOException if an I/O error occurs
     */
    void complete() throws IOException;
}
