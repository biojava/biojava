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

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

/**
 * Writer for FASTQ formatted sequences.
 *
 * @since 3.0.3
 */
public interface FastqWriter
{

    /**
     * Append the specified FASTQ formatted sequences to the specified appendable.
     *
     * @param <T> extends Appendable
     * @param appendable appendable to append the specified FASTQ formatted sequences to, must not be null
     * @param fastq variable number of FASTQ formatted sequences to append, must not be null
     * @return the specified appendable with the specified FASTQ formatted sequences appended
     * @throws IOException if an I/O error occurs
     */
    <T extends Appendable> T append(T appendable, Fastq... fastq) throws IOException;

    /**
     * Append the specified FASTQ formatted sequences to the specified appendable.
     *
     * @param <T> extends Appendable
     * @param appendable appendable to append the specified FASTQ formatted sequences to, must not be null
     * @param fastq zero or more FASTQ formatted sequences to append, must not be null
     * @return the specified appendable with the specified FASTQ formatted sequences appended
     * @throws IOException if an I/O error occurs
     */
    <T extends Appendable> T append(T appendable, Iterable<Fastq> fastq) throws IOException;

    /**
     * Write the specified FASTQ formatted sequences to the specified file.
     *
     * @param file file to write to, must not be null
     * @param fastq variable number of FASTQ formatted sequences to write, must not be null
     * @throws IOException if an I/O error occurs
     */
    void write(File file, Fastq... fastq) throws IOException;

    /**
     * Write the specified FASTQ formatted sequences to the specified file.
     *
     * @param file file to write to, must not be null
     * @param fastq zero or more FASTQ formatted sequences to write, must not be null
     * @throws IOException if an I/O error occurs
     */
    void write(File file, Iterable<Fastq> fastq) throws IOException;

    /**
     * Write the specified FASTQ formatted sequences to the specified output stream.
     *
     * @param outputStream output stream to write to, must not be null
     * @param fastq variable number of FASTQ formatted sequences to write, must not be null
     * @throws IOException if an I/O error occurs
     */
    void write(OutputStream outputStream, Fastq... fastq) throws IOException;

    /**
     * Write the specified FASTQ formatted sequences to the specified output stream.
     *
     * @param outputStream output stream to write to, must not be null
     * @param fastq zero or more FASTQ formatted sequences to write, must not be null
     * @throws IOException if an I/O error occurs
     */
    void write(OutputStream outputStream, Iterable<Fastq> fastq) throws IOException;
}