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

import java.net.URL;

import java.io.Closeable;
import java.io.File;
import java.io.InputStream;
import java.io.IOException;

import com.google.common.io.InputSupplier;

/**
 * Reader for FASTQ formatted sequences.
 *
 * @since 3.0.3
 */
public interface FastqReader
{

    /**
     * Parse the specified input supplier.
     *
     * @param supplier input supplier, must not be null
     * @param listener low-level event based parser callback, must not be null
     * @throws IOException if an I/O error occurs
     */
    <R extends Readable & Closeable> void parse(InputSupplier<R> supplier, ParseListener listener) throws IOException;

    /**
     * Stream the specified input supplier.
     *
     * @param supplier input supplier, must not be null
     * @param listener event based reader callback, must not be null
     * @throws IOException if an I/O error occurs
     */
    <R extends Readable & Closeable> void stream(InputSupplier<R> supplier, StreamListener listener) throws IOException;

    /**
     * Read zero or more FASTQ formatted sequences from the specified file.
     *
     * @param file file to read from, must not be null
     * @return zero or more FASTQ formatted sequences read from the specified file
     * @throws IOException if an I/O error occurs
     */
    Iterable<Fastq> read(File file) throws IOException;

    /**
     * Read zero or more FASTQ formatted sequences from the specified url.
     *
     * @param url URL to read from, must not be null
     * @return zero or more FASTQ formatted sequences read from the specified url
     * @throws IOException if an I/O error occurs
     */
    Iterable<Fastq> read(URL url) throws IOException;

    /**
     * Read zero or more FASTQ formatted sequences from the specified input stream.
     *
     * @param inputStream input stream to read from, must not be null
     * @return zero or more FASTQ formatted sequences read from the specified input stream
     * @throws IOException if an I/O error occurs
     */
    Iterable<Fastq> read(InputStream inputStream) throws IOException;
}