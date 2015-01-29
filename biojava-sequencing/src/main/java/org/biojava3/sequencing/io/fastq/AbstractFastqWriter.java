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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Writer;

import java.util.Arrays;

/**
 * Abstract writer implementation for FASTQ formatted sequences.
 *
 * @since 3.0.3
 */
abstract class AbstractFastqWriter
    implements FastqWriter
{

    /**
     * Validate the specified FASTQ formatted sequence for writing.
     *
     * @param fastq FASTQ formatted sequence to validate, will not be null
     * @throws IOException if the specified FASTQ formatted sequence is not valid for writing
     */
    protected abstract void validate(final Fastq fastq) throws IOException;

    @Override
    public final <T extends Appendable> T append(final T appendable, final Fastq... fastq) throws IOException
    {
        return append(appendable, Arrays.asList(fastq));
    }

    @Override
    public final <T extends Appendable> T append(final T appendable, final Iterable<Fastq> fastq) throws IOException
    {
        if (appendable == null)
        {
            throw new IllegalArgumentException("appendable must not be null");
        }
        if (fastq == null)
        {
            throw new IllegalArgumentException("fastq must not be null");
        }
        for (Fastq f : fastq)
        {
            validate(f);
            if (f != null)
            {
                appendable.append("@");
                appendable.append(f.getDescription());
                appendable.append("\n");
                appendable.append(f.getSequence());
                appendable.append("\n");
                appendable.append("+\n");
                appendable.append(f.getQuality());
                appendable.append("\n");
            }
        }
        return appendable;
    }

    @Override
    public final void write(final File file, final Fastq... fastq) throws IOException
    {
        write(file, Arrays.asList(fastq));
    }

    @Override
    public final void write(final File file, final Iterable<Fastq> fastq) throws IOException
    {
        if (file == null)
        {
            throw new IllegalArgumentException("file must not be null");
        }
        if (fastq == null)
        {
            throw new IllegalArgumentException("fastq must not be null");
        }
        Writer writer = null;
        try
        {
            writer = new BufferedWriter(new FileWriter(file));
            append(writer, fastq);
        }
        finally
        {
            if (writer != null)
            {
                try
                {
                    writer.close();
                }
                catch (IOException e)
                {
                    // ignore
                }
            }
        }
    }

    @Override
    public final void write(final OutputStream outputStream, final Fastq... fastq) throws IOException
    {
        write(outputStream, Arrays.asList(fastq));
    }

    @Override
    public final void write(final OutputStream outputStream, final Iterable<Fastq> fastq) throws IOException
    {
        if (outputStream == null)
        {
            throw new IllegalArgumentException("outputStream must not be null");
        }
        if (fastq == null)
        {
            throw new IllegalArgumentException("fastq must not be null");
        }
        Writer writer = null;
        try
        {
            writer = new BufferedWriter(new OutputStreamWriter(outputStream));
            append(writer, fastq);
        }
        finally
        {
            if (writer != null)
            {
                try
                {
                    writer.flush();
                }
                catch (IOException e)
                {
                    // ignore
                }
            }
        }
    }
}