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

package org.biojava.bio.seq.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.IllegalSymbolException;

/**
 * Defines what a sequence format does.
 *
 * <p>Sequence formats are responsible for both reading and writing a
 * sequence in a format, presumably in such a way as the written
 * record can be read back in by the same formatter.</p>
 *
 * <p>Where possible, the methods are parameterised so that they don't
 * need any knowledge of the specific implementation of Sequence they
 * are reading or writing. E.g. it should be possible to parameterise
 * readSequence to read from a Genbank stream and construct Ensembl
 * CORBA objects, just by specifying an Ensembl SequenceFactory.</p>
 * 
 * <p>More functionality is offered by {@link org.biojavax.bio.seq.io.RichSequenceFormat RichSequenceFormat},
 * Use of this interface is prefered.</p>
 *
 * @author Matthew Pocock
 * @author Thomas Down
 * @author Keith James
 */

public interface SequenceFormat
{
    /**
     * Read a sequence and pass data on to a SeqIOListener.
     *
     * @param reader The stream of data to parse.
     * @param symParser A SymbolParser defining a mapping from
     * character data to Symbols.
     * @param listener A listener to notify when data is extracted
     * from the stream.
     *
     * @return a boolean indicating whether or not the stream contains
     * any more sequences.
     *
     * @throws IOException if an error occurs while reading from the
     * stream.
     * @throws IllegalSymbolException if it is not possible to
     * translate character data from the stream into valid BioJava
     * symbols.
     * @throws BioException if there is an error in the format of the
     * stream.
     */
    public boolean readSequence(BufferedReader     reader,
                                SymbolTokenization symParser,
                                SeqIOListener      listener)
        throws BioException, IllegalSymbolException, IOException;

    /**
     * <code>writeSequence</code> writes a sequence to the specified
     * PrintStream, using the default format.
     *
     * @param seq the sequence to write out.
     * @param os the printstream to write to.
     */
    public void writeSequence(Sequence seq, PrintStream os)
        throws IOException;

    /**
     * <code>writeSequence</code> writes a sequence to the specified
     * <code>PrintStream</code>, using the specified format.
     *
     * @param seq a <code>Sequence</code> to write out.
     * @param format a <code>String</code> indicating which sub-format
     * of those available from a particular
     * <code>SequenceFormat</code> implemention to use when
     * writing.
     * @param os a <code>PrintStream</code> object.
     *
     * @exception IOException if an error occurs.
     * @deprecated use writeSequence(Sequence seq, PrintStream os)
     */
    public void writeSequence(Sequence seq, String format, PrintStream os)
        throws IOException;

    /**
     * <code>getDefaultFormat</code> returns the String identifier for
     * the default sub-format written by a <code>SequenceFormat</code>
     * implementation.
     *
     * @return a <code>String</code>.
     * @deprecated new implementations should only write a single
     * format.
     */
    public String getDefaultFormat();
}
