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

package org.biojava.bio.seq.db.emblcd;

import java.io.IOException;
import java.io.InputStream;

/**
 * <code>AcnumTrgReader</code> reads the "acnum.trg" file of an EMBL
 * CD-ROM format binary index.
 *
 * @author Keith James
 * @since 1.2
*/
public class AcnumTrgReader extends EmblCDROMIndexReader
{
    /**
     * Creates a new <code>AcnumTrgReader</code>.
     *
     * @param input an <code>InputStream</code>.
     *
     * @exception IOException if an error occurs.
     */
    public AcnumTrgReader(InputStream input)
        throws IOException
    {
        super(input);
    }

    /**
     * <code>readRecord</code> creates an array of Objects from the
     * raw byte array of a single record.
     *
     * @return an <code>Object []</code> array.
     *
     * @exception IOException if an error occurs.
     */
    public Object [] readRecord() throws IOException
    {
        return recParser.parseAcnumTrgRecord(readRawRecord());
    }
}
