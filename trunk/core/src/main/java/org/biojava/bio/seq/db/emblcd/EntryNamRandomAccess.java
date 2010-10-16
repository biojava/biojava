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

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

/**
 * <code>EntryNamRandomAccess</code> objects provide random access to
 * records within the "entrynam.idx" file of an EMBL CD-ROM format
 * binary index. Records are retrieved by their sequence ID.
 *
 * @author Keith James
 * @since 1.2
 */
public class EntryNamRandomAccess extends EmblCDROMRandomAccess
{
    public EntryNamRandomAccess(File indexFile,
                                int  headerLength,
                                int  recordLength,
                                long recordCount)
        throws FileNotFoundException
    {
        super(indexFile, headerLength, recordLength, recordCount);
    }

    /**
     * <code>readRecord</code> creates an array of Objects from the
     * raw byte array of a single record. For this file type the array
     * contains <code>String seqID, Long rPosition, Long sPosition,
     * Integer fileNumber</code>. See EMBOSS documentation for a full
     * description.
     *
     * @return an <code>Object []</code> array.
     *
     * @exception IOException if an error occurs.
     */
    protected Object [] readRecord() throws IOException
    {
        int eof = raIndexFile.read(recBytes);
        if (eof == -1)
            raIndexFile.close();

        return recParser.parseEntryNamRecord(recBytes);
    }

    /**
     * <code>getRecordKey</code> returns the field from the record on
     * which the records were sorted in the index. (i.e. sequence ID
     * or accession number).
     *
     * @return a <code>String</code>.
     */
    protected String getRecordKey(Object [] record)
    {
        return (String) record[0];
    }
}
