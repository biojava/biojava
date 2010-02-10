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
import java.io.RandomAccessFile;

/**
 * <p><code>EmblCDROMRandomAccess</code> is an abstract class whose
 * concrete subclasses can perform fast lookups in EMBL CD-ROM format
 * index files. As the format of the records varies between file
 * types, subclasses should implement two methods;
 * <code>readRecord()</code>, which should parse the record into an
 * array of objects and <code>getRecordKey()</code> which should
 * retrieve the the field from the parsed record on which the records
 * were sorted in the index. This is used during the binary search in
 * the <code>findRecord()</code> method.</p>
 *
 * <p>Implementing <code>readRecord()</code> is easy because it simply
 * means delegating to the supplied <code>RecordParser</code> and
 * calling the appropriate method on it.</p>
 *
 * @author Keith James
 * @since 1.2
 */
public abstract class EmblCDROMRandomAccess
{
    private   File             indexFile;
    protected RandomAccessFile raIndexFile;

    private int headerLength;
    private int recordLength;
    private long recordCount;

    /**
     * A <code>recParser</code> for implementing
     * <code>readRecord()</code> specific to each concrete subclass.
     */
    protected RecordParser recParser;

    protected byte [] recBytes;

    /**
     * Creates a new <code>EmblCDROMRandomAccess</code> object.
     *
     * @param indexFile a <code>File</code> to wrap.
     * @param headerLength an <code>int</code> (normally 300 bytes).
     * @param recordLength an <code>int</code> indicating the length
     * of a single record.
     * @param recordCount an <code>long</code> indicating the total
     * number of records.
     *
     * @exception FileNotFoundException if indexFile cannot be found.
     */
    public EmblCDROMRandomAccess(File indexFile,
                                 int  headerLength,
                                 int  recordLength,
                                 long recordCount)
        throws FileNotFoundException
    {
        this.indexFile = indexFile;
        raIndexFile = new RandomAccessFile(indexFile, "r");

        this.headerLength = headerLength;
        this.recordLength = recordLength;
        this.recordCount  = recordCount;

        recBytes  = new byte [recordLength];
        recParser = new RecordParser();
    }

    /**
     * <code>getFile</code> returns the <code>File</code> wrapped.
     *
     * @return a <code>File</code>.
     */
    public File getFile()
    {
        return indexFile;
    }

    /**
     * <code>findRecord</code> performs a binary search within the
     * file for a record specified by an identifier String.
     *
     * @param identifier a <code>String</code> identifier (sequence ID
     * or accession number).
     *
     * @return an <code>Object []</code> array containing the
     * record. If there is no such record an empty array is returned.
     *
     * @exception IOException if an error occurs.
     */
    public Object [] findRecord(String identifier)
        throws IOException
    {
        long startRecord = 0;
        long  endRecord  = recordCount - 1;

        while (startRecord <= endRecord)
        {
            long midPoint = (startRecord + endRecord) / 2;
            raIndexFile.seek(headerLength + (midPoint * recordLength));

            Object [] record = readRecord();
            String recordKey = getRecordKey(record).trim();

            if (recordKey.equals(identifier))
                return record;
            else if (recordKey.compareTo(identifier) < 0)
                startRecord = midPoint + 1;
            else
                endRecord = midPoint - 1;
        }

        // No such record
        return new Object [0];
    }

    /**
     * <code>close</code> closes the underlying
     * <code>RandomAccessFile</code>.
     *
     * @exception IOException if an error occurs.
     */
    public void close() throws IOException
    {
        raIndexFile.close();
    }

    /**
     * <code>readRecord</code> returns an array of objects parsed from
     * a single record. Its content will depend on the type of index
     * file. Concrete subclasses must provide an implementation of
     * this method.
     *
     * @return an <code>Object []</code> array.
     *
     * @exception IOException if an error occurs.
     */
    protected abstract Object [] readRecord() throws IOException;

    /**
     * <code>getRecordKey</code> returns the field from the record on
     * which the records were sorted in the index. (i.e. sequence ID
     * or accession number).
     *
     * @return a <code>String</code>.
     */
    protected abstract String getRecordKey(Object [] record);
}
