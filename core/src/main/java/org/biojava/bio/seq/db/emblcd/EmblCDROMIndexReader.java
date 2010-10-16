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
 * <p><code>EmblCDROMIndexReader</code> is an abstract class whose
 * concrete subclasses read EMBL CD-ROM format indices from an
 * underlying <code>InputStream</code>. This format is used by the
 * EMBOSS package for database indexing (see programs dbiblast,
 * dbifasta, dbiflat and dbigcg). Indexing produces four binary files
 * with a simple format:</p>
 * 
 * <ul>
 *   <li>division.lkp : master index</li>
 *   <li>entrynam.idx : sequence ID index</li>
 *   <li>   acnum.trg : accession number index</li>
 *   <li>   acnum.hit : accession number auxiliary index</li>
 * </ul>
 *
 * <p>Internally EMBOSS checks for Big-endian architechtures and
 * switches the byte order to Little-endian. This means trouble if you
 * try to read the file using <code>DataInputStream</code>, but at
 * least the binaries are consistent across architechtures. This class
 * carries out the necessary conversion.</p>
 *
 * <p>The EMBL CD-ROM format stores the date in 4 bytes. One byte is
 * unused (the first one), leaving one byte for the day, one for the
 * month and one (!) for the year.</p>
 *
 * <p> For further information see the EMBOSS documentation, or for a
 * full description, the source code of the dbi programs and the Ajax
 * library.</p>
 *
 * @author Keith James
 * @since 1.2
 */
public abstract class EmblCDROMIndexReader
{
    protected InputStream  input;
    protected StringBuffer sb;
    protected RecordParser recParser;

    // Header fields
    private byte []      int4 = new byte [4];
    private byte []      int2 = new byte [2];
    private byte []    dbName = new byte [20];
    private byte [] dbRelease = new byte [10];
    private byte []    dbDate = new byte [4];
    // Record field
    private byte []    record;

    private long   fileLength;
    private long   recordCount;
    private int    recordLength;
    private String name;
    private String release;
    private String date;

    /**
     * Creates a new <code>EmblCDROMIndexReader</code> instance. A
     * <code>BufferedInputStream</code> is probably the most suitable.
     *
     * @param input an <code>InputStream</code>.
     *
     * @exception IOException if an error occurs.
     */
    public EmblCDROMIndexReader(InputStream input)
        throws IOException
    {
        this.input = input;
        sb = new StringBuffer(512);
        recParser = new RecordParser();

        parseHeader();
    }

    /**
     * <code>readFileLength</code> returns the file length in bytes
     * (stored within the file's header by the indexing program). This
     * may be called more than once as the value is cached.
     *
     * @return a <code>long</code>.
     */
    public long readFileLength()
    {
        return fileLength;
    }

    /**
     * <code>readRecordCount</code> returns the number of records in
     * the file. This may be called more than once as the value is
     * cached.
     *
     * @return a <code>long</code>.
     */
    public long readRecordCount()
    {
        return recordCount;
    }

    /**
     * <code>readRecordLength</code> returns the record length
     * (bytes). This may be called more than once as the value is
     * cached.
     *
     * @return an <code>int</code>.
     */
    public int readRecordLength()
    {
        return recordLength;
    }

    /**
     * <code>readDBName</code> returns the database name from the
     * index header. This may be called more than once as the value is
     * cached.
     *
     * @return a <code>String</code>.
     */
    public String readDBName()
    {
        return name;
    }

    /**
     * <code>readDBRelease</code> returns the database release from
     * the index header. This may be called more than once as the
     * value is cached.
     *
     * @return a <code>String</code>.
     */
    public String readDBRelease()
    {
        return release;
    }

    /**
     * <code>readDBDate</code> reads the date from the index
     * header. The date is stored in 4 bytes: 0, unused; 1, year; 2,
     * month; 3, day. With a 1 byte year it's not very much use and
     * I'm not sure that the EMBOSS programs set the value correctly
     * anyway.
     *
     * @return a <code>String</code>.
     */
    public String readDBDate()
    
    {
        return date;
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
    public abstract Object [] readRecord() throws IOException;

    /**
     * <code>readRawRecord</code> returns the raw bytes of a single
     * record from the index.
     *
     * @return a <code>byte []</code> array.
     *
     * @exception IOException if an error occurs.
     */
    public byte [] readRawRecord() throws IOException
    {
        int eof = input.read(record);
        if (eof == -1)
            input.close();

        return record;
    }

    /**
     * <code>close</code> closes the underlying
     * <code>InputStream</code>.
     *
     * @exception IOException if an error occurs.
     */
    public void close() throws IOException
    {
        input.close();
    }

    /**
     * <code>parseHeader</code> carries out a full parse of the 300
     * byte header (common to all the index types) when first
     * encountered.
     *
     * @exception IOException if an error occurs.
     */
    private void parseHeader() throws IOException
    {
        int eof = 0;

        eof = input.read(int4);
        if (eof == -1)
            input.close();

        fileLength = recParser.parseInt4(int4);

        eof = input.read(int4);
        if (eof == -1)
            input.close();

        recordCount = recParser.parseInt4(int4);

        eof = input.read(int2);
        if (eof == -1)
            input.close();

        recordLength = recParser.parseInt2(int2);

        // Set up array for reading records now that we know their
        // length
        record = new byte [recordLength];

        eof = input.read(dbName);
        if (eof == -1)
            input.close();

        sb.setLength(0);
        name = recParser.parseString(sb, dbName);

        eof = input.read(dbRelease);
        if (eof == -1)
            input.close();

        sb.setLength(0);
        release = recParser.parseString(sb, dbRelease);

        eof = input.read(dbDate);
        if (eof == -1)
            input.close();

        sb.setLength(0);
        date = recParser.parseDate(sb, dbDate);

        // Skip the remainder of the header (padding)
        input.skip(256);
    }
}
