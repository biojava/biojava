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

package org.biojava.bio.program.abi;

import java.io.DataInput;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.io.CachingInputStream;
import org.biojava.utils.io.Seekable;

/**
 * A general base parser for files produced by ABI software.  This includes
 * chromatograms derived from ABI sequencers and potentially other data files
 * as well. The format was described by Clark Tibbetts in his paper "Raw Data
 * File Formats, and the Digital and Analog Raw Data Streams of the ABI PRISM
 * 377 DNA Sequencer."  Available online
 * <kbd><a href="http://www-2.cs.cmu.edu/afs/cs/project/genome/WWW/Papers/clark.html">
 * http://www-2.cs.cmu.edu/afs/cs/project/genome/WWW/Papers/clark.html</a></kbd>
 * <p>
 * Briefly, the format consists of a set of named fixed-length "tagged data
 * records" which may contain data themselves, or pointers to data elsewhere
 * in the file.  This class reads these records and exposes them to subclasses
 * through the {@link #getDataRecord} method.  The attributes of the records as
 * described in Tibbets' paper are exposed through public (final) fields of
 * {@link TaggedDataRecord} instances.
 * </p>
 * <p>
 * If a record only contains a pointer to the desired data (see 
 * {@link TaggedDataRecord#hasOffsetData}, subclasses may get
 * at the raw data by using {@link TaggedDataRecord#offsetData}:
 * </p>
 * <p>
 * This parser provides methods and classes for dealing with the files as
 * streams or local files (local files being more memory-efficient).
 * </p>
 *
 * @author Rhett Sutphin (<a href="http://genome.uiowa.edu/">UI CBCB</a>)
 * @author Richard Holland
 */
public class ABIFParser {
    private ABIFParser.DataAccess din;
    private boolean parsed = false;
    private Map records;

    private final int RECORD_COUNT_OFFSET  = 18;
    private final int RECORD_OFFSET_OFFSET = 26;

    /** Creates a new ABIFParser for a file. */
    public ABIFParser(File f) throws IOException {
        this(new ABIFParser.RandomAccessFile(f));
    }

    /**
     * Creates a new ABIFParser for an input stream.  Note that the stream
     * will be wrapped in a {@link CachingInputStream} if it isn't one already.
     * If it is, it will be seeked to 0.
     */
    public ABIFParser(InputStream in) throws IOException {
        this(new ABIFParser.DataStream(in));
    }

    /**
     * Creates a new ABIFParser for the specified {@link DataAccess} object.
     * If you need to read from something other than a file or a stream, you'll
     * have to implement a {@link DataAccess}-implementing class wrapping your
     * source and then pass an instance to this constructor.
     */
    public ABIFParser(ABIFParser.DataAccess toParse) throws IOException {
        din = toParse;
        readDataRecords();
    }

    /**
     * Returns the accessor for the raw data being parsed by this parser.
     */
    public final ABIFParser.DataAccess getDataAccess() {
        return din;
    }

    private final void readDataRecords() throws IOException {
        parsed = false;
        din.seek(RECORD_COUNT_OFFSET);
        long recordCount  = 0xffffffff & din.readInt();
        din.seek(RECORD_OFFSET_OFFSET);
        long recordOffset = 0xffffffff & din.readInt();
        din.seek(recordOffset);
        TaggedDataRecord tdr;
        StringBuffer label;
        records = new HashMap();
        for (int i = 0 ; i < recordCount ; i++) {
            tdr = new TaggedDataRecord(din);
            label = new StringBuffer(6).append(tdr.tagName).append(tdr.tagNumber);
            records.put(label.substring(0), tdr);
        }
        for (Iterator i = records.values().iterator(); i.hasNext(); ) {
        	TaggedDataRecord record = (TaggedDataRecord)i.next();
        	if (record.hasOffsetData) {
        		din.seek(record.dataRecord);
        		din.readFully(record.offsetData);
        	}
        }
        parsed = true;
        din.finishedReading();
    }

    /**
     * Decodes a character into a {@link Symbol} in the DNA alphabet.
     * Uses a definition of characters that is compatible with the ABI format.
     * @param token the character to decode
     * @throws IllegalSymbolException when token isn't in
     *         <code>{ a, A, c, C, g, G, t, T, n, N, - }</code>
     */
    public static Symbol decodeDNAToken(char token) throws IllegalSymbolException {
        switch (token) {
            case 'a': case 'A':
                return DNATools.a();
            case 'c': case 'C':
                return DNATools.c();
            case 'g': case 'G':
                return DNATools.g();
            case 't': case 'T':
                return DNATools.t();
            case 'n': case 'N':
                return DNATools.n();
            case '-':
                return DNATools.getDNA().getGapSymbol();
            default:
                throw new IllegalSymbolException("Can't decode token " + token + " into DNA");
        }
    }

    /**
     * Get the entry from the file TOC with the given name and tag number.
     * @param tagName the four-character string name of the desired data record
     * @param tagNumber which one of the tags with this name to return (must be positive)
     * @throws IllegalArgumentException if tagName is the wrong length or tagNumber
     *         is 0 or negative
     * @throws IllegalStateException if the initial parsing is not complete
     * @return the requested data record, or null if no such record exists
     */
    public ABIFParser.TaggedDataRecord getDataRecord(String tagName, int tagNumber)
    throws IllegalArgumentException, IllegalStateException {
        if (!parsed)
            throw new IllegalStateException("parsing is not complete");
        if (tagNumber < 1)
            throw new IllegalArgumentException("tagNumber must be positive");
        if (tagName.length() != 4)
            throw new IllegalArgumentException("tagName must be 4 characters long");
        return (ABIFParser.TaggedDataRecord) records.get(tagName + tagNumber);
    }
    
    /**
     * Obtain all data records. Keys of the map are strings consisting of
     * tag names with tag numbers concatenated immediately afterwards. Values
     * are TaggedDataRecord objects. The map has no particular order and so 
     * cannot be relied on to iterate over records in the same order they
     * were read from the file.
     * @return the map of all data records.
     */
    public Map getAllDataRecords() {
    	return Collections.unmodifiableMap(records);
    }

    /**
     * An aggregate immutable type for an ABIF tagged data record.  See the
     * Tibbets paper (referenced in the javadoc for {@link ABIFParser}) for
     * more information.
     */
    public static class TaggedDataRecord {
        public static final int DATA_TYPE_ASCII_ARRAY = 2;
        public static final int DATA_TYPE_INTEGER = 4;
        public static final int DATA_TYPE_FLOAT   = 7;
        public static final int DATA_TYPE_DATE    = 10;
        public static final int DATA_TYPE_TIME    = 11;
        public static final int DATA_TYPE_PSTRING = 18;

        public final char[] tagName;
        public final long   tagNumber;
        public final int    dataType;
        public final int    elementLength;
        public final long   numberOfElements;
        public final long   recordLength;
        public final long   dataRecord;
        public final long   crypticVariable;
        public final boolean hasOffsetData;
        public final byte[] offsetData;
        
        /**
         * Creates a new TaggedDataRecord from the next 28 bytes of
         * <code>din</code>.
         * @param din the source of the raw data to be parsed
         * @throws IOException if there's a problem with <code>din</code>
         */
        public TaggedDataRecord(ABIFParser.DataAccess din) throws IOException {
            tagName = new char[4];
            tagName[0] = (char) din.readByte();
            tagName[1] = (char) din.readByte();
            tagName[2] = (char) din.readByte();
            tagName[3] = (char) din.readByte();

            tagNumber        = 0xffffffff & din.readInt();
            dataType         = 0xffff & din.readShort();
            elementLength    = 0xffff & din.readShort();
            numberOfElements = 0xffffffff & din.readInt();
            recordLength     = 0xffffffff & din.readInt();
            dataRecord       = 0xffffffff & din.readInt();
            crypticVariable  = 0xffffffff & din.readInt();
            
            hasOffsetData = recordLength>4L;
            if (hasOffsetData)
            	offsetData = new byte[(int)recordLength];
            else 
            	offsetData = new byte[0];
        }

        /**
         * A very verbose <code>toString</code> that dumps all of the
         * data in this record in a human-readable format.
         */
        public String toString() {
            StringBuffer sb = new StringBuffer(super.toString()).append("[\n");
            sb.append("  tagName         = ").append(tagName).append('\n');
            sb.append("  tagNumber       = ").append(tagNumber).append('\n');
            sb.append("  dataType        = ");
            switch (dataType) {
                case DATA_TYPE_ASCII_ARRAY: sb.append("ASCII"); break;
                case DATA_TYPE_INTEGER: sb.append("INTEGER"); break;
                case DATA_TYPE_FLOAT:   sb.append("FLOAT");   break;
                case DATA_TYPE_DATE:    sb.append("DATE");    break;
                case DATA_TYPE_TIME:    sb.append("TIME");    break;
                case DATA_TYPE_PSTRING: sb.append("PSTRING"); break;
                default: sb.append(dataType);
            }
            sb.append('\n');
            sb.append("  elementLength   = ").append(elementLength).append('\n');
            sb.append("  numberOfElements= ").append(numberOfElements).append('\n');
            sb.append("  recordLength    = ").append(recordLength).append('\n');
            sb.append("  dataRecord      = ");
            if (recordLength <= 4) {
                switch (dataType) {
                case DATA_TYPE_ASCII_ARRAY:
                    if (recordLength > 3)
                        sb.append((char) ((dataRecord >>> 24) & 0xFF));
                    if (recordLength > 2)
                        sb.append((char) ((dataRecord >>> 16) & 0xFF));
                    if (recordLength > 1)
                        sb.append((char) ((dataRecord >>> 8 ) & 0xFF));
                    sb.append((char) ((dataRecord) & 0xFF));
                    break;
                case DATA_TYPE_DATE:
                    sb.append((dataRecord >>> 16) & 0xffff).append('/');
                    sb.append((dataRecord >>> 8 ) & 0xff).append('/');
                    sb.append((dataRecord) & 0xff);
                    break;
                case DATA_TYPE_TIME:
                    sb.append((dataRecord >>> 24) & 0xff).append(':');
                    sb.append((dataRecord >>> 16) & 0xff).append(':');
                    sb.append((dataRecord >>> 8 ) & 0xff);
                    break;
                case DATA_TYPE_INTEGER:
                    sb.append(dataRecord >>> (4 - recordLength)*8);
                    break;
                default:
                    hexStringify((int)dataRecord, sb);
                }
            }
            else {
                hexStringify((int)dataRecord, sb);
            }
            sb.append("  hasOffsetData   = ").append(hasOffsetData).append('\n');
            sb.append('\n');
            sb.append("  crypticVariable = ").append(crypticVariable).append('\n');
            sb.append(']');
            return sb.toString();
        }

        private void hexStringify(int l, StringBuffer sb) {
            sb.append("0x");
            String hex = Integer.toHexString(l).toUpperCase();
            for (int i = 8 ; i > hex.length() ; i--)
                sb.append('0');
            sb.append(hex);
        }
    }

    /**
     * Concatenation of the {@link Seekable} and {@link DataInput} interfaces.
     */
    public static interface DataAccess extends Seekable, DataInput { 
    	/**
    	 * Called when the parser has finished reading. The access
    	 * may choose to close itself at this point, e.g. if it is
    	 * using a RandomAccessFile.
    	 * @throws IOException if it could not do what it needs to.
    	 */
    	public void finishedReading() throws IOException;
    }

    private static class RandomAccessFile
    extends java.io.RandomAccessFile implements DataAccess {
        public RandomAccessFile(File f) throws FileNotFoundException {
            super(f, "r");
        }
        public void finishedReading() throws IOException {
        	this.close();
        }
    }

    /** Implements DataAccess by delegating to a CachingStream and a
     *  DataInputStream */
    private static class DataStream implements DataAccess {
        CachingInputStream cin;
        DataInputStream din;

        public DataStream(InputStream src) throws IOException {
            if (src instanceof CachingInputStream)
                cin = (CachingInputStream) src;
            else
                cin = new CachingInputStream(src);
            cin.seek(0);
            din = new DataInputStream(cin);
        }

        public DataStream(CachingInputStream cin) throws IOException {
            this((InputStream) cin);
        }
        
        public void finishedReading() throws IOException {
        	// We don't care.
        }

        public boolean readBoolean() throws IOException { return din.readBoolean(); }
        public byte    readByte()    throws IOException { return din.readByte();    }
        public char    readChar()    throws IOException { return din.readChar();    }
        public short   readShort()   throws IOException { return din.readShort();   }
        public int     readInt()     throws IOException { return din.readInt();     }
        public long    readLong()    throws IOException { return din.readLong();    }
        public float   readFloat()   throws IOException { return din.readFloat();   }
        public double  readDouble()  throws IOException { return din.readDouble();  }
        public String  readUTF()     throws IOException { return din.readUTF();     }

        public int readUnsignedByte()  throws IOException { return din.readUnsignedByte();  }
        public int readUnsignedShort() throws IOException { return din.readUnsignedShort(); }

        public void readFully(byte[] values) throws IOException {
            din.readFully(values);
        }

        public void readFully(byte[] values, int start, int len) throws IOException {
            din.readFully(values, start, len);
        }

        public String readLine() throws IOException {
            throw new UnsupportedOperationException("DataInputStream#readLine is deprecated.  Use readUTF instead");
        }

        public int skipBytes(int count) throws IOException { return din.skipBytes(count); }

        public void seek(long pos) throws IOException {
            cin.seek(pos);
        }
    }
}
