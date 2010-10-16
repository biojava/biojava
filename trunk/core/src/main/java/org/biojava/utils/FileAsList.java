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

package org.biojava.utils;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.AbstractList;
import java.util.Comparator;
import java.util.Iterator;

/**
 * <code>FileAsList</code> creates a writable <code>List</code>
 * implementation backed by a random access file. There is a
 * restriction on the record length that the string representation of
 * that integer may not be longer than 4 bytes. This is because a
 * fixed 4 byte leader is used to encode the record length in the
 * file.
 *
 * @author Matthew Pocock
 * @author Keith James
 * @author Greg Cox
 */

// fixme: throughout this class, we are raising assertions for things that
// are legitimiate exceptions. This needs re-factoring.
public abstract class FileAsList
    extends
        AbstractList
    implements
        Commitable
{
    private static final int LEADER = 4;

    private RandomAccessFile mappedFile;
    private int commitedRecords;
    private int lastIndx = -1;
    private Object lastRec;
    private byte[] buffer;
    private int sizeCache = -1;

    /**
     * Creates a new <code>FileAsList</code> and corresponding backing
     * file.
     *
     * @param mappedFile a <code>File</code> used to back the
     * list. This file must not already exist.
     * @param recordLength an <code>int</code> byte record length.
     *
     * @exception IOException if an error occurs.
     */
    public FileAsList(File mappedFile, int recordLength)
        throws IOException {
        if(mappedFile.exists()) {
            throw new IOException("Can't create file as it already exists: " + mappedFile);
        }
        mappedFile.createNewFile();
        this.mappedFile = new RandomAccessFile(mappedFile, "rw");
        buffer = new byte[recordLength];
        this.mappedFile.seek(0L);
        byte[] rl = String.valueOf(recordLength).getBytes();
        if(rl.length > LEADER) {
            throw new IOException("Length of record too long"); // FIXME: ugg
        }
        for(int i = 0; i < rl.length; i++) {
            this.mappedFile.write(rl[i]);
        }
        for(int i = rl.length; i < LEADER; i++) {
            this.mappedFile.write(' ');
        }

        this.mappedFile.close();
    }

    /**
     * Creates a new <code>FileAsList</code> instance from an existing
     * backing file.
     *
     * @param mappedFile a <code>File</code> used to back the
     * list. This file must already exist.
     * @param mutable true if this list should support edits, false otherwise
     *
     * @exception IOException if an error occurs.
     */
    public FileAsList(File mappedFile, boolean mutable)
        throws IOException {
        if(!mappedFile.exists()) {
            throw new IOException("Can't load mapped list as the file does not exist: " + mappedFile);
        }

        if(mutable) {
          this.mappedFile = new RandomAccessFile(mappedFile, "rw");
        } else {
          this.mappedFile = new RandomAccessFile(mappedFile, "r");
        }
        StringBuffer sbuff = new StringBuffer();
        this.mappedFile.seek(0L);
        for(int i = 0; i < Math.min(LEADER, mappedFile.length()); i++) {
            char c = (char) this.mappedFile.readByte();
            sbuff.append(c);
        }

        buffer = new byte[Integer.parseInt(sbuff.substring(0).trim())];
    }

    /**
     * <code>rawGet</code> reads the record at the specified index as
     * a raw byte array.
     *
     * @param indx an <code>int</code> list index.
     *
     * @return a <code>byte []</code> array containing the raw record
     * data.
     */
    public byte[] rawGet(int indx) {
        if(indx < 0 || indx >= size()) {
            throw new IndexOutOfBoundsException("Can't access element: " + indx + " of " + size());
        }

        if(indx != lastIndx) {
            long offset = fixOffset(indx * buffer.length);
            try {
                mappedFile.seek(offset);
                mappedFile.readFully(buffer);
            } catch (IOException ioe) {
                throw new AssertionFailure("Failed to seek for record", ioe);
            }
        }

        return buffer;
    }

    public Object get(int indx) {
        if(indx == lastIndx) {
            return lastRec;
        }

        byte[] buffer = rawGet(indx);

        lastRec = parseRecord(buffer);
        lastIndx = indx;

        return lastRec;
    }

    public int size() {
        if(sizeCache < 0) {
            try {
                sizeCache = (int) (unFixOffset(mappedFile.length()) / (long) buffer.length);
            } catch (IOException ioe) {
                throw new AssertionFailure("Can't read file length", ioe);
            }
        };

        return sizeCache;
    }

    public boolean add(Object o) {
        sizeCache = -1;

        try {
            generateRecord(buffer, o);
        } catch (IOException e) {
            throw new AssertionFailure("Failed to write index", e);
        }

        try {
            mappedFile.seek(mappedFile.length());
            mappedFile.write(buffer);
        } catch (IOException ioe) {
            throw new AssertionFailure("Failed to write index", ioe);
        }

        return true;
    }

    /**
     * This always returns null, not the previous object.
     */
    public Object set(int indx, Object o) {
        try {
            generateRecord(buffer, o);
        } catch (IOException e) {
            throw new AssertionFailure("Failed to write index", e);
        }

        try {
            mappedFile.seek(fixOffset(indx * buffer.length));
            mappedFile.write(buffer);
        } catch (IOException ioe) {
            throw new AssertionFailure("Failed to write index", ioe);
        }

        return null;
    }

    public void clear() {
        try {
            mappedFile.setLength(fixOffset(0));
        } catch (IOException ioe) {
            throw new AssertionFailure("Could not truncate list", ioe);
        }
        commitedRecords = 0;
    }

    public void commit() {
        commitedRecords = this.size();
    }

    public void rollback() {
        try {
            mappedFile.setLength(fixOffset((long) commitedRecords * (long) buffer.length));
        } catch (Throwable t) {
            throw new AssertionFailure(
              "Could not roll back. "
              + "The index store will be in an inconsistent state "
              + "and should be discarded. File: "
              + mappedFile,
              t );
        }
    }

    private long fixOffset(long offset) {
        return offset + (long) LEADER;
    }

    private long unFixOffset(long offset) {
        return offset - (long) LEADER;
    }

    protected abstract Object parseRecord(byte[] buffer);

    protected abstract void generateRecord(byte[] buffer, Object item)
        throws IOException;

    public abstract Comparator getComparator();

    public Iterator iterator() {
        return new Iterator() {
                int i = 0;

                public Object next() {
                    return get(i++);
                }

                public boolean hasNext() {
                    return i < size();
                }

                public void remove() {}
            };
    }
}
