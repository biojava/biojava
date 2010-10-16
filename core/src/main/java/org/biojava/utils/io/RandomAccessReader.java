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

package org.biojava.utils.io;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.io.Reader;

/**
 * <code>RandomAccessReader</code> extends <code>Reader</code> to
 * provide a means to create buffered <code>Reader</code>s from
 * <code>RandomAccessFile</code>s.
 *
 * @author Keith James
 * @since 1.2
 */
public class RandomAccessReader extends Reader
{
    private static final int DEFAULT_BUFFER_SIZE = 1 << 13;

    private RandomAccessFile raf;

    private char [] buffer;
    private byte [] bytes;

    private int bufferPos = 0;
    private int bufferEnd = 0;
    private long raPtrPos = 0;

    private boolean atEOF = false;

    /**
     * Creates a new <code>RandomAccessReader</code> wrapping the
     * <code>RandomAccessFile</code> and using a default-sized buffer
     * (8192 bytes).
     *
     * @param raf a <code>RandomAccessFile</code> to wrap.
     *
     * @exception IOException if an error occurs.
     */
    public RandomAccessReader(RandomAccessFile raf)
        throws IOException
    {
        this(raf, DEFAULT_BUFFER_SIZE);
    }

    /**
     * Creates a new <code>RandomAccessReader</code> wrapping the
     * <code>RandomAccessFile</code> and using a buffer of the
     * specified size.
     *
     * @param raf a <code>RandomAccessFile</code> to wrap.

     * @param sz an <code>int</code> buffer size.
     */
    public RandomAccessReader(RandomAccessFile raf, int sz)
        throws IOException
    {
        super();
        this.raf = raf;

        buffer = new char [sz];
        bytes  = new byte [sz];

        resetBuffer();
    }

    /**
     * <code>close</code> closes the underlying
     * <code>RandomAccessFile</code>.
     *
     * @exception IOException if an error occurs.
     */
    public void close() throws IOException
    {
        raf.close();
        raf = null;
    }

    /**
     * <code>length</code> returns the length of the underlying
     * <code>RandomAccessFile</code>.
     *
     * @return a <code>long</code>.
     *
     * @exception IOException if an error occurs.
     */
    public long length() throws IOException
    {
        return raf.length();
    }

    /**
     * <code>read</code> reads one byte from the underlying
     * <code>RandomAccessFile</code>.
     *
     * @return an <code>int</code>, -1 if the end of the stream has
     * been reached.
     *
     * @exception IOException if an error occurs.
     */
    public final int read() throws IOException
    {
	if (atEOF)
	    return -1;

        if (bufferPos >= bufferEnd)
            if (fill() < 0)
                return -1;

        if (bufferEnd == 0)
            return -1;
        else
            return buffer[bufferPos++];
    }

    /**
     * <code>read</code> reads from the underlying
     * <code>RandomAccessFile</code> into an array.
     *
     * @param cbuf a <code>char []</code> array to read into.
     * @param off an <code>int</code> offset in the array at which to
     * start storing chars.
     * @param len an <code>int</code> maximum number of char to read.
     *
     * @return an <code>int</code> number of chars read, or -1 if the
     * end of the stream has been reached.
     *
     * @exception IOException if an error occurs.
     */
    public int read(char [] cbuf, int off, int len) throws IOException
    {
	if (atEOF)
	    return -1;

        int remainder = bufferEnd - bufferPos;

        // If there are enough chars in the buffer to handle this
        // call, use those
        if (len <= remainder)
        {
            System.arraycopy(buffer, bufferPos, cbuf, off, len);
            bufferPos += len;

            return len;
        }

        // Otherwise start getting more chars from the delegate
        for (int i = 0; i < len; i++)
        {
            // Read from our own method which checks the buffer
            // first
            int c = read();

            if (c != -1)
            {
                cbuf[off + i] = (char) c;
            }
            else
	    {
		// Need to remember that EOF was reached to return -1
		// next read
		atEOF= true;

                return i;
	    }
        }

        return len;
    }

    /**
     * <code>getFilePointer</code> returns the effective position of
     * the pointer in the underlying <code>RandomAccessFile</code>.
     *
     * @return a <code>long</code> offset.
     *
     * @exception IOException if an error occurs.
     */
    public long getFilePointer() throws IOException
    {
        return raPtrPos - bufferEnd + bufferPos;
    }

    /**
     * <code>seek</code> moves the pointer to the specified position.
     *
     * @param pos a <code>long</code> offset.
     *
     * @exception IOException if an error occurs.
     */
    public void seek(long pos) throws IOException
    {
	// If we seek backwards after reaching EOF, we are no longer
	// at EOF.
	if (pos < raf.length())
	    atEOF = false;

        int p = (int) (raPtrPos - pos);

        // Check if we can seek within the buffer
        if (p >= 0 && p <= bufferEnd)
        {
            bufferPos = bufferEnd - p;
        }
        // Otherwise delegate to do a "real" seek and clean the
        // dirty buffer
        else
        {
            raf.seek(pos);
            resetBuffer();
        }
    }

    /**
     * <code>fill</code> fills the buffer from the
     * <code>RandomAccessFile</code>.
     *
     * @return an <code>int</code>.
     *
     * @exception IOException if an error occurs.
     */
    private int fill() throws IOException
    { 
        if (raf == null)
            throw new IOException("Random access file closed");

        // Read bytes from random access delegate
        int b = raf.read(bytes, 0, DEFAULT_BUFFER_SIZE);

        // Copy and cast bytes read to char buffer
        for (int i = b; --i >= 0;)
            buffer[i] = (char) bytes[i];

        // If read any bytes
        if (b >= 0)
        {
            raPtrPos += b;
            bufferPos = 0;
            bufferEnd = b;
        }

        // Return number bytes read
        return b;
    }

    /**
     * <code>resetBuffer</code> resets the buffer when the pointer
     * leaves its boundaries.
     *
     * @exception IOException if an error occurs.
     */
    private void resetBuffer() throws IOException
    {
        bufferPos = 0;
        bufferEnd = 0;
        raPtrPos = raf.getFilePointer();
    }
}
