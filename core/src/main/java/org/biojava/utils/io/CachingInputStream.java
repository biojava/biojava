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

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;

/** 
 * A wrapper around {@link java.io.InputStream} that provides in-memory
 * caching of the input data.  This allows it to provide a {@link #seek}
 * method, which lets the user use an {@link java.io.InputStream} like a 
 * {@link java.io.RandomAccessFile} (with appropriate caveats about memory 
 * footprint, security, and performance).
 * <p>
 * This class has not been tested with very long input streams.  It might choke.
 * </p>
 *
 * @author Rhett Sutphin (<a href="http://genome.uiowa.edu/">UI CBCB</a>)
 */
public class CachingInputStream extends InputStream implements Seekable {
    private final static int INIT_CACHE_SIZE = 1024;
    private final static int RESIZE_FACTOR = 3;
    
    /** The byte cache itself. */
    protected byte[] cache;
    /** The 0-based index into cache of the _next_ byte to return.  If
     *  ptr == validLen, data must be read from the stream into the cache. */
    protected int ptr;
    /** A count of the number of bytes in {@link #cache} that contain
     *  data read from the stream. */
    protected int validLen;
    /** The underlying input stream whose data we're caching */
    protected InputStream in;
    
    public CachingInputStream(InputStream in) {
        this.in = in;
        cache = new byte[INIT_CACHE_SIZE];
        ptr = validLen = 0;
    }
    
    public void seek(long pos) throws IOException {
        if (pos > Integer.MAX_VALUE || pos < 0) {
            throw new IllegalArgumentException("Cannot seek to " 
                    + pos + ": can only do 0 <= seek < " + Integer.MAX_VALUE);
        }
        int newPtr = (int) pos;
        if (newPtr <= validLen) {
            ptr = newPtr;
        }
        else {
            skip(newPtr - ptr);
        }
    }
    
    public int read() throws IOException {
        if (ptr < validLen) {
            int out = 0xFF & cache[ptr];
            ptr++;
            return out;
        }
        else {
            int read = in.read();
            if (read >= 0) {
                expandCache(1);
                cache[ptr] = (byte) read;
                ptr++;
            }
            return read;
        }
    }
    
    public int read(byte[] b, int start, int len) throws IOException {
        int cachedLen = Math.min( Math.max(validLen - ptr, 0) , len );
        // copy the cached bytes to b, if any
        System.arraycopy(cache, ptr, b, start, cachedLen);
        ptr += cachedLen;
        // read additional bytes from the stream, if any
        int bytesRead = in.read(b, start + cachedLen, len - cachedLen);
        // copy newly read bytes into cache
        expandCache(bytesRead);
        System.arraycopy(b, start+cachedLen, cache, ptr, bytesRead);
        ptr += bytesRead;
        return bytesRead + cachedLen;
    }
    
    // FIXME: assumes ptr == validLen - 1
    public long skip(long num) throws IOException {
        if (ptr + num > Integer.MAX_VALUE)
            return 0;
        int n = (int) num;
        // skip through as much cache as exists up to n
        int availCache = Math.min(validLen - ptr, n);
        n -= availCache;
        ptr += availCache;
        // read any additional "skipped" bytes into cache
        expandCache(n);
        int i = 0, count;
        IOException ioEx = null;
        try {
            while (i < n) {
                count = in.read(cache, ptr + i, n - i);
                if (count < 0)
                    break;
                else
                    i += count;
            }
        }
        catch (EOFException e) { }
        catch (IOException e) { ioEx = e; }
        // if we couldn't skip enough bytes from the stream, 
        // mark those bytes invalid in the cache
        validLen -= (n - i);
        // update the pointer to indicate the skipped bytes
        ptr += i;
        // we save and rethrow the IOException in case the user of the code
        // tries to recover from the IOException -- this way there's a 
        // nonzero chance that ptr and validLen will still be valid
        if (ioEx != null)
            throw ioEx;
        // return the total number of bytes skipped with both methods
        return i + availCache;
    }
    
    /** Expands the cache to hold some number of <code>additionalBytes</code>.
     *  Expansion is done multiplicatively for efficiency. Immediately after
     *  calling this method, you must fill the additional bytes from the stream
     *  because this method also updates validLen.  */
    protected void expandCache(int additionalBytes) {
        if (cache.length < validLen + additionalBytes) {
            int newLen = cache.length;
            while (newLen < validLen + additionalBytes)
                newLen *= RESIZE_FACTOR;
            byte[] newCache = new byte[newLen];
            System.arraycopy(cache, 0, newCache, 0, cache.length);
            cache = newCache;
        }
        validLen += additionalBytes;
    }
}
