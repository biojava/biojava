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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;

/**
 * @author Thomas Down
 */
    public class CountedBufferedReader extends BufferedReader {
	private final static int DEFAULT_BUFFER_SIZE = 1 << 14;

	private long position;

	private Reader stream;
	private char[] buffer;
	private int buffPos;
	private int buffFill;

	private boolean reachedEOF = false;

	private int mark = -1;
	private int markLimit = -1;

	public long getFilePointer() {
	    return position;
	}

	public CountedBufferedReader(Reader stream) {
	    super(new Reader() {
		public void close() {}
		public int read(char[] cbuf, int off, int len) { return 0; }
	    });

	    this.stream = stream;
	    this.buffer = new char[DEFAULT_BUFFER_SIZE];
	    position = 0;
	    buffPos = 0;
	    buffFill = 0;
	}

	public void close()
	    throws IOException
	{
	    stream.close();
	    stream = null;
	}

	public int read() 
	    throws IOException
	{
	    if (buffPos >= buffFill)
		fillBuffer();

	    if (reachedEOF) {
		return -1;
	    } else {
		position++;
		return buffer[buffPos++];
	    }
	}

	private int peek()
	    throws IOException
	{
	    if (buffPos >= buffFill)
		fillBuffer();

	    if (reachedEOF) {
		return -1;
	    } else {
		return buffer[buffPos];
	    }
	}

	public int read(char[] cbuf)
	    throws IOException
	{
	    return read(cbuf, 0, cbuf.length);
	}

	public int read(char[] cbuf, int off, int len) 
	    throws IOException
	{
	    if (buffPos >= buffFill)
		fillBuffer();

	    if (reachedEOF) {
		return -1;
	    } else {
		int toReturn = Math.min(len, buffFill - buffPos);
		System.arraycopy(buffer, buffPos, cbuf, off, toReturn);
		buffPos += toReturn;
		position += toReturn;

		return toReturn;
	    }
	}

	public boolean ready() 
	    throws IOException
	{
	    if (buffPos < buffFill)
		return true;
	    return stream.ready();
	}

	public long skip(long n)
	    throws IOException
	{
	    int skipInBuffer;
	    if (n < buffer.length) {
		skipInBuffer = Math.min((int) n, buffFill - buffPos);
	    } else {
		skipInBuffer = buffFill - buffPos;
	    }
	    position += skipInBuffer;
	    buffPos += skipInBuffer;

	    if (n > skipInBuffer) {
		long skippedInStream;

		if (mark >= 0) {
		    // Yuck, fix this...
		    char[] dummy = new char[(int) (n - skipInBuffer)];
		    skippedInStream = read(dummy); 
		} else {
		    skippedInStream = stream.skip(n - skipInBuffer);
		}

		position += skippedInStream;
		return skippedInStream + skipInBuffer;
	    } else {
		return skipInBuffer;
	    }
	}

	public boolean markSupported() {
	    return true;
	}

	public void mark(int limit)
	    throws IOException
	{
	    //	    System.err.println("*** Mark");

	    if (limit + 1> buffer.length) {
		char[] newBuffer = new char[limit + 1];
		System.arraycopy(buffer, buffPos, newBuffer, 0, buffFill - buffPos);
		buffer = newBuffer;
		buffFill = buffFill - buffPos;
		buffPos = 0;
	    } else if (buffPos + limit > buffer.length) {
		System.arraycopy(buffer, buffPos, buffer, 0, buffFill - buffPos);
		buffFill = buffFill - buffPos;
		buffPos = 0;
	    }

	    mark = buffPos;
	    markLimit = limit;
	}

	public void reset()
	    throws IOException
	{
	    //	    System.err.println("*** Reset");

	    if (mark < 0)
		throw new IOException("The mark is not currently in scope");

	    position = position - buffPos + mark;
	    buffPos = mark;
	}

	public String readLine()
	    throws IOException 
	{
	    String line = null;
	    
      while(!reachedEOF) {
        for(int i = buffPos; i < buffFill; i++) {
          char c = buffer[i];
          if(c == '\n' || c == '\r') {
            int len = i - buffPos;
            String bit = new String(buffer, buffPos, len);
            position += len;
            if(line == null) {
              line = bit;
            } else {
              line += bit;
            }
            buffPos = i;
            read();
            char d = (char) peek();
            if(c == '\r' && d == '\n') {
              read();
            }
            
            return line;
          }
        }
        int len = buffFill - buffPos;
        String bit = new String(buffer, buffPos, len);
        position += len;
        buffPos = buffFill;
        if(line == null) {
          line = bit;
        } else {
          line += bit;
        }
        fillBuffer();
      }
      
      return line;
	}

	private void fillBuffer()
	    throws IOException
	{
	    // System.err.println("*** Fill buffer");

	    if (mark < 0) {
		buffFill = stream.read(buffer);
		if (buffFill == -1) {
		    buffFill = 0;
		    reachedEOF = true;
		} 
		// System.out.println("Filled buffer: " + buffFill);
		
		buffPos = 0;
	    } else {
		if (buffPos >= (markLimit + mark)) {
		    // Mark's gone out of scope -- wheee!
		    mark = -1;
		    markLimit = -1;
		    fillBuffer();
		    return;
		}

		System.arraycopy(buffer, mark, buffer, 0, buffFill - mark);
		buffFill = buffFill - mark;
		mark = 0;
		buffPos = buffFill;
		int newChars = stream.read(buffer, buffFill, buffer.length - buffFill);
		if (newChars == -1) {
		    reachedEOF = true;
		} else {
		    buffFill = buffFill + newChars;
		}
	    }
	}
    }

