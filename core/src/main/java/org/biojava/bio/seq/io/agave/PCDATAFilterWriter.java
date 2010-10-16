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
package org.biojava.bio.seq.io.agave;
import java.io.FilterWriter;
import java.io.IOException;
import java.io.Writer;

/**
 * Write XML PCDATA using entity substitutions.
 *
 * @author Hanning Ni
 * @author Brian King
 */
class PCDATAFilterWriter extends FilterWriter
{
    /** entity substitution for <     */
    protected static char [] LT = {'&', 'l', 't', ';'};

    /** entity substitution for >    */
    protected static char [] GT = {'&', 'g', 't', ';'};

    /** entity substitution for >    */
    protected static char [] AMP = {'&', 'a', 'm', 'p', ';'};

    /**
     * Constructor.
     */
    public PCDATAFilterWriter(Writer out)
    {
            super(out);
    }

    /**
     * Write a single character.
     *
     * @exception  IOException  If an I/O error occurs
     */
    public void write(int c) throws IOException
    {
        synchronized(lock)
        {
            if(c == '<')
            {
                out.write(LT);
            }
            else if(c == '>')
            {
                out.write(GT);
            }
            else if(c == '&')
            {
                out.write(AMP);
            }
            else
            {
                out.write(c);
            }
        }
    }

    /**
     * Write a portion of an array of characters.
     *
     * @param  cbuf  Buffer of characters to be written
     * @param  off   Offset from which to start reading characters
     * @param  len   Number of characters to be written
     *
     * @exception  IOException  If an I/O error occurs
     */
    public void write(char cbuf[], int off, int len) throws IOException
    {
        synchronized(lock)
        {
            for (int i = 0; i < len; i++)
            {
                write(cbuf[off + i]);
            }
        }
    }

    /**
     * Write a portion of a string.
     *
     * @param  str  String to be written
     * @param  off  Offset from which to start reading characters
     * @param  len  Number of characters to be written
     *
     * @exception  IOException  If an I/O error occurs
     */
    public void write(String str, int off, int len) throws IOException
    {
        synchronized(lock)
        {
            for (int i = 0; i < len; i++)
            {
                    write(str.charAt(off+i));
            }
        }
    }
    /**
     * Write a  string.
     *
     * @param  str  String to be written
     *
     * @exception  IOException  If an I/O error occurs
     */
    public void write(String str) throws IOException
    {
        synchronized(lock)
        {
            int len = str.length();
            for (int i = 0; i < len; i++)
            {
                    int c = str.charAt(i);
                if(c == '<')
                {
                    out.write(LT);
                }
                else if(c == '>')
                {
                    out.write(GT);
                }
                else if(c == '&')
                {
                    out.write(AMP);
                }
                else
                {
                    out.write(c);
                }
            }
        }
    }
}
