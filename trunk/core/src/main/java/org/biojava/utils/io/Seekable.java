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

/** 
 * This interface provides a collective name for IO classes that implement a
 * seek function (e.g., {@link java.io.RandomAccessFile}).
 *
 * @author Rhett Sutphin (<a href="http://genome.uiowa.edu/">UI CBCB</a>)
 */
public interface Seekable {
    /** 
     * Moves the pointer in the inputstream such that the byte starting  at 
     * <code>pos</code> are returned by the next read.
     * @param pos the position to which to seek
     * @throws IOException when there's an I/O problem
     */
    public void seek(long pos) throws IOException;
}
