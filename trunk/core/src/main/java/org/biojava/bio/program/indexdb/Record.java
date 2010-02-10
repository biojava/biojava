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

package org.biojava.bio.program.indexdb;

import org.biojava.utils.io.RAF;

/**
 * <code>Record</code> represents a record within an indexed flat file
 * databank as defined by the OBDA standard.
 *
 * @author Matthew Pocock
 * @author Keith James
 */
public interface Record {

    /**
     * <code>getID</code> returns the primary identifier of the
     * record.
     *
     * @return a <code>String</code> ID.
     */
    public String getID();

    /**
     * <code>getFile</code> returns the random access file in which
     * the record belongs.
     *
     * @return a <code>RAF</code>.
     */
    public RAF getFile();

    /**
     * <code>getOffset</code> returns the byte offset in the file at
     * which the record begins.
     *
     * @return a <code>long</code> offset.
     */
    public long getOffset();

    /**
     * <code>getLength</code> returns the length of the record in
     * bytes.
     *
     * @return an <code>int</code>.
     */
    public int getLength();
  
    /**
     * <code>Impl</code> is the default implementation of Record.
     *
     * @author Matthew Pocock
     */
    public static class Impl
        implements Record {
        private final String id;
        private final RAF file;
        private final long offset;
        private final int length;
    
        /**
         * Creates a new <code>Impl</code> record.
         *
         * @param id a <code>String</code> primary ID.
         * @param file a <code>RAF</code> file.
         * @param offset a <code>long</code> byte offset.
         * @param length an <code>int</code> byte record length.
         */
        public Impl(String id, RAF file, long offset, int length) {
            if (id == null) {
                throw new NullPointerException("Can't have null ID");
            }
            this.id = id;
            this.file = file;
            this.offset = offset;
            this.length = length;
        }

        public String  getID()     { return id; }
        public RAF     getFile()   { return file; }
        public long    getOffset() { return offset; }
        public int     getLength() { return length; }
    }
}
