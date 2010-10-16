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

package org.biojava.bio.chromatogram;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

import org.biojava.bio.program.abi.ABIFChromatogram;
import org.biojava.bio.program.scf.SCF;
import org.biojava.utils.io.CachingInputStream;

/**
 * A factory that creates {@link Chromatogram} objects from files or streams.
 * In either case, the type of object to create is determined from the first
 * four bytes of the stream (the magic number).
 *
 * @author Rhett Sutphin (<a href="http://genome.uiowa.edu/">UI CBCB</a>)
 * @author Matthew Pocock
 * @since 1.3
 */
public class ChromatogramFactory {
    /**
     *  The magic number for SCF files.
     */
    public static final int SCF_MAGIC = (((byte) '.') << 24)
                                      + (((byte) 's') << 16)
                                      + (((byte) 'c') << 8)
                                      + (((byte) 'f'));
    /**
     *  The magic number for ABIF files.
     */
    public static final int ABI_MAGIC = (((byte) 'A') << 24)
                                      + (((byte) 'B') << 16)
                                      + (((byte) 'I') << 8)
                                      + (((byte) 'F'));

    /**
     * Creates a new <code>Chromatogram</code> object from the named file.
     * @param f the file to read
     * @return a new Chromatogram
     *
     * @throws IOException when the file can't be read or some other I/O error occurs
     * @throws UnsupportedChromatogramFormatException when the file doesn't
     *         contain a chromatogram in a supported format
     */
    public static Chromatogram create(File f)
    throws IOException, UnsupportedChromatogramFormatException {
        FileInputStream fin = new FileInputStream(f);
        int magic = magicFromStream(fin);
        fin.close();

        switch (magic) {
        case SCF_MAGIC:
            return SCF.create(f);
        case ABI_MAGIC:
            return ABIFChromatogram.create(f);
        default:
            throw new UnsupportedChromatogramFormatException("File "+f+" with magic "+magic+" has an unsupported format");
        }
    }

    /**
     * Creates a new <code>Chromatogram</code> object from the supplied stream.
     * Note that for some chromatogram formats, this can be much more
     * memory-intensive than reading from a file.
     * <p>
     * Note also that if the provided stream is a
     * {@link org.biojava.utils.io.CachingInputStream}, it will be seeked
     * back to 0 before being passed to the parser.  This is because the
     * parsers that use <code>CachingInputStream</code> assume that the
     * "file" starts at 0.
     * </p>
     *
     * @param in the stream from which to read the chromatogram.
     * @return a new Chromatogram
     * @throws IOException when there's a problem with the stream
     * @throws UnsupportedChromatogramFormatException when the file doesn't
     *         contain a chromatogram in a supported format
     */
    public static Chromatogram create(InputStream in)
    throws IOException, UnsupportedChromatogramFormatException {
        CachingInputStream cin;
        if (in instanceof CachingInputStream)
            cin = (CachingInputStream) in;
        else
            cin = new CachingInputStream(in);
        // parsers assume that the image of the file in the stream starts at
        // the beginning of the stream-as-provided.  If the stream
        // was a CachingInputStream, it needs to go to zero.
        cin.seek(0);
        int magic = magicFromStream(cin);
        cin.seek(0);
        switch (magic) {
        case SCF_MAGIC:
            // for SCF, we don't need the cache, so don't use it
            return SCF.create(in, 4);
        case ABI_MAGIC:
            return ABIFChromatogram.create(cin);
        default:
            throw new UnsupportedChromatogramFormatException("The provided input stream with magic "+magic+" has an unsupported format");
        }

    }

  /**
   * Extract the magic number as an integer from a byte-array.
   *
   * <p>
   * This assumes the magic array has at least 4 elements.
   * </p>
   *
   * @param magic  the byte array of magic values
   * @return the magic number integer
   */
    private static int makeMagic(byte[] magic) {
        return (magic[0] << 24) | (magic[1] << 16) | (magic[2] << 8) | (magic[3]);
    }

    /**
     * Reads the next four bytes from a stream to build a 32-bit magic number.
     *
     * @param src the source InputStream
     * @return an integer representing the magic number
     * @throws IOException if data could not be read from src
     */
    private static int magicFromStream(InputStream src) throws IOException {
        byte[] magicBytes = new byte[4];
        src.read(magicBytes);
        return makeMagic(magicBytes);
    }
}
