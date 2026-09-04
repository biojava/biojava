/**
 * BioJava development code
 *
 * This code may be freely distributed and modified under the terms of the GNU
 * Lesser General Public Licence. This should be distributed with the code. If
 * you do not have a copy, see:
 *
 * http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual authors. These
 * should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims, or to join the
 * biojava-l mailing list, visit the home page at:
 *
 * http://www.biojava.org/
 */
package org.biojava.nbio.structure.io.density;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.zip.GZIPInputStream;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Recognises CCP4/MRC map files by their header.
 * <p>
 * This is a cheap but effective guard against a cached file that is not
 * actually a map. A misbehaving or overloaded server can answer with an HTML
 * error page and an HTTP 200, in which case nothing about the status code or the
 * content length reveals the problem &mdash; but the missing <code>MAP&nbsp;</code>
 * stamp does, turning a silent failure into a clean cache miss.
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public class Ccp4Header {

	private static final Logger logger = LoggerFactory.getLogger(Ccp4Header.class);

	/**
	 * Byte offset of the four-character format stamp within a CCP4/MRC header. It
	 * sits in word 53 of the 256-word header.
	 */
	public static final int MAP_STAMP_OFFSET = 208;

	/** The stamp itself: the three letters of "MAP" followed by a space. */
	public static final String MAP_STAMP = "MAP ";

	/** Number of header bytes that must be readable for the check to be possible. */
	private static final int HEADER_BYTES = MAP_STAMP_OFFSET + 4;

	private Ccp4Header() {
	}

	/**
	 * Checks whether a file is a CCP4/MRC map. Gzipped files are decompressed on
	 * the fly, so <code>.map.gz</code> works as well as <code>.ccp4</code>.
	 *
	 * @param file the file to check
	 * @return <code>true</code> if the CCP4 stamp is present
	 * @throws IOException if the file could not be read
	 */
	public static boolean isCcp4(File file) throws IOException {
		if (file == null || !file.isFile()) {
			return false;
		}
		// Deliberately no shortcut on file.length(): a gzipped map compresses to far
		// less than the size of the header it contains, so a length test here would
		// reject perfectly good small maps. Reading decides it instead.
		try (InputStream in = openPossiblyGzipped(file)) {
			return isCcp4(in);
		}
	}

	/**
	 * Checks whether a stream carries a CCP4/MRC map header. The stream is read
	 * from its current position and is not closed; it is not un-read afterwards,
	 * so pass a fresh stream or one that supports marking.
	 *
	 * @param in the stream to check, already decompressed
	 * @return <code>true</code> if the CCP4 stamp is present
	 * @throws IOException if the stream could not be read
	 */
	public static boolean isCcp4(InputStream in) throws IOException {
		byte[] header = new byte[HEADER_BYTES];
		int read = 0;
		while (read < HEADER_BYTES) {
			int n = in.read(header, read, HEADER_BYTES - read);
			if (n < 0) {
				return false; // shorter than a CCP4 header, so certainly not one
			}
			read += n;
		}
		String stamp = new String(header, MAP_STAMP_OFFSET, 4, StandardCharsets.US_ASCII);
		return MAP_STAMP.equals(stamp);
	}

	/**
	 * Same as {@link #isCcp4(File)} but reports a problem rather than propagating
	 * it, for use in cache-validity checks where an unreadable file and an invalid
	 * one lead to the same action.
	 *
	 * @param file the file to check
	 * @return <code>true</code> if the file is readable and carries the CCP4 stamp
	 */
	public static boolean isCcp4Quietly(File file) {
		try {
			return isCcp4(file);
		} catch (IOException e) {
			logger.debug("Could not read [{}] to check for a CCP4 header: {}", file, e.getMessage());
			return false;
		}
	}

	private static InputStream openPossiblyGzipped(File file) throws IOException {
		InputStream in = new BufferedInputStream(Files.newInputStream(file.toPath()));
		in.mark(2);
		int b1 = in.read();
		int b2 = in.read();
		in.reset();
		if (b1 == 0x1F && b2 == 0x8B) {
			return new GZIPInputStream(in);
		}
		return in;
	}
}
