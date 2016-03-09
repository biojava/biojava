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
package org.biojava.nbio.survival.data;

import java.io.Serializable;
import java.io.UnsupportedEncodingException;

/**
 *http://www.javamex.com/tutorials/memory/ascii_charsequence.shtml
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class CompactCharSequence implements CharSequence, Serializable {

	static final long serialVersionUID = 1L;
	private static final String ENCODING = "ISO-8859-1";
	private final int offset;
	private final int end;
	private final byte[] data;
	private final boolean nullstring;

	private CompactCharSequence(byte[] data, int offset, int end) {
		this.data = data;
		this.offset = offset;
		this.end = end;
		nullstring = false;
	}

	/**
	 *
	 * @param str
	 */
	public CompactCharSequence(String str) {
		try {
			if (str != null) {
				data = str.getBytes(ENCODING);
				offset = 0;
				end = data.length;
				nullstring = false;
			} else {
				data = new byte[0];
				offset = 0;
				end = 0;
				nullstring = true;
			}
		} catch (UnsupportedEncodingException e) {
			throw new RuntimeException("Unexpected: " + ENCODING + " not supported!");
		}
	}

	@Override
	public char charAt(int index) {
		int ix = index + offset;
		if (ix >= end) {
			throw new StringIndexOutOfBoundsException("Invalid index "
					+ index + " length " + length());
		}
		return (char) (data[ix] & 0xff);
	}

	@Override
	public int length() {
		return end - offset;
	}

	@Override
	public CharSequence subSequence(int start, int end) {
		if (start < 0 || end >= (this.end - offset)) {
			throw new IllegalArgumentException("Illegal range "
					+ start + "-" + end + " for sequence of length " + length());
		}
		return new CompactCharSequence(data, start + offset, end + offset);
	}

	@Override
	public String toString() {
		try {
			if (nullstring) {
				return null;
			} else {
				if(length() == 0)
					return "";
				else
					return new String(data, offset, end - offset, ENCODING);
			}
		} catch (UnsupportedEncodingException e) {
			throw new RuntimeException("Unexpected: " + ENCODING + " not supported");
		}
	}
}
