/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.survival.data;

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

    public char charAt(int index) {
        int ix = index + offset;
        if (ix >= end) {
            throw new StringIndexOutOfBoundsException("Invalid index "
                    + index + " length " + length());
        }
        return (char) (data[ix] & 0xff);
    }

    public int length() {
        return end - offset;
    }

    public CharSequence subSequence(int start, int end) {
        if (start < 0 || end >= (this.end - offset)) {
            throw new IllegalArgumentException("Illegal range "
                    + start + "-" + end + " for sequence of length " + length());
        }
        return new CompactCharSequence(data, start + offset, end + offset);
    }

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
