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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.NoSuchElementException;

import org.biojava.utils.FileAsList;

/**
 * <code>SearchableFileAsList</code> is an abstract base class which
 * implements binary ID searching over the backing random access file.
 *
 * @author Matthew Pocock
 * @author Keith James
 */
abstract class SearchableFileAsList
    extends
        FileAsList
    implements
        SearchableList
{
    /**
     * Creates a new <code>SearchableFileAsList</code> and
     * corresponding backing file.
     *
     * @param file a <code>File</code> used to back the list. This
     * file must not already exist.
     * @param recordLen an <code>int</code> byte record length.
     *
     * @exception IOException if an error occurs.
     */
    public SearchableFileAsList(File file, int recordLen)
        throws IOException {
        super(file, recordLen);
    }

    /**
     * Creates a new <code>SearchableFileAsList</code> instance from
     * an existing backing file.
     *
     * @param file a <code>File</code>  used to back the
     * list. This file must already exist.
     * @param mutable  true if the list can be edited, false otherwise
     *
     * @exception IOException if an error occurs
     */
    public SearchableFileAsList(File file, boolean mutable)
        throws IOException {
        super(file, mutable);
    }

    public Object search(String id) {
        // binary search by id
        byte[] idBytes = id.getBytes();
        byte[] bytes;
    
        int min = 0;
        int max = size()-1;
        do {
            int mid = (min + max) / 2;
      
            bytes = rawGet(mid);
            int cmp = cmp(bytes, idBytes);

            if(cmp < 0) {
                 if(min != mid) {
                     min = mid;
                 } else {
                    min = mid+1;
                 }
            } else if(cmp > 0) {
                 if(max != mid) {
                     max = mid;
                 } else {
                    max = mid-1;
                 }
            } else if(cmp == 0) {
                return parseRecord(bytes);
            }
        } while(min <= max);

        throw new NoSuchElementException("No element with ID: " + id);
    }

    public List searchAll(String id) {
        // binary search by id
        byte[] idBytes = id.getBytes();
        byte[] bytes;

        int min = 0;
        int max = size()-1;
        int mid = -1;
        do {
            mid = (min + max) / 2;

            bytes = rawGet(mid);
            int cmp = cmp(bytes, idBytes);

            if(cmp < 0) {
                if(min != mid) {
                    min = mid;
                } else {
                    min = mid+1;
                }
            } else if(cmp > 0) {
                if(max != mid) {
                    max = mid;
                } else {
                    max = mid-1;
                }
            } else if(cmp == 0) {
                break;
            }
        } while(min <= max);

        if(min > max) {
            throw new NoSuchElementException("No element with ID: " + id);
        }

        ArrayList items = new ArrayList();

        // scan back through file for all items with the same ID
        for(int i = mid-1; i >= 0; i--) {
            bytes = rawGet(i);
            if(cmp(bytes, idBytes) != 0) {
                break;
            }
            items.add(parseRecord(bytes));
        }

        // scan forward through file for all items with the same ID
        for(int i = mid; i < size(); i++) {
            bytes = rawGet(i);
            if(cmp(bytes, idBytes) != 0) {
                System.out.println("Stopped at: " + i);
                break;
            }
            items.add(parseRecord(bytes));
        }

        return items;
    }

    private int cmp(byte[] a, byte[] b) {
        int iMax = Math.min(a.length, b.length);
        for(int i = 0; i < iMax; i++) {
            if(a[i] < b[i]) return -1;
            if(a[i] > b[i]) return +1;
        }

        return 0;
    }
}

