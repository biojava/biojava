// CandyEntry.java
//
//    senger@ebi.ac.uk
//    February 2001
//

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
package org.biojava.utils.candy;

import java.util.Hashtable;

/**
 * <p>
 * This is a basic container for a vocabulary entry. It consists only
 * of the basic attributes which is sufficient for the vocabularies
 * providing string-type contents.
 * </p>
 *
 * <p>
 * However, it may still accomodate more complex data types using
 * the {@link #extras} member.
 * </p>
 * 
 * @author <A HREF="mailto:senger@ebi.ac.uk">Martin Senger</A>
 * @version $Id$
 */

public class CandyEntry {

    /** A unique identifier of this entry. */
    public String entry = "";

    /** A value of this entry. */
    public String description = "";

    /** A container for the additional properties represented by this entry. */
    public Hashtable extras = new Hashtable();

    /** An empty constructor. */
    public CandyEntry() {
    }

    /** It creates an entry instance with given name and empty value. */
    public CandyEntry (String entry) {
	this (entry, "", null);
    }

    /** It creates an entry instance with given name and value. */
    public CandyEntry (String entry, String description) {
        this (entry, description, null);
    }

    /**
     * It creates an entry instance with given name, value and
     * additional properties.
     */
    public CandyEntry (String entry, String description, Hashtable extras) {
	this.entry = entry;
	this.description = description;
	if (extras != null)
	    this.extras = extras;
    }

    /** It prints the entry contents. */
    public String toString() {
	return entry + "\t" + description +
	    (extras != null && extras.size() > 0 ?
	     "\n\t" + extras.toString() : "");
    }
}
