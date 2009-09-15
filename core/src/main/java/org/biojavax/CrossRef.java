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

package org.biojavax;

import org.biojava.utils.Changeable;

/**
 * Represents a cross reference to another database. Relates to the 
 * dbxref table in BioSQL. The interface is immutable, the fields are
 * intended to be set by the constructor.
 * @author Mark Schreiber
 * @author Richard Holland
 * @see RankedCrossRef
 * @since 1.5
 */
public interface CrossRef extends RichAnnotatable,Comparable,Changeable {
    
   /**
     * Returns the name of the database the cross reference refers to. This
     * would normally be a namespace name, eg. 'gb' for GenBank.
     * @return Value of property dbname.
     */
    public String getDbname();

    /**
     * Returns the accession of the object that the crossref refers to.
     * @return Value of property accession.
     */
    public String getAccession();

    /**
     * Returns the version of the object that the crossref refers to.
     * @return Value of property version.
     */
    public int getVersion();
    
}

