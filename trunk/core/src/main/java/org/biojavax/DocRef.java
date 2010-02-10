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
 *      http://www.biojava.orDocRef
 */

package org.biojavax;
import java.util.List;

import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Changeable;


/**
 * Represents a documentary reference. Relates to the reference table 
 * in BioSQL.
 * @author Mark Schreiber
 * @author Richard Holland
 * @see RankedDocRef
 * @since 1.5
 */
public interface DocRef extends Comparable,Changeable {
    
    public static final ChangeType CROSSREF = new ChangeType(
            "This reference's crossref has changed",
            "org.biojavax.DocRef",
            "CROSSREF"
            );
    public static final ChangeType REMARK = new ChangeType(
            "This reference's remark has changed",
            "org.biojavax.DocRef",
            "REMARK"
            );
    
    /**
     * The document reference may refer to an object in another database. If so,
     * this method will return that reference.
     * @return Value of property crossref.
     */
    public CrossRef getCrossref();
    
    /**
     * The document reference may refer to an object in another database. Use this
     * method to set that reference. Null will unset it.
     * @param crossref New value of property crossref.
     * @throws ChangeVetoException in case of objections.
     */
    public void setCrossref(CrossRef crossref) throws ChangeVetoException;
    
    /**
     * Returns a textual description of the document reference. This field is 
     * immutable so should be set using the constructor of the implementing class.
     * @return Value of property location.
     */
    public String getLocation();
    
    /**
     * Returns the title of the document reference.
     * @return Value of property title.
     */
    public String getTitle();
    
    /**
     * Returns the authors of the document reference. 
     * It will usually be in the form "Jones H., Bloggs J et al" or similar -
     * a human-readable text value. Editors will have (ed.) appended, 
     * consortiums will have (consortium) appended.
     * @return Value of property authors.
     */
    public String getAuthors();
    
    /**
     * Returns the authors of the document reference as a set of DocRefAuthor
     * implementation instances. This field is immutable so should be set using 
     * the constructor of the implementing class.
     * @return The set of authors.
     */
    public List getAuthorList();
    
    /**
     * Returns a CRC64 checksum of this document reference, allowing for easy
     * comparisons with other document references.
     * @return Value of property CRC.
     */
    public String getCRC();
    
    /**
     * If remarks have been made about this document reference, this method
     * will return them.
     * @return Value of property Remark.
     */
    public String getRemark();
    
    /**
     * Set the remarks for this document reference using this method. Remarks
     * can be anything, it is derived from the equivalent field in the GenBank
     * format.
     * @param Remark New value of property Remark.
     * @throws ChangeVetoException in case of objections.
     */
    public void setRemark(String Remark) throws ChangeVetoException;
    
}
