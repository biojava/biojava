// BiblioPerson.java
//
//    senger@ebi.ac.uk
//    March 2001
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
package org.biojava.bibliography;

/**
 * <p>
 * It represents a person dealing with the bibliographic resources.
 * </p>
 *
 * @author <A HREF="mailto:senger@ebi.ac.uk">Martin Senger</A>
 * @version $Id$
 * @since 1.3
 */

public class BiblioPerson
    extends BiblioProvider {

  /**
   * The person's sirname (family name).
   */
    public String surname;

  /**
   * The person's christian name (first name).
   */
    public String firstName;

  /**
   * Any middle innitials.
   */
    public String midInitials;

  /**
   * Their e-mail address.
   */
    public String email;

  /**
   * Their postal address.
   */
    public String postalAddress;

  /**
   * Their affiliation.
   */
    public String affiliation;

}
