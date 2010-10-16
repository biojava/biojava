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
/*
 *                   Eponine development code
 *
 * Eponine is developed by Thomas Down at the Sanger Centre.
 * For more information, see:
 *
 *    http://www.sanger.ac.uk/~td2/eponine/
 *
 * This code is currently not intended for public distribution.
 * If you wish to use or distribute it, or if you have any queries,
 * please contact the author:
 *
 *    td2@sanger.ac.uk
 *
 */

package org.biojava.utils.xml;

/**
 * @author Thomas Down
 * @author Matthew Pocock
 */
public class AppException extends Exception {
    public AppException(String reason) {
	super(reason);
    }
}
