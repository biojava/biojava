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

package org.biojava.bio.seq.io;

/**
 * This class is necessitated by the deprecation of
 * writeSequence(Sequence seq, String format, PrintStream os)
 * method.  Now the format object must be intimately
 * associated to the format it is handling.
 * 
 * @author David Huen
 * @deprecated Use org.biojavax.bio.seq.io framework instead
 */
public class GenpeptFormat
    extends GenbankFormat
{

    public String getDefaultFormat()
    {
        return "GENPEPT";
    }


}

