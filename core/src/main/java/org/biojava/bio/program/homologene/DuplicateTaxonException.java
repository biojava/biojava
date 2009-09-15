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

package org.biojava.bio.program.homologene;

import org.biojava.bio.BioException;

public class DuplicateTaxonException
    extends BioException
{

    public DuplicateTaxonException()
    {
        super();
    }

    public DuplicateTaxonException(String message)
    {
        super(message);
    }

    public DuplicateTaxonException(Throwable t)
    {
        super(t);
    }

    public DuplicateTaxonException(Throwable t, String message)
    {
        super(message, t);
    }
}
