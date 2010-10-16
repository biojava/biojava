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

import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeVetoException;

public abstract class AbstractOrthologueSet 
    extends AbstractChangeable
    implements OrthologueSet
{

    public OrthologueSet filter(OrthologueFilter filter)
    {
        OrthologueSet results = new SimpleOrthologueSet();

        for (Iterator orthoI = iterator();
               orthoI.hasNext(); )
        {
            Orthologue ortho = orthoI.nextOrthologue();

            if (filter.accept(ortho)) {
                try {
                    results.addOrthologue(ortho);
                }
                catch (ChangeVetoException cve) {
                    // should be impossible as this group was created by me
                }
            }
        }
        return results;
    }
}

