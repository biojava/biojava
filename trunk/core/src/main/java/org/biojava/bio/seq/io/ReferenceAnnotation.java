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
 */

package org.biojava.bio.seq.io;

import java.util.HashMap;
import java.util.Map;

import org.biojava.bio.AbstractAnnotation;
import org.biojava.utils.ChangeVetoException;

/**
 * @author Lorna Morris
 * @since 1.3.1
 * @deprecated Use org.biojavax.bio.seq.io framework instead
 */
public class ReferenceAnnotation extends AbstractAnnotation {

     /**
   * The properties map. This may be null if no property values have
   * yet been set.
   */
    private Map properties;

    public ReferenceAnnotation() {

            super();
        try {
            //System.out.println("Calling refAnnot");
            this.setProperty(EmblLikeFormat.SEPARATOR_TAG, "");//all references have an epty XX line
        } catch (ChangeVetoException e) {
            e.printStackTrace();
        }
    }

    protected Map getProperties() {
        if(!propertiesAllocated()) {
            properties = new HashMap();
        }
        return properties;
    }

    protected boolean propertiesAllocated() {
        return properties != null;
    }
}
