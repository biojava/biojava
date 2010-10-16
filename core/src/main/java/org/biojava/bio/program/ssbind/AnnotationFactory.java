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

package org.biojava.bio.program.ssbind;

import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.SimpleAnnotation;
import org.biojava.bio.SmallAnnotation;
import org.biojava.utils.ChangeVetoException;

/**
 * <code>AnnotationFactory</code> is a utility class for making
 * <code>Annotation</code>s from <code>Map</code>s. Shared by the
 * search and homology builders. Now public to allow use by anyone
 * making custom handlers.
 *
 * @author Keith James
 * @author Thomas Down
 * @since 1.2
 */
public class AnnotationFactory
{
    /**
     * <code>makeAnnotation</code> creates the annotation.
     *
     * @param m a <code>Map</code> of raw data.
     * @return an <code>Annotation</code>.
     */
    public static Annotation makeAnnotation(Map m)
    {
        int elements = m.size();

        if (elements == 0)
        {
            return Annotation.EMPTY_ANNOTATION;
        }
        else
        {
            Annotation annotation;

            if (elements < 15)
                annotation = new SmallAnnotation();
            else
                annotation = new SimpleAnnotation();

            Set keySet = m.keySet();

            try
            {
                for (Iterator ksi = keySet.iterator(); ksi.hasNext();)
                {
                    Object key = ksi.next();
                    annotation.setProperty(key, m.get(key));
                }
            }
            catch (ChangeVetoException cve)
            {
                throw new BioError("Assert failed: couldn't modify newly created Annotation",cve);
            }

            return annotation;
        }
    }
}
