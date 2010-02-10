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

package org.biojava.utils.net;

import java.net.URL;

/**
 * <p><code>URLFactory</code> defines a means of obtaining a URL
 * associated with an object. The URL returned may be based on any
 * property of the object, for example its Java class, methods or
 * fields or its <code>Annotation</code>. As the criteria by which the
 * URL are created will be highly variable it is left to the
 * implementation to cast the <code>Object</code> argument and perform
 * any necessary checks. An implementation may make any additional
 * checks such as applying <code>PropertyConstraint</code>s or
 * checking an <code>AnnotationType</code>.</p>
 *
 * <p>An example use case is in obtaining hyperlink target to
 * associate with a sequence hit in a database search which will then
 * be placed in an image map.</p>
 *
 * @author <a href="mailto:kdj@sanger.ac.uk">Keith James</a>
 */
public interface URLFactory
{
    /**
     * <code>createURL</code> returns a URL which is relevant to the
     * object in a way specified by the implementation.
     *
     * @param object an <code>Object</code>.
     *
     * @return a <code>URL</code>.
     */
    public URL createURL(Object object);
}
