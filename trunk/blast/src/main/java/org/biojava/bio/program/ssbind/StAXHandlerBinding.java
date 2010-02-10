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

import org.biojava.bio.program.xff.ElementRecognizer;

/**
 * <code>StAXHandlerBinding</code>s associates an
 * <code>ElementRecognizer</code> with a factory which creates
 * <code>StAXContentHandler</code>s for elements which it the
 * <code>ElementRecognizer</code> accepts.
 *
 * @author Thomas Down
 * @author Keith James
 * @since 1.2
 */
public class StAXHandlerBinding
{
    ElementRecognizer recognizer;
    StAXHandlerFactory factory;

    /**
     * Creates a new <code>StAXHandlerBinding</code>.
     *
     * @param recognizer an <code>ElementRecognizer</code>.
     * @param factory a <code>StAXHandlerFactory</code>.
     */
    StAXHandlerBinding(ElementRecognizer recognizer,
                       StAXHandlerFactory factory)
    {
        this.recognizer = recognizer;
        this.factory    = factory;
    }
}
