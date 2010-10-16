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
package org.biojava.bio.seq.io.agave;
import org.xml.sax.Attributes;
import org.xml.sax.Locator;
import org.xml.sax.SAXException;

/**
 * Simple implementation of the <code>StAXContentHandler</code>
 * interface, with empty implementations for all the methods.
 *
 * <p>
 * This class is provided as a base for content handlers where
 * the implementor does not wish to provide all the methods.
 * </p>
 *
 * @author copied from Thomas Down
 */
public class StAXContentHandlerBase implements StAXContentHandler {
    public void startTree()
        throws SAXException
    {
    }

    public void endTree()
        throws SAXException
    {
    }

    /**
     * Signal a span of character data in the XML input.
     *
     * @param ch an array of characters
     * @param start index of the first significant character for this event.
     * @param length number of characters significant to this event.
     */

    public void characters(char[] ch,
                           int start,
                           int length)
        throws SAXException
    {
    }

    public void ignorableWhitespace(char[] ch,
                                    int start,
                                    int length)
        throws SAXException
    {
    }

    public void startPrefixMapping(String prefix, String uri)
        throws SAXException
    {
    }

    public void endPrefixMapping(String prefix)
        throws SAXException
    {
    }

    public void processingInstruction(String target, String data)
        throws SAXException
    {
    }

    public void setDocumentLocator(Locator locator)
    {
    }

    public void skippedEntity(String name)
        throws SAXException
    {
    }

    public void startElement(String nsURI,
                             String localName,
                             String qName,
                             Attributes attrs,
                             DelegationManager dm)
        throws SAXException
    {
    }

    public void endElement(String nsURI,
                           String localName,
                           String qName,
                           StAXContentHandler delegate)
        throws SAXException
    {
    }
}
