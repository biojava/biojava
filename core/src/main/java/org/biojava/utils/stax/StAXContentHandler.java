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
package org.biojava.utils.stax;

import org.xml.sax.Attributes;
import org.xml.sax.Locator;
import org.xml.sax.SAXException;

/**
 * Interface for StAX content handlers.  This interface is
 * very similar in spirit and design to the SAX content handler.
 * Differences are:
 *
 * <ol>
 * <li>start/endDocument methods are replaced by start/endTree.  This
 *     recognises the fact that a StAX content handler may only see
 *     a sub-tree of an XML document, rather than the whole document.</li>
 * <li>the startElement method takes a <code>DelegationManager</code>,
 *     allowing delegation of sub-trees to other content handlers.</li>
 * </ol>
 *
 * @author Thomas Down
 */

public interface StAXContentHandler {
    public void startTree()
	throws SAXException;

    public void endTree()
	throws SAXException;

    public void characters(char[] ch,
			   int start,
			   int length)
	throws SAXException;

    public void ignorableWhitespace(char[] ch,
				    int start,
				    int length)
	throws SAXException;

    public void startPrefixMapping(String prefix, String uri)
        throws SAXException;

    public void endPrefixMapping(String prefix)
        throws SAXException;

    public void processingInstruction(String target, String data)
        throws SAXException;

    public void setDocumentLocator(Locator locator);

    public void skippedEntity(String name)
        throws SAXException;

    public void startElement(String nsURI,
			     String localName,
			     String qName,
			     Attributes attrs,
			     DelegationManager dm)
	throws SAXException;

    public void endElement(String nsURI,
			   String localName,
			   String qName,
			   StAXContentHandler delegate)
        throws SAXException;
}
