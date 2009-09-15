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
import java.util.ArrayList;
import java.util.List;

import org.xml.sax.Attributes;
import org.xml.sax.ContentHandler;
import org.xml.sax.Locator;
import org.xml.sax.SAXException;
/**
 * Lightweight adaptor which translates SAX content events into
 * StAX form, and provides delegation services.
 *
 * @author Thomas Down
 * @author modified by Hanning Ni
 */
public class SAX2StAXAdaptor implements ContentHandler {
    private List stack;
    private HandlerBinding current;
    private boolean isRecursive ;

    {
        stack = new ArrayList();
        isRecursive = true ;
    }

    /**
     * Construct a new SAX Content handler which wraps a StAX
     * handler.
     */

    public SAX2StAXAdaptor(StAXContentHandler rootHandler) {
        current = new HandlerBinding(rootHandler);
        stack.add(current);
    }

    public void startDocument() throws SAXException {
        current.handler.startTree();
    }

    public void endDocument() throws SAXException {
        current.handler.endTree();
    }

    public void characters(char[] ch,
                           int start,
                           int end)
        throws SAXException
    {
        current.handler.characters(ch, start, end);
    }

    public void ignorableWhitespace(char[] ch,
                                    int start,
                                    int end)
        throws SAXException
    {
        current.handler.ignorableWhitespace(ch, start, end);
    }

    public void startPrefixMapping(String prefix, String uri)
        throws SAXException
    {
        current.handler.startPrefixMapping(prefix, uri);
    }

    public void endPrefixMapping(String prefix)
        throws SAXException
    {
        current.handler.endPrefixMapping(prefix);
    }

    public void processingInstruction(String target, String data)
        throws SAXException
    {
        current.handler.processingInstruction(target, data);
    }

    public void setDocumentLocator(Locator locator) {
        current.handler.setDocumentLocator(locator);
    }

    public void skippedEntity(String name)
        throws SAXException
    {
        current.handler.skippedEntity(name);
    }

    public void startElement(final String nsURI,
                             final String localName,
                             final String qName,
                             final Attributes attrs)
        throws SAXException
    {
        S2SDelegationManager dm = new S2SDelegationManager();
        dm.setRecursive( isRecursive ) ;
        current.handler.startElement(nsURI,
                                     localName,
                                     qName,
                                     attrs,
                                     dm);

        if (dm.getDelegate() != null) {
            current = new HandlerBinding(dm.getDelegate());
            stack.add(current);
            isRecursive = false ;
            current.handler.startTree();
            startElement(nsURI, localName, qName, attrs); // Recurse until someone
        } else {
            current.count++;
            isRecursive = true ;
        }
    }

    private static class S2SDelegationManager implements DelegationManager {
        private StAXContentHandler delegate;
        private boolean recursive = true ;

        public void delegate(StAXContentHandler handler)
            throws SAXException
        {
            if (this.delegate != null) {
                throw new SAXException("Tried to multiply delegate a single StAX element");
            }

            this.delegate = handler;
        }

        public void setRecursive(boolean recursive )
        {
            this.recursive = recursive ;
        }

        public boolean getRecursive()
        {
            return recursive ;
        }

        private StAXContentHandler getDelegate() {
            return delegate;
        }
    }

    public void endElement(String nsURI, String localName, String qName)
        throws SAXException
    {
        current.handler.endElement(nsURI, localName, qName, null);
        current.count--;
        while (current.count == 0 && stack.size() > 1) {
            current.handler.endTree();
            stack.remove(stack.size() - 1);
            current = (HandlerBinding) stack.get(stack.size() - 1);
         //   current.handler.endElement(nsURI, localName, qName, oldHandler);
        }
    }

    private class HandlerBinding {
        public StAXContentHandler handler;
        public int count;

        private HandlerBinding(StAXContentHandler handler) {
            this.handler = handler;
            this.count = 0;
        }
    }
}
