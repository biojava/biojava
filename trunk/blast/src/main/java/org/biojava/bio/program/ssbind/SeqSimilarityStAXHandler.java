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

import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.program.xff.ElementRecognizer;
import org.biojava.utils.stax.DelegationManager;
import org.biojava.utils.stax.StAXContentHandler;
import org.biojava.utils.stax.StAXContentHandlerBase;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * <code>SeqSimilarityStAXHandler</code> is a base class for creating
 * modular StAX handlers which send callbacks to a
 * <code>SeqSimilarityStAXAdapter</code>.
 *
 * @author Keith James
 * @since 1.3
 */
public class SeqSimilarityStAXHandler extends StAXContentHandlerBase
{
    // Available handler bindings
    private List bindings;
    // Incremented on startElement, decremented on endElement. Used to
    // identify which method calls to handle here and which to
    // delegate.
    private int depth;

    // The target handler
    protected SeqSimilarityStAXAdapter ssContext;

    /**
     * Creates a new <code>SeqSimilarityStAXHandler</code> which
     * simply maintains a list of <code>StAXHandlerBinding</code>s and
     * delegates to any suitable <code>StAXContentHandler</code> bound
     * by one of them.
     */
    public SeqSimilarityStAXHandler(SeqSimilarityStAXAdapter ssContext)
    {
        this.ssContext = ssContext;
        bindings = new ArrayList();
    }

    public void startElement(String            nsURI,
                             String            localName,
                             String            qName,
                             Attributes        attrs,
                             DelegationManager dm)
        throws SAXException
    {
        depth++;

        if (depth == 1)
        {
            handleStartElement(nsURI, localName, qName, attrs);
        }
        else
        {
            for (int i = bindings.size(); --i >= 0;)
            {
                StAXHandlerBinding b = (StAXHandlerBinding) bindings.get(i);
            
                if (b.recognizer.filterStartElement(nsURI, localName, qName, attrs))
                {
                    dm.delegate(b.factory.getHandler(ssContext));
                    return;
                }
            }
        }
    }

    public void endElement(String             nsURI,
                           String             localName,
                           String             qName,
                           StAXContentHandler handler)
        throws SAXException
    {
        depth--;

        if (depth == 0)
        {
            handleEndElement(nsURI, localName, qName);
        }
    }

    protected void addHandler(ElementRecognizer recognizer,
                              StAXHandlerFactory factory)
    {
        bindings.add(new StAXHandlerBinding(recognizer, factory));
    }

    protected void handleStartElement(String     nsURI,
                                      String     localName,
                                      String     qName,
                                      Attributes attrs)
        throws SAXException { }

    protected void handleEndElement(String nsURI,
                                    String localName,
                                    String qName)
        throws SAXException { }
}
