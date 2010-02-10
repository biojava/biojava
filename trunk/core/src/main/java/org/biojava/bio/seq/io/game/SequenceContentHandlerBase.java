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

package org.biojava.bio.seq.io.game;

import org.biojava.bio.seq.io.StreamParser;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.utils.stax.DelegationManager;
import org.biojava.utils.stax.StAXContentHandler;
import org.biojava.utils.stax.StAXContentHandlerBase;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * StAX handler for elements containing sequence
 *
 * @author David Huen
 * @since 1.8
 */
public class SequenceContentHandlerBase extends StAXContentHandlerBase 
{
    private int level = 0;
    private StreamParser parser;

    public void startElement(String nsURI,
			     String localName,
			     String qName,
			     Attributes attrs,
			     DelegationManager dm)
	 throws SAXException
    {
	level++;
	if (level > 1) {
	    throw new SAXException("Found child element when expecting character data");
	}
    }

    public void endElement(String nsURI,
			   String localName,
			   String qName,
			   StAXContentHandler handler)
	throws SAXException
    {
	level--;
    }

/**
 * assign a StreamParser object to instance.
 */
    public void setStreamParser(StreamParser parser)
    {
       this.parser = parser;
    }

    public void characters(char[] ch, int start, int length) 
        throws SAXException
    {
       int cnt=0;
       int bstart = start;
       int bcnt=0;
       // call parser to pass Symbols to SeqIOListener
       try {
//         System.out.println("SequenceContentHandlerBase: " + start + " " + length);
         while (cnt < length) {
           // clear whitespace
           while (cnt < length && (!(Character.isLetter(ch[start + cnt]))) ) cnt++;

           // map length of non-whitespace
           bstart = start + cnt; bcnt = 0;
           while (cnt < length && (Character.isLetter(ch[start + cnt])) ) {
             cnt++;
             bcnt++;
           }

           // process current block
           parser.characters(ch, bstart, bcnt);
         }
         parser.close();
         
       }
       catch (IllegalSymbolException ise) {
         throw new SAXException("SequenceContentHandlerBase: illegal symbol encountered.");
       }
    }
}
