/*
 *                  BioJava development code
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
 * Created on Mar 9, 2005
 *
 */
package org.biojava.dasobert.das;

import org.biojava.utils.stax.StAXContentHandlerBase;
import org.biojava.utils.stax.*;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.biojava.bio.program.das.dasstructure.*;
import org.biojava.bio.structure.Structure;

/** a Class to adapt the Biojava StructureXMLResponseParser class for usage with StAX parser.
 * shouldbe moved to biojava eventually ...
 * 
 * @author Andreas Prlic
 *
 */
public class StructureXMLStAXAdaptor extends StAXContentHandlerBase {
    
    DASStructureXMLResponseParser parser;
    /**
     * 
     */
    public StructureXMLStAXAdaptor() {
        super();
        parser = new DASStructureXMLResponseParser();
        // TODO Auto-generated constructor stub
    }
    
    public void startElement(String nsURI,
            String localName,
            String qName,
            Attributes attrs,
            DelegationManager dm)
    throws SAXException
    {
        //System.out.println("startElement PDB "+nsURI + " localName "+ qName+ " attrs ");
        parser.startElement(nsURI, localName,qName,attrs);
    }
    
    public void endElement(String nsURI,
            String localName,
            String qName,
            StAXContentHandler delegate)
    throws SAXException
    {
        parser.endElement(nsURI,localName,qName);
    }
    
    public void characters(char[] ch,
            int start,
            int length)
    throws SAXException
    {
        parser.characters(ch,start,length);
    }
    
    public void startTree()
    throws SAXException
    {
        System.out.println("Structure StAX adaptor called");
        parser.startDocument();
    }
    
    public void endTree()
    throws SAXException
    {
        parser.endDocument();
    }
    
    
    public Structure getStructure() {
        // wrong name of function, sould change this at some point ...
        return parser.get_structure();
    }
    
}
