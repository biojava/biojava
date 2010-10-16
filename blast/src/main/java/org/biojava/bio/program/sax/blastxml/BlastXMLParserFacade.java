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
package org.biojava.bio.program.sax.blastxml;

import java.io.IOException;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParserFactory;

import org.biojava.bio.BioException;
import org.biojava.utils.stax.SAX2StAXAdaptor;
import org.xml.sax.ContentHandler;
import org.xml.sax.DTDHandler;
import org.xml.sax.EntityResolver;
import org.xml.sax.ErrorHandler;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.SAXNotRecognizedException;
import org.xml.sax.SAXNotSupportedException;
import org.xml.sax.XMLReader;
import org.xml.sax.helpers.DefaultHandler;

/**
 * A facade class that wraps the NCBI Blast XML 
 * parsing framework in a more user-friendly form.
 * It is identical to BlastlikeSAXParser in use.
 *
 * @author David Huen
 * @since 1.3
 */
public class BlastXMLParserFacade
    implements XMLReader
{
    // these are defined for the handlers
    // DOWNSTREAM of this one.
    private ContentHandler contentHandler;

    // these are internal handlers
    final BlastXMLParser blasthandler = new BlastXMLParser();

    // this is the SAX parser
    XMLReader parser;

    // this is a default base URI so SAX does not complain
    // when user doesn't give an absolute URI.
    private String baseURI;

    private class Resolver
        implements EntityResolver
    {
        public InputSource resolveEntity(String publicID, String systemID)
            throws SAXException
        {
//            try {
                // resolve the NCBI URN
//                System.out.println("resolve " + publicID + ":" + systemID);

                String resourceName = "org/biojava/bio/program/sax/blastxml/";

                if (publicID.equals("-//NCBI//NCBI BlastOutput/EN")) {
                    resourceName = resourceName + "NCBI_BlastOutput.dtd";
                }
                else if (publicID.equals("-//NCBI//NCBI Entity Module//EN")) {
                    resourceName = resourceName + "NCBI_Entity.mod";
                }
                else if (publicID.equals("-//NCBI//NCBI BlastOutput Module//EN")) {
                    resourceName = resourceName + "NCBI_BlastOutput.mod";
                }
                else
                    return null;

                InputSource is = new InputSource(this.getClass().getClassLoader().getResourceAsStream(resourceName));
                is.setSystemId(baseURI);

                return is;
        }
    }

    public BlastXMLParserFacade()
        throws BioException
    {
        // just initialise content handler
        // to avoid fubar if undefined
        DefaultHandler handler = new DefaultHandler();
        contentHandler = handler;

        try {
            // create the parser
            parser = SAXParserFactory.newInstance().newSAXParser().getXMLReader();

            // assign the EntityResolver
            parser.setEntityResolver(new Resolver());
            parser.setContentHandler(new SAX2StAXAdaptor(blasthandler)); 

            // assign default sane settings
            // namespaces must be true if the SAX2StAX parser isn't to fubar.
            parser.setFeature("http://xml.org/sax/features/namespaces", true);
            parser.setFeature("http://xml.org/sax/features/validation", false);

            // make a base URI just in case the user doesn't
            baseURI = this.getClass().getClassLoader().getResource("org/biojava/bio/program/sax/blastxml/").toString();
        }
        catch (SAXException se) { throw new BioException (se); }
        catch (ParserConfigurationException sce) { throw new BioException(sce); }
    } 

    /**
     * correct this later
     */
    public ContentHandler getContentHandler()
    {
        return contentHandler;
    }

    public DTDHandler getDTDHandler()
    {
        return parser.getDTDHandler();
    }

    /**
     * This class has an EntityResolver that
     * resolves the public ID specifying the
     * NCBI DTDs to resource files within the
     * BioJava libraries.  This call will return
     * that resolver.  It you should set your
     * own resolver, ensure you resolve that
     * URN yourself or the parser will blow up
     * on you!.
     */
    public EntityResolver getEntityResolver()
    {
        return parser.getEntityResolver();
    }

    public ErrorHandler getErrorHandler()
    {
        return parser.getErrorHandler();
    }

    public boolean getFeature(String name)
        throws SAXNotRecognizedException, SAXNotSupportedException
    {
        return parser.getFeature(name);
    }

    public Object getProperty(String name)
        throws SAXNotRecognizedException, SAXNotSupportedException
    {
        return parser.getProperty(name);
    }

    public void parse(InputSource is)
        throws IOException, SAXException
    {
        if (is.getSystemId() == null)
            is.setSystemId(baseURI);
        parser.parse(is);
    }

    public void parse(String systemId)
        throws IOException, SAXException
    {
        parser.parse(systemId);
    }
/**
 * this sets the ContentHandler that receives
 * SAX events from the internal Blast XML parser which
 * is the actual ContentHandler.  <b> It will not
 * change the internal Blast XML parser. </b>
 */
    public void setContentHandler(ContentHandler handler)
    {
        contentHandler = handler;
        blasthandler.setContentHandler(handler);
    }

    public void setDTDHandler(DTDHandler handler)
    {
        parser.setDTDHandler(handler);        
    }

    /**
     * This class has an EntityResolver that
     * resolves the public ID specifying the
     * NCBI DTDs to resource files within the
     * BioJava libraries.  This call will return
     * that resolver.  It you should set your
     * own resolver, ensure you resolve that
     * URN yourself or the parser will blow up
     * on you!.
     */
    public void setEntityResolver(EntityResolver resolver)
    {
        parser.setEntityResolver(resolver);
    }

    public void setErrorHandler(ErrorHandler handler)
    {
        parser.setErrorHandler(handler);
    }

    /**
     * by default, we set the parser to non-validating.
     * change it if you wish/dare!  The parser is 
     * also set to be namespace aware. <b>DO NOT
     * CHANGE THAT!!!</b>
      */
    public void setFeature(String key, boolean value)
        throws SAXNotRecognizedException, SAXNotSupportedException
    {
        parser.setFeature(key, value);
    }

    public void setProperty(String key, Object value)
        throws SAXNotRecognizedException, SAXNotSupportedException
    {
        parser.setProperty(key, value);
    }
}
 
