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
package org.biojava.bio.program.sax;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.HashMap;

import org.xml.sax.Attributes;
import org.xml.sax.ContentHandler;
import org.xml.sax.DTDHandler;
import org.xml.sax.EntityResolver;
import org.xml.sax.ErrorHandler;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.SAXNotRecognizedException;
import org.xml.sax.SAXNotSupportedException;

/**
 * An abstract convenience baseclass implementing the
 * org.xml.sax.XMLReader interface (i.e. SAX2) that application
 * writers can inherit from when writing SAX2 compliant
 * parsers that take non-XML input, and send SAX2 messages
 * to XML ContentHandlers.
 * <p>
 * Note to application writers - support for namespaces may be flaky
 * in some instances. Namespace support is expected to work properly
 * for generating events according to DTDs that use explicit namespaces
 * for all elements.  This is because start/end of Prefix Mapping
 * events are currently not handled properly (at all!).
 * <p>
 * The class provides the following:
 * <ul>
 * <li> An implmentation of the org.xml.sax.XMLReader interface.
 * <li> Utilty methods for generating SAX messages.
 * </ul>
 * <p>
 * The implementation of the org.xml.sax.XMLReader. interface
 * provides working implementation of the following method(s)
 * <ul>
 * <li>setContentHandler(ContentHandler handler)
 * <li>setFeature
 * <li>getFeature
 * </ul>
 * and default, do-nothing implementations of other methods
 * in the interface.
 * <p>
 * Non-interface methods are functional utility methods to help
 * application writers.  For example, see the source of
 * <code>org.biojava.bio.program.sax.BlastLikeSAXParser</code>
 * for an example of how these utility methods can be used.
 * <p>
 * Copyright &copy; 2000 Cambridge Antibody Technology Group plc.
 * <p>
 * Primary author -<ul>
 * <li>Simon Brocklehurst (CAT)
 * </ul>
 * Other authors  -<ul>
 * <li>Tim Dilks          (CAT)
 * <li>Colin Hardman      (CAT)
 * <li>Stuart Johnston    (CAT)
 * <li>Mathieu Wiepert    (Mayo Foundation)
 *</ul>
 *
 * This code was first released to the biojava.org project, July 2000.
 *
 * @author Cambridge Antibody Technology Group plc (CAT)
 * @version 1.1
 *
 * @see BlastLikeSAXParser
 */
abstract class AbstractNativeAppSAXParser
    implements org.xml.sax.XMLReader,
               NamespaceConfigurationIF {

    private   HashMap             oPrefixMap            = new HashMap();
    protected ContentHandler      oHandler              =  null;
    protected boolean             tNamespaces           =  true;
    protected boolean             tNamespacePrefixes    =  false;
    protected String              oNamespacePrefix      = "";
    protected String              oFullNamespacePrefix  = "";

    protected int                 iState;


    /**
     * Allow an application to register a content event handler.
     * If the application does not register a content handler, all content
     * events reported by the SAX parser will be silently ignored.
     * <p>
     * Applications may register a new or different handler in the
     * middle of a parse, and the SAX parser must begin using the
     * new handler immediately.
     *
     * @param poHandler a <code>ContentHandler</code> The XML content handler
     *
     * @throws java.lang.NullPointerException If the handler argument is null
     */
    public void setContentHandler(ContentHandler poHandler) {
    if (poHandler == null) {
        throw new
             java.lang.NullPointerException("ContentHandler is null");
    }
    oHandler = poHandler;
    }
    /**
     * Return the content handler.
     *
     * @return a <code>ContentHandler</code> The current content handler,
     * or null if none has been registered.
     */
    public ContentHandler getContentHandler() {
    return oHandler;
    }
    /**
     * Do-nothing implementation of interface method
     *
     */
    public void parse(InputSource input)
           throws java.io.IOException,
    SAXException {

    }
    /**
     * Full implementation of interface method.
     */
    public void parse(java.lang.String poSystemId)
    throws java.io.IOException,
    SAXException {

    this.parse(new InputSource(poSystemId));
    }
    /**
     * Do-nothing implementation of interface method
     *
     */
    public boolean getFeature(String poName)
    throws SAXNotRecognizedException,SAXNotSupportedException {

    //if get here, throw exception
    throw (new SAXNotSupportedException("The feature \"" + poName +
            "\" is not supported "+
            "in the biojava native SAX2 parsing framework."));
   }
    /**
     * Handles support for ReasoningDomain and Namespace-prefixes
     *
     */
    public void setFeature(java.lang.String poName,
               boolean value)
    throws SAXNotRecognizedException,
    SAXNotSupportedException {

    //handle namespaces
    if (poName.equals("http://xml.org/sax/features/namespaces")) {
        this.setNamespaces(value);
        //check if features combination is supported
        if ((!this.getNamespaces()) && (!this.getNamespacePrefixes())) {
        throw (new SAXNotSupportedException(
               "Illegal feature combination"));
        }

        return;
    }


    if (poName.equals("http://xml.org/sax/features/namespace-prefixes")) {
        this.setNamespacePrefixes(value);
        if ((!this.getNamespaces()) && (!this.getNamespacePrefixes())) {
        throw (new SAXNotSupportedException(
               "Illegal feature combination"));
        }
        return;
    }

    //if get here, throw exception
    throw (new SAXNotSupportedException("The feature \"" + poName +
            "\" is not supported "+
            "in the biojava native SAX2 parsing framework."));
    }
    /**
     * Do-nothing implementation of interface method
     *
     */
    public Object getProperty(String name)
    throws SAXNotRecognizedException,
    SAXNotSupportedException {

    throw (new SAXNotSupportedException("This method is not supported" +
                "in the biojava native SAX2 parser."));
    }

    /**
     * Do-nothing implementation of interface method
     *
     */
    public void setProperty(java.lang.String name,
                        java.lang.Object value)
                 throws SAXNotRecognizedException,
    SAXNotSupportedException {

    throw (new SAXNotSupportedException("This method is not supported"+
                    "in the biojava native SAX2 parser."));
    }


    /**
     * Do-nothing implementation of interface method
     *
     */
    public void setEntityResolver(EntityResolver resolver) {

    }
    /**
     * Do-nothing implementation of interface method
     *
     */
    public EntityResolver getEntityResolver() {
    return null;
    }
    /**
     * Do-nothing implementation of interface method
     *
     */
    public void setDTDHandler(DTDHandler handler) {

    }
    /**
     * Do-nothing implementation of interface method
     *
     */
    public DTDHandler getDTDHandler() {
    return null;
    }

    /**
     * Do-nothing implementation of interface method
     *
     */
    public void setErrorHandler(ErrorHandler handler) {
    }

    /**
     * Do-nothing implementation of interface method
     *
     */
    public ErrorHandler getErrorHandler() {
    return null;
    }
    //----------------------------------------
    /**
     * Utility method to centralize sending of a SAX
     * startElement message to document handler
     *
     * @param poQName a <code>QName</code> value
     * @param atts an <code>Attributes</code> value
     * @exception SAXException if an error occurs
     */
    protected void startElement(QName poQName, Attributes atts)
    throws SAXException{

	oHandler.startElement(poQName.getURI(),
			      poQName.getLocalName(),
			      poQName.getQName(),atts);
    }
    /**
     * Utility method to centralize the sending of a SAX endElement
     * message a document handler.
     *
     * @param poQName   -
     * @exception SAXException thrown if
     * @exception  thrown if
     */
    protected void endElement(QName poQName)
    throws SAXException {

    oHandler.endElement(poQName.getURI(),
                poQName.getLocalName(),
                poQName.getQName());
    }
    /**
     * Utility method to centralize the sending of a SAX characters
     * message a document handler.
     *
     * @param ch   -
     * @param start  -
     * @param length     -
     * @exception SAXException thrown if
     * @exception  thrown if
     */
    protected void characters(char []ch, int start, int length)
    throws SAXException {

    oHandler.characters(ch,start,length);
    }

    /**
     * Support SAX2 configuration of namespace support of parser.
     */
    private void setNamespaces(boolean ptNamespaces) {
    tNamespaces = ptNamespaces;
    }

    /**
     * Support SAX2 configuration of namespace support of parser.
     */
    public boolean getNamespaces() {
    return tNamespaces;
    }
    /**
     * Support SAX2 configuration of namespace support of parser.
     */
    private void setNamespacePrefixes(boolean ptNamespacePrefixes) {
    tNamespacePrefixes = ptNamespacePrefixes;
    }
    /**
     * Support SAX2 configuration of namespace support of parser.
     */
    public boolean getNamespacePrefixes() {
    return tNamespacePrefixes;
    }
    /**
     * Adds a namespace prefix to URI mapping as (key,value) pairs.
     * This mapping can be looked up later to get URIs on request
     * using the getURIFromPrefix method.
     *
     * @param poPrefix a <code>String</code> representation of the
     * namespace prefix
     * @param poURI a <code>String</code> representation of the URI
     * for the namespace prefix.
     *
     */
    public void addPrefixMapping(String poPrefix, String poURI) {

    oPrefixMap.put(poPrefix,poURI);
    }

    /**
     * Gets the URI for a namespace prefix, given that prefix,
     * or null if the prefix is not recognised.
     *
     * @param poPrefix a <code>String</code> The namespace prefix.
     */
    public String getURIFromPrefix(String poPrefix) {

    if (oPrefixMap.containsKey(poPrefix)) {
        return (String)oPrefixMap.get(poPrefix);
    } else {
        //return an empty string, null makes XSLT processors unhappy
        return "";
    }
    }
    /**
     *
     * @param poPrefix a <code>String</code> value
     */
    public void setNamespacePrefix(String poPrefix) {
    oNamespacePrefix = poPrefix;
    oFullNamespacePrefix = oNamespacePrefix.concat(":");
    }
    /**
     * Describe <code>getNamespacePrefix</code> method here.
     *
     * @return a <code>String</code> value
     */
    public String getNamespacePrefix() {
    return oNamespacePrefix;
    }
    /**
     * Given an unprefixed element name, returns
     * a new element name with a namespace prefix
     *
     * @return a <code>String</code> value
     */
    public String prefix(String poElementName) {
    return oFullNamespacePrefix.concat(poElementName);
    }

    /**
     * Create a stream from an an InputSource, picking the
     * correct stream according to order of precedance.
     *
     * @param poSource an <code>InputSource</code> value
     * @return a <code>BufferedReader</code> value
     */
    protected BufferedReader getContentStream(InputSource poSource) {


    InputSource               oSource;
    BufferedReader            oContents;
    int                       iBuffSize          = 8192;


    oSource = poSource;

    //Check contents InputSource in order of precedence

    //Highest - Character stream

    if (oSource.getCharacterStream() != null) {

        oContents = new BufferedReader(oSource.getCharacterStream(),
                       iBuffSize);
        return oContents;
    }

    //Next to lowest -  Byte stream


    if ( (oSource.getByteStream() != null) ) {

        oContents = new BufferedReader(
           new InputStreamReader(oSource.getByteStream()),
                     iBuffSize);

        return oContents;
    }

    //Lowest precedence - System URI

    if ( (oSource.getSystemId() != null) ) {

        try {
        //URL handles file URI's and URL's, which might be sent in
        //by an XSLT processor.  Handles things like file:\C:\file.txt
        URL IORUrl = new java.net.URL(oSource.getSystemId());
        BufferedInputStream inStr =
                new BufferedInputStream(IORUrl.openStream());
        // Create input stream
        oContents = new
            BufferedReader(new InputStreamReader(inStr),
                   iBuffSize);

        return oContents;
        //catch Malfromed URL and IOException
        } catch (MalformedURLException x) {
        System.out.println(x.getMessage());
        System.out.println("Couldn't open file");
        System.exit(0);
        } catch (IOException x) {
        System.out.println(x.getMessage());
        System.out.println("Couldn't open file");
        System.exit(0);
        }
    }

    //should throw an exception to say couldn't get stream
    return null;
    }

    /**
     * Centralise chaining of iState field to help
     * with debugging. E.g. printing out value etc.
     * All changes to iState should be made through this method.
     *
     * @param piState an <code>int</code> value
     */
    protected void changeState(int piState) {
    iState = piState;
    }



}
