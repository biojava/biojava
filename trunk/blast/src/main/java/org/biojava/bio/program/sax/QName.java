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

import java.util.StringTokenizer;

/**
 * A utility class with package level visibility for dealing 
 * element and attribute names that are potentially
 * namespace-qualified - qNames e.g. <code>biojava:MyElement</code>
 * 
 * It takes a qName, and allows the SAXParser applications writers to
 * obtain prefix and local names.  Note, whether these
 * values are reported or not, or reported as unknown, should be
 * controlled by the namespace and namespace-prefix properties
 * of the parser.
 *
 * For attributes, the LocalName will be reported as an empty string
 * if the attribute should not be reported by the parser due to
 * a namespace configuration.
 *
 * <p>
 * Copyright &copy; 2000 Cambridge Antibody Technology Group plc.
 * 
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
 */
final class QName  {

    private String                    oLocalName          = "";      
    private String                    oPrefixName         = "";     
    private String                    oQName              = "";     
    private String                    oURI                = "";     
    private NamespaceConfigurationIF  oNamespace;
    private StringTokenizer           oTokenizer;
    /**
     * Creates a new <code>QName</code> instance.
     *
     * @param poNamespace a <code>NamespaceConfigurationIF</code> value
     */
    public QName(NamespaceConfigurationIF poNamespace) {
    oNamespace = poNamespace;
    }
    /**
     * Sets the qName for the object. When the qName is set,
     * local names and namespace prefixes are automatically
     * accessible from calling relevant methods on th object.
     *
     * @param poQName a <code>String</code> value of qName. For
     * example, it could be &quot;biojava:MyElement&quot;
     */
    
    public QName(NamespaceConfigurationIF poNamespace,String poQName) {
    oNamespace = poNamespace;
    this.setQName(poQName);
    }
    /**
     * Sets the qName for the object. When the qName is set,
     * local names and namespace prefixes are automatically
     * accessible from calling relevant methods on the object.
     *
     * @param poQName a <code>String</code> value of qName. For
     * example, it could be &quot;biojava:MyElement&quot;
     */
    public void setQName(String poQName) {
    this.reset();
    oQName = poQName;
    if (poQName.indexOf(':') != -1) {
        //if name contains a colon, get prefix and local name
        oTokenizer  = new StringTokenizer(oQName,":");
        oPrefixName = oTokenizer.nextToken();
        if ( (oPrefixName.equals("xmlns")) && 
          (!oNamespace.getNamespacePrefixes()) ) {  
            //prefix != "xmlns" and getNameSpacePrefixes=true
            //leave LocalName as empty string
        } else {
            //special case for xml:space 
            oLocalName  = oTokenizer.nextToken();
        }
        if (oPrefixName.equals("xml")) {
            //Set up namespace correctly for XSLT
            oURI = "http://www.w3.org/XML/1998/namespace"; 
        }
        else if (oPrefixName.equals("xmlns")) {
            //means that oNamespace.getNamespacePrefixes() = true
            if (oTokenizer.hasMoreTokens()) {
                if (oTokenizer.nextToken().equals("biojava")) {
                    //looking for xmlns:biojava and processing 
                    //namespace prefixes
                    oURI = oNamespace.getURIFromPrefix("biojava");
                }
            }
        } else {
            oURI = oNamespace.getURIFromPrefix(oPrefixName);
        }
    } else {
        //if name or attribute does not have a namespace,
        //local name is QName.
        oPrefixName = "";
        if ( (oQName.equals("xmlns")) &&
          (!oNamespace.getNamespacePrefixes()) ) {
            //leave LocalName as empty string
        } else {
            oLocalName =  oQName;
        }
    }
    }
    /**
     * Gets the raw QName - Parser writers should
     * determine if the parser should report it or not
     *
     * @return a <code>String</code> value
     */
    public String  getQName() {
    return oQName;
    }
    /**
     * Gets the raw Prefix - Parser writers should
     * determine if the parser should report it or not
     * (and whether it should be used to get the namespace
     * URI from the parser). 
     *
     * @return a <code>String</code> value
     */
    public String  getPrefixName() {
    return oPrefixName;
    }
    /**
     * Gets the raw LocalName - Parser writers should
     * determine if the parser should report it or not.
     *
     * @return a <code>String</code> value
     */
    public String  getLocalName() {
    return oLocalName;
    }


    public String getURI() {
    return oURI;
    }

    /**
     * Reset private fields to allow partial object reuse
     *
     */
    private void reset() {
    oLocalName          = "";      
    oPrefixName         = "";     
    oURI                = "";     
    }
}
