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
package org.biojava.bio.program.xml;

import java.util.Stack;

import org.xml.sax.Attributes;

/**
 *
 * Base XMLWriter class for writing XML representations of Java Value
 * Objects with bespoke architectures.
 *
 * Currently reporting QNames for all attribute values as a fix
 * for needing to defining namespaces by attributes with QNames.
 * This may not be ideal.
 *
 * <p>
 * Copyright &copy; 2000 Cambridge Antibody Technology.
 *
 * <p>
 * Primary author -<ul>
 * <li>Simon Brocklehurst (CAT)
 * </ul>
 * Other authors  -<ul>
 * <li>Tim Dilks          (CAT)
 * <li>Colin Hardman      (CAT)
 * <li>Stuart Johnston    (CAT)
 *</ul>
 *
 * @author Cambridge Antibody Technology (CAT)
 * @author Greg Cox
 * @version 1.01
 *
 */
public class BaseXMLWriter {
    private Stack oElementStack;
    private StringBuffer oStr;
    private StringBuffer oIndent;
    private String oLineSeparator;
    private StringBuffer oSpaceStore;
    private boolean tLastWritePCData;
    private boolean tLastWriteStartElement;
    private String oAttName;
    private String oAttValue;
    static final String oHeader = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>";

    public BaseXMLWriter()
    {
        oElementStack = new Stack();
        oStr = new StringBuffer();
        oIndent = new StringBuffer();
        oSpaceStore = new StringBuffer();
        oLineSeparator = System.getProperty("line.separator");
        oIndent.setLength(0);
        tLastWritePCData = false;
	tLastWriteStartElement = false;
    }

    private void decreaseIndent()
    {
        oIndent.setLength(oIndent.length() - 2);
    }

    public String endElement()
    {
        decreaseIndent();
        oStr.setLength(0);
        if (!tLastWriteStartElement && !tLastWritePCData)
        {
            oStr.append(nl());
            oStr.append(indent());
        }
        oStr.append("</");
        oStr.append((String)oElementStack.pop());
        oStr.append(">");
        tLastWritePCData = false;
	tLastWriteStartElement = false;
        return oStr.substring(0);
    }

    private void increaseIndent()
    {
        oIndent.append("  ");
    }

    public String indent()
    {
        return oIndent.substring(0);
    }

    public String nl()
    {
        return oLineSeparator;
    }

    public String startElement(String string)
    {
        oElementStack.push(string);
        oStr.setLength(0);
        if (!tLastWritePCData)
            oStr.append(nl());
        oStr.append(indent());
        oStr.append("<");
        oStr.append(string);
        oStr.append(">");
        increaseIndent();
        tLastWritePCData = false;
	tLastWriteStartElement = true;
        return oStr.substring(0);
    }

    public String startElement(String string, Attributes attributes)
    {
        oElementStack.push(string);
        oStr.setLength(0);
        if (!tLastWritePCData)
            oStr.append(nl());
        oStr.append(indent());
        oStr.append("<");
        oStr.append(string);
        int i = 0;
        oSpaceStore.setLength(0);
        oSpaceStore.append(indent());
        for (i = 0; i <= string.length(); i++)
            oSpaceStore.append(" ");
        for (i = 0; i < attributes.getLength() - 1; i++)
        {
            //oAttName = attributes.getLocalName(i);
            oAttName = attributes.getQName(i);
            oAttValue = attributes.getValue(i);
            oStr.append(" ");
            if (i > 0)
                oStr.append(oSpaceStore.substring(0));
            oStr.append(oAttName);
            oStr.append("=\"");
            oStr.append(oAttValue);
            oStr.append("\"");
            oStr.append(nl());
        }

	//oAttName = attributes.getLocalName(i);
	oAttName = attributes.getQName(i);
        oAttValue = attributes.getValue(i);
        oStr.append(" ");
        if (attributes.getLength() > 1)
            oStr.append(oSpaceStore.substring(0));
        oStr.append(oAttName);
        oStr.append("=\"");
        oStr.append(oAttValue);
        oStr.append("\"");
        oStr.append(">");
        increaseIndent();
        tLastWritePCData = false;
	tLastWriteStartElement = true;
        return oStr.substring(0);
    }

    public String writeEmptyElement(String string)
    {
        oStr.setLength(0);
        if (!tLastWritePCData)
            oStr.append(nl());
        oStr.append(indent());
        oStr.append("<");
        oStr.append(string);
        oStr.append("/>");
        tLastWritePCData = false;
	tLastWriteStartElement = false;
        return oStr.substring(0);
    }

    public String writeEmptyElement(String string, Attributes attributes)
    {
        oStr.setLength(0);
        if (!tLastWritePCData)
            oStr.append(nl());
        oStr.append(indent());
        oStr.append("<");
        oStr.append(string);
        for (int i = 0; i < attributes.getLength(); i++)
        {
            //oAttName = attributes.getLocalName(i);
            oAttName = attributes.getQName(i);
            oAttValue = attributes.getValue(i);
            oStr.append(" ");
            oStr.append(oAttName);
            oStr.append("=\"");
            oStr.append(oAttValue);
            oStr.append("\"");
        }
        oStr.append("/>");
        tLastWritePCData = false;
	tLastWriteStartElement = false;
        return oStr.substring(0);
    }

    public String writeHeader()
    {
        oStr.setLength(0);
        oStr.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
        tLastWritePCData = false;
        return oStr.substring(0);
    }

    public String writePCData(String poPCData)
    {
        oStr.setLength(0);
        oStr.append(poPCData);
        tLastWritePCData = true;
	tLastWriteStartElement = false;
        return oStr.substring(0);
    }
}
