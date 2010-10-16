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

import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

/**
 * A simple XML DocumentHandler that processes SAX2 events
 * to create a sensibly formatted XML as it parsed
 * without populating objects with data. Note, need to
 * tidy up parameter names which are poor due to decompiling
 * a .class file from an accidently deleted .java file!
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
 * @version 1.0
 *
 * @see BaseXMLWriter
 */
public class SimpleXMLEmitter extends DefaultHandler {
    private BaseXMLWriter oXMLWriter;
    private boolean tEmitQNames;
    private String oName;

    public SimpleXMLEmitter()
    {
        oXMLWriter = new BaseXMLWriter();
        tEmitQNames = true;
        System.out.println("<?xml version=\"1.0\"?>");
    }

    public SimpleXMLEmitter(boolean flag)
    {
        oXMLWriter = new BaseXMLWriter();
        tEmitQNames = true;
        setEmitQNames(flag);
        System.out.println("<?xml version=\"1.0\"?>");
    }

    public void characters(char ach[], int i, int j)
        throws SAXException
    {
        System.out.print(oXMLWriter.writePCData(new String(ach, i, j)));
    }

    public void endElement(String string1, String string2, String string3)
    {
        System.out.print(oXMLWriter.endElement());
    }

    private boolean isEmitQNames()
    {
        return tEmitQNames;
    }

    private void setEmitQNames(boolean flag)
    {
        tEmitQNames = flag;
    }

    public void startElement(String string1, String string2, String string3, Attributes attributes)
    {
        if (isEmitQNames())
            oName = string3;
        else
            oName = string2;
        if (attributes.getLength() != 0)
            System.out.print(oXMLWriter.startElement(oName, attributes));
        else
            System.out.print(oXMLWriter.startElement(oName));
    }
}
