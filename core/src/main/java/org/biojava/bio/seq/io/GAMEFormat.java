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

package org.biojava.bio.seq.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;

import javax.xml.parsers.SAXParserFactory;

import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.io.game12.GAMEHandler;
import org.biojava.utils.stax.SAX2StAXAdaptor;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.XMLReader;

/**
 * A rudimentary read-only GAME 1.2 Format object.
 *
 * @author David Huen
 */
public class GAMEFormat implements SequenceFormat
{
    public static final String DEFAULT = "GAME1.2";

    /**
     * this version only reads annotations (no symbols)
     */
    public boolean readSequence(BufferedReader reader, SymbolTokenization symParser, SeqIOListener listener)
        throws IOException
    {
        try {
            // set up processing pipeline
            InputSource is = new InputSource(reader);

            GAMEHandler handler = new GAMEHandler(listener);

            XMLReader parser;
            try {
                SAXParserFactory spf = SAXParserFactory.newInstance();
                spf.setValidating(false);
                spf.setNamespaceAware(true);
                parser = spf.newSAXParser().getXMLReader();
            } catch (Exception ex) {
                throw new IOException("Error creating SAX parser");
            }
            parser.setContentHandler(new SAX2StAXAdaptor(handler));

            parser.parse(is);

            return false;
        }
        catch (SAXException se) {
            se.printStackTrace();
            throw new IOException("SAXException encountered during parsing");
        }
    }

    public void writeSequence(Sequence seq, PrintStream os)
    {

    }

    public void writeSequence(Sequence seq, String format, PrintStream os)
    {

    }

    /**
     * <code>getDefaultFormat</code> returns the String identifier for
     * the default format.
     *
     * @return a <code>String</code>.
     * @deprecated
     */
    public String getDefaultFormat()
    {
        return DEFAULT;
    }
}

