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

package org.biojavax.utils;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import java.util.regex.Pattern;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

/**
 * Utility class for reading chunks of XML files and feeding them to SAX.
 * @author Richard Holland
 * @since 1.5
 */
public class XMLTools {
    
    // Static methods so should never be instantiated.
    private XMLTools() {}
    
    /**
     * Attempts to read XML file in chunks, passing each chunk to a SAX parser.
     * As each chunk is read into memory in a buffer, you need to ensure that each chunk
     * is small enough to fit into available memory. Only one chunk is held in memory
     * at any one time, and then only long enough for it to be parsed.
     * When checking for the presence of further chunks, it'll only read up to 1000 chars
     * further into the file, after which results will be unpredictable.
     * @param reader the reader to read the XML from
     * @param m_handler the SAX parser to feed the XML to
     * @param chunkToken the token to read. The parser will locate the first instance of
     * &lt;chunkToken and will buffer all content, including the opening tag and up to
     * and including the closing &lt;/chunkToken> tag. It will not currently handle
     * &lt;chunkToken/> instances, nor instances where more than one tag appears per line,
     * or extra spaces appear between the angle brackets, slashes, and tag name of the
     * tag we are searching for.
     * @return true if there is another chunk left to read after this one, false if not.
     * @throws ParserConfigurationException if there was a problem setting up the SAX parser.
     * @throws SAXException if there was a problem parsing the XML.
     * @throws IOException if there was a problem reading the XML from the reader.
     */
    public static boolean readXMLChunk(BufferedReader reader, DefaultHandler m_handler, String chunkToken) throws ParserConfigurationException, SAXException, IOException {
        // read next chunk from <chunkToken> to <chunkToken/> inclusive into buffer
        // process buffer through XML parser
        StringBuffer buffer = new StringBuffer();

        Pattern start = Pattern.compile(".*<"+chunkToken+".*");
        Pattern end = Pattern.compile(".*</"+chunkToken+">.*");
        
        boolean begunChunk = false;
        boolean filledBuffer = false;
        String line = null;
        while (!filledBuffer && (line=reader.readLine())!=null) {
            line = line.trim();
            if (!begunChunk && !start.matcher(line).matches()) continue;
            else begunChunk = true;
            buffer.append(line+"\n");
            if (end.matcher(line).matches()) filledBuffer = true;
        }
        if (!filledBuffer) throw new SAXException("Unexpectedly reached end of file");
        reader.mark(10000);
        boolean hasAnotherChunk = false;
        while (!hasAnotherChunk && (line=reader.readLine())!=null) {
            line = line.trim();
            if (start.matcher(line).matches()) hasAnotherChunk = true;
        }
        reader.reset();
        
        SAXParser m_xmlParser;
        SAXParserFactory factory = SAXParserFactory.newInstance();
        factory.setValidating(true);
        m_xmlParser = factory.newSAXParser();
        
        InputSource source = new InputSource(new StringReader(buffer.toString()));
        m_xmlParser.parse(source, m_handler);
        
        // return true if there are more in our buffer
        return hasAnotherChunk;
    }
}
