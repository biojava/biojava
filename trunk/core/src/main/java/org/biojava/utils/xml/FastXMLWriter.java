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
package org.biojava.utils.xml;

import java.io.IOException;
import java.io.PrintWriter;

/**
 * Simple implementation of XMLWriter, optimized for speed.  The output is
 * not necessarily human-readable, but is fine for automated parsing.
 *
 * @author Thomas Down
 */

public class FastXMLWriter /* implements XMLWriter */ {
    private PrintWriter writer;
    private boolean isOpeningTag = false;

    public FastXMLWriter(PrintWriter writer) {
	this.writer = writer;
    }

    public void openTag(String qName)
        throws IOException
    {
	if (isOpeningTag) {
	    writer.print('>');
	}
	writer.print('<');
	writer.print(qName);
	isOpeningTag = true;
    }

    public void attribute(String qName, String value)
        throws IOException
    {
	if (! isOpeningTag) {
	    throw new IOException("attributes must follow an openTag");
	}

	writer.print(' ');
	writer.print(qName);
	writer.print("=\"");
	printAttributeValue(value);
	writer.print('"');
    }

    public void closeTag(String qName) 
        throws IOException
    {
	if (isOpeningTag) {
	    writer.println(" />");
	} else {
	    writer.print("</");
	    writer.print(qName);
	    writer.print('>');
	}

	isOpeningTag = false;
    }

    public void println(String data)
        throws IOException
    {
	if (isOpeningTag) {
	    writer.println('>');
	    isOpeningTag = false;
	}
	printChars(data);
	writer.println();
    }

    public void print(String data)
        throws IOException
    {
	if (isOpeningTag) {
	    writer.print('>');
	    isOpeningTag = false;
	}
	printChars(data);
    }


    public void printRaw(String data)
        throws IOException
    {
	writer.println(data);
    }
    
    protected void printChars(String data) 
        throws IOException
    {
	if (data == null) {
	    printChars("null");
	    return;
	}

	for (int pos = 0; pos < data.length(); ++pos) {
	    char c = data.charAt(pos);
	    if (c == '<' || c == '>' || c == '&') {
		numericalEntity(c);
	    } else {
		writer.write(c);
	    }
	}
    }

    protected void printAttributeValue(String data) 
        throws IOException
    {
	if (data == null) {
	    printAttributeValue("null");
	    return;
	}

	for (int pos = 0; pos < data.length(); ++pos) {
	    char c = data.charAt(pos);
	    if (c == '<' || c == '>' || c == '&' || c == '"') {
		numericalEntity(c);
	    } else {
		writer.write(c);
	    }
	}
    }

    protected void numericalEntity(char c)
        throws IOException
    {
	writer.print("&#");
	writer.print((int) c);
	writer.print(';');
    }
}
