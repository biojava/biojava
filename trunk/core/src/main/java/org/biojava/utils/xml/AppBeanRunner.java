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

import java.io.FileInputStream;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.w3c.dom.Document;
import org.xml.sax.InputSource;

/**
 * Create a bean from an XML file, then attempt to enter it.
 *
 * @author Thomas Down
 */

public class AppBeanRunner {
    public static void main(String[] args) {
	try {
	    if (args.length != 1) {
		throw new RuntimeException("Usage: java eponine.AppBeanRunner app.xml");
	    }
	    String name = args[0];

	    String[] altArgs = new String[args.length - 1];
	    for (int i = 1; i < args.length; ++i)
		altArgs[i-1] = args[i];
	    
	    InputSource is = new InputSource(new FileInputStream(name));
	    DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
	    dbf.setValidating(false);
	    DocumentBuilder parser = dbf.newDocumentBuilder();
	    Document doc = parser.parse(is);
	    Object bean = XMLBeans.INSTANCE.instantiateBean(doc.getDocumentElement());
	    if (! (bean instanceof AppEntry)) {
		System.out.println("Application can't be entered");
		return;
	    }
	    
	    ((AppEntry) bean).start(altArgs);
	} catch (Exception ex) {
	    ex.printStackTrace();
	}
    }
}
