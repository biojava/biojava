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
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import org.biojava.utils.ClassTools;
import org.xml.sax.EntityResolver;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;

/**
 * SAX EntityResolve which looks up system IDs as resources
 * from a Java ClassLoader.
 *
 * @author Thomas Down
 * @since 1.2
 */

public class ResourceEntityResolver implements EntityResolver {
    private String[] path;
    private ClassLoader classLoader;

    private String[] parsePath(String path) {
	List pathElements = new ArrayList();
	StringTokenizer toke = new StringTokenizer(path, ":");
	while (toke.hasMoreTokens()) {
	    pathElements.add(toke.nextToken());
	}
	
	return (String[]) pathElements.toArray(new String[0]);
    }

    /**
     * Construct a resolver which searches for resources in the specified
     * path relative to the current classloader.
     */

    public ResourceEntityResolver(String path) {
	super();
	this.path = parsePath(path);
	this.classLoader = ClassTools.getClassLoader(this);
    }
    
    /**
     * Construct a resolver which searches for resources in the specified
     * list of directories relative to the current classloader.
     */

    public ResourceEntityResolver(String[] path) {
	super();
	this.path = path;
	this.classLoader = ClassTools.getClassLoader(this);
    }

    /**
     * Construct a resolver which searches for resources in the specified
     * list of directories relative to the supplied classloader.
     */

    public ResourceEntityResolver(String[] path, ClassLoader classLoader) {
	this.path = path;
	this.classLoader = classLoader;
    }

    /**
     * Construct a resolver which searches for resources in the specified
     * path relative to the supplied classloader.
     */

    public ResourceEntityResolver(String path, ClassLoader classLoader) {
	this.path = parsePath(path);
	this.classLoader = classLoader;
    }

    public InputSource resolveEntity(String publicId,
				     String systemId)
	throws SAXException, IOException
    {
	int index = systemId.lastIndexOf('/');
	if (index >= 0) {
	    systemId = systemId.substring(index + 1);
	}

	for (int i = 0; i < path.length; ++i) {
	    InputStream is = classLoader.getResourceAsStream(path[i] + "/" + systemId);
	    if (is != null) {
		InputSource source = new InputSource(is);
		source.setPublicId(publicId);
		source.setSystemId(systemId);
		return source;
	    }
	}

	return null;
    }
}

