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
package org.biojava.nbio.aaproperties.xml;

import javax.xml.bind.SchemaOutputResolver;
import javax.xml.transform.Result;
import javax.xml.transform.stream.StreamResult;
import java.io.File;
import java.io.IOException;

public class SchemaGenerator extends SchemaOutputResolver{
	private String fileName;

	public SchemaGenerator(String filename){
		this.fileName = filename;
	}

	@Override
	public Result createOutput(String namespaceUri, String suggestedFileName) throws IOException {
		File f = new File(this.fileName);
		f.createNewFile();
		return new StreamResult(f);
	}
}


