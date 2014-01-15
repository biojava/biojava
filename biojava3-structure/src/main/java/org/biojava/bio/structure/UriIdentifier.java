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
 * Created on December 19, 2013
 * Author: Douglas Myers-Turnbull
 */

package org.biojava.bio.structure;

import java.io.File;
import java.io.FileNotFoundException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.List;

/**
 * A {@link StructureIdentifier} that uses URIs to identify structures.
 * TODO Finish
 * @author dmyersturnbull
 */
public class UriIdentifier implements StructureIdentifier {

	private final URI uri;
	
	public UriIdentifier(String uri) throws URISyntaxException {
		this.uri = new URI(uri);
	}
	
	public UriIdentifier(File file) throws FileNotFoundException {
		uri = file.toURI();
	}
	
	public UriIdentifier(URI uri) {
		this.uri = uri;
	}

	public URI getUri() {
		return uri;
	}

	@Override
	public String getIdentifier() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getPdbId() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<ResidueRange> getResidueRanges() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<String> getRanges() {
		// TODO Auto-generated method stub
		return null;
	}

}
