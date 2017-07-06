/*
 *                  BioJava development code
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
 * Created on Dec 21, 2005
 *
 */
package org.biojava.nbio.structure;

import org.biojava.nbio.structure.io.PDBFileParser;

import java.io.InputStream;
import java.util.HashMap;
import java.util.Map;
import java.util.zip.GZIPInputStream;


/** 
 * A class that provides a set of standard amino acids.
 *
 *
 * @author Andreas Prlic
 * @author Tamas Horvath provided the standard amino acids.
 *
 *
 */
public final class StandardAminoAcid {

	private static final String STANDARD_AMINOS_FILE = "org/biojava/nbio/structure/standardaminos.pdb.gz";

	static private Map<String,AminoAcid> aminoAcids;

	/**
	 * Cannot be instantiated.
	 */
	private StandardAminoAcid() {
	}

	/**
	 * <p>
	 * Initialize the static StandardAminoAcid resource.
	 * </p>
	 *
	 * <p>
	 * This parses the resource
	 * <code>{@value #STANDARD_AMINOS_FILE}</code>
	 * and builds a basic set of amino acids.
	 *</p>
	 * @author Tamas Horvath provided the standard amino acids
	 */
	static {
		aminoAcids = new HashMap<String,AminoAcid>();


		InputStream fileStream = StandardAminoAcid.class.getClassLoader().getResourceAsStream(STANDARD_AMINOS_FILE);
		if (fileStream == null) {
			throw new RuntimeException("Could not find resource "+STANDARD_AMINOS_FILE+".  This probably means that your biojava.jar file is corrupt or incorrectly built.");
		}



		try {
			GZIPInputStream gzipIS = new GZIPInputStream(fileStream);
			PDBFileParser parser = new PDBFileParser();
			Structure s = parser.parsePDBFile(gzipIS);


			GroupIterator iter = new GroupIterator(s);
			while (iter.hasNext()){
				Group g = iter.next();

				if ( g instanceof AminoAcid){
					AminoAcid aa = (AminoAcid)g;

					aminoAcids.put(aa.getPDBName(),aa);
					aminoAcids.put(aa.getAminoType().toString(),aa);

				}
			}

		} catch (Exception t) {
			throw new RuntimeException( "Unable to initialize standard aminoacids", t);
		}
	}

	/** get a standard amino acid.
	 *
	 * @param name the 3- or 1-letter representation of the amino acid.
	 * @return the amino acids, or null if the name can not be matched
	 */
	public static AminoAcid getAminoAcid(String name){

		return aminoAcids.get(name);
	}

}
