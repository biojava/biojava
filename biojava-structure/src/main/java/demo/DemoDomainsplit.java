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
package demo;


import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.domain.LocalProteinDomainParser;
import org.biojava.nbio.structure.domain.pdp.Domain;
import org.biojava.nbio.structure.domain.pdp.Segment;
import org.biojava.nbio.structure.io.FileParsingParameters;

import java.util.List;

public class DemoDomainsplit {

	public static void main(String[] args){

		DemoDomainsplit split = new DemoDomainsplit();

		//String pdbId = "3gly";
		String pdbId = "4hhb";

		split.basicLoad(pdbId);

	}

	public void basicLoad(String pdbId){

		try {

			// This utility class can automatically download missing PDB files.
			AtomCache cache = new AtomCache();

			//
			// configure the parameters of file parsing (optional)

			FileParsingParameters params = new FileParsingParameters();

			// should the ATOM and SEQRES residues be aligned when creating the internal data model?
			params.setAlignSeqRes(true);
			// should secondary structure get parsed from the file
			params.setParseSecStruc(false);

			// and set the params in the cache.
			cache.setFileParsingParams(params);

			// end of optional part

			Structure struc = cache.getStructure(pdbId);

			System.out.println("structure loaded: " + struc);

			List<Domain> domains = LocalProteinDomainParser.suggestDomains(struc);

			System.out.println("RESULTS: =====");
			for ( Domain dom : domains){
				System.out.println("DOMAIN:" + dom.getSize() + " " +  dom.getScore());
				List<Segment> segments = dom.getSegments();
				for ( Segment s : segments){
					System.out.println("   Segment: " + s);
				}
			}
		} catch (Exception e){
			e.printStackTrace();
		}

	}
}
