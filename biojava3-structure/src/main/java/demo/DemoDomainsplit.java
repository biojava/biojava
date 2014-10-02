package demo;

import java.util.List;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.domain.LocalProteinDomainParser;
import org.biojava.bio.structure.domain.pdp.Domain;
import org.biojava.bio.structure.domain.pdp.Segment;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class DemoDomainsplit {

	private static final Logger logger = LoggerFactory.getLogger(DemoDomainsplit.class);

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
			params.setLoadChemCompInfo(true);
			// should secondary structure get parsed from the file
			params.setParseSecStruc(false);

			// and set the params in the cache.			
			cache.setFileParsingParams(params);

			// end of optional part
			
			Structure struc = cache.getStructure(pdbId);
			
			logger.info("structure loaded: {}", struc);
			
			List<Domain> domains = LocalProteinDomainParser.suggestDomains(struc);

			logger.info("RESULTS: =====");
			for ( Domain dom : domains){
				logger.info("DOMAIN:" + dom.getSize() + " " +  dom.getScore());
				List<Segment> segments = dom.getSegments();
				for ( Segment s : segments){
					logger.info("Segment: {}", s);
				}
			}
		} catch (Exception e){
			logger.error("Exception: ", e);
		}
	}
}