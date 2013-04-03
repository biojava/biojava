package demo;


import java.util.List;

import org.biojava.bio.structure.Structure;

import org.biojava.bio.structure.domain.LocalProteinDomainParser;
import org.biojava.bio.structure.domain.pdp.Domain;
import org.biojava.bio.structure.domain.pdp.Segment;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.PDBFileReader;

public class DemoDomainsplit {

	public static void main(String[] args){

		DemoDomainsplit split = new DemoDomainsplit();

		//String pdbId = "3gly";
		String pdbId = "4hhb";
		
		split.basicLoad(pdbId);

	}

	public void basicLoad(String pdbId){

		try {

			PDBFileReader reader = new PDBFileReader();

			// the path to the local PDB installation
			reader.setPath("/tmp/");

			// are all files in one directory, or are the files split,
			// as on the PDB ftp servers?
			reader.setPdbDirectorySplit(true);

			// should a missing PDB id be fetched automatically from the FTP servers?
			reader.setAutoFetch(true);

			// configure the parameters of file parsing

			FileParsingParameters params = new FileParsingParameters();

			// should the ATOM and SEQRES residues be aligned when creating the internal data model?
			params.setAlignSeqRes(true);
			params.setLoadChemCompInfo(true);
			// should secondary structure get parsed from the file
			params.setParseSecStruc(false);

			reader.setFileParsingParameters(params);

			Structure struc = reader.getStructureById(pdbId);
			
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
