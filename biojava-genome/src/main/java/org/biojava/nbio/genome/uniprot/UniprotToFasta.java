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
package org.biojava.nbio.genome.uniprot;


import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.io.FastaWriterHelper;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author Scooter
 */
public class UniprotToFasta {

	private static final Logger logger = LoggerFactory.getLogger(UniprotToFasta.class);

	public static void main( String[] args ){
		try{
			String uniprotDatFileName = "uniprot_trembl_fungi.dat";
			String fastaFileName = "uniprot__trembel_fungi.faa";
			UniprotToFasta uniprotToFasta = new UniprotToFasta();
			uniprotToFasta.process(uniprotDatFileName, fastaFileName);
		}catch(Exception e){
			logger.error("Exception: ", e);
		}
	}

	/**
	 * Convert a Uniprot sequence file to a fasta file. Allows you to download all sequence data for a species
	 * and convert to fasta to be used in a blast database
	 * @param uniprotDatFileName
	 * @param fastaFileName
	 * @throws Exception
	 */

	public void process( String uniprotDatFileName,String fastaFileName ) throws Exception{

			FileReader fr = new FileReader(uniprotDatFileName);
			BufferedReader br = new BufferedReader(fr);
			String line = br.readLine();
			String id = "";
			StringBuffer sequence = new StringBuffer();
			ArrayList<ProteinSequence> seqCodingRegionsList = new ArrayList<ProteinSequence>();
			int count = 0;
			HashMap<String,String> uniqueGenes = new HashMap<String,String>();
			HashMap<String,String> uniqueSpecies = new HashMap<String,String>();
			while(line != null){
				if(line.startsWith("ID")){
					String[] data = line.split(" ");
					id = data[3];
				}else if(line.startsWith("SQ")){
					line = br.readLine();
					while(!line.startsWith("//")){

						for(int i = 0; i < line.length(); i++){
							char aa = line.charAt(i);
							if((aa >= 'A' && aa <= 'Z') || (aa >= 'a' && aa <= 'z' )){
								sequence.append(aa);
							}
						}
						line = br.readLine();
					}

				 //   System.out.println(">" + id);
				 //   System.out.println(sequence.toString());

					ProteinSequence seq = new ProteinSequence(sequence.toString() );
					seq.setAccession(new AccessionID(id));

					seqCodingRegionsList.add(seq);
					sequence = new StringBuffer();
					count++;
					if(count % 100 == 0)
						logger.info("Count: ", count);
					String[] parts = id.split("_");
					uniqueGenes.put(parts[0], "");
					uniqueSpecies.put(parts[1],"");
				}
				line = br.readLine();
			}
	   //     System.out.println("Unique Genes=" + uniqueGenes.size());
	   //     System.out.println("Unique Species=" + uniqueSpecies.size());
	   //     System.out.println("Total sequences=" + seqCodingRegionsList.size());
			FastaWriterHelper.writeProteinSequence(new File(fastaFileName), seqCodingRegionsList);

			br.close();
			fr.close();

	  //      System.out.println(uniqueGenes.keySet());
	  //      System.out.println("====================");
	  //      System.out.println(uniqueSpecies.keySet());


	}

}
