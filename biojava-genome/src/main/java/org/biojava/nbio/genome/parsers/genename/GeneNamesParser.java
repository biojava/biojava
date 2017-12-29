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
 * created at 28 Jan 2014
 * Author: Andreas Prlic
 */

package org.biojava.nbio.genome.parsers.genename;

import org.biojava.nbio.core.util.InputStreamProvider;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

/** 
 * Parses a file from the www.genenames.org website that contains a mapping of human gene names to other databases
 *
 * @author Andreas Prlic
 *
 */
public class GeneNamesParser {

	private static final Logger logger = LoggerFactory.getLogger(GeneNamesParser.class);

	public static final String DEFAULT_GENENAMES_URL = "https://www.genenames.org/cgi-bin/download?title=HGNC+output+data&hgnc_dbtag=on&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_prev_sym&col=gd_prev_name&col=gd_aliases&col=gd_pub_chrom_map&col=gd_pub_acc_ids&col=md_mim_id&col=gd_pub_refseq_ids&col=md_ensembl_id&col=md_prot_id&col=gd_hgnc_id" +
			 "&status=Approved&status_opt=2&where=((gd_pub_chrom_map%20not%20like%20%27%patch%%27%20and%20gd_pub_chrom_map%20not%20like%20%27%ALT_REF%%27)%20or%20gd_pub_chrom_map%20IS%20NULL)%20and%20gd_locus_group%20%3d%20%27protein-coding%20gene%27&order_by=gd_app_sym_sort&format=text&limit=&submit=submit&.cgifields=&.cgifields=chr&.cgifields=status&.cgifields=hgnc_dbtag";

	/** parses a file from the genenames website
	 *
	 * @param args
	 */
	public static void main(String[] args) {

		try {

			List<GeneName> geneNames = getGeneNames();

			logger.info("got {} gene names", geneNames.size());

			for ( GeneName g : geneNames){
				if ( g.getApprovedSymbol().equals("FOLH1"))
					logger.info("Gene Name: {}", g);
			}
			// and returns a list of beans that contains key-value pairs for each gene name

		} catch (Exception e) {
			// TODO Auto-generated catch block
			logger.error("Exception: ", e);
		}

	}


	public static List<GeneName> getGeneNames() throws IOException{
		URL url = new URL(DEFAULT_GENENAMES_URL);

		InputStreamProvider prov = new InputStreamProvider();

		InputStream inStream = prov.getInputStream(url);

		return getGeneNames(inStream);
	}

	/** Get a list of GeneNames from an input stream.
	 *
	 * @param inStream
	 * @return list of geneNames
	 * @throws IOException
	 */
	public static List<GeneName> getGeneNames(InputStream inStream) throws IOException{

		ArrayList<GeneName> geneNames = new ArrayList<GeneName>();
		BufferedReader reader = new BufferedReader(new InputStreamReader(inStream));

		// skip reading first line (it is the legend)
		String line = reader.readLine();

		while ((line = reader.readLine()) != null) {
			// process line...
			//System.out.println(Arrays.toString(line.split("\t")));

			GeneName  geneName = getGeneName(line);
			if ( geneName != null)
				geneNames.add(geneName);
				//System.out.println(geneName);

		}

		// since this is a large list, let's free up unused space...
		geneNames.trimToSize();
		return geneNames;
	}

	private static GeneName getGeneName(String line) {
		// data is in this order:
		//[HGNC ID, Approved Symbol, Approved Name, Status, Previous Symbols,
		// Previous Names, Synonyms, Chromosome, Accession Numbers, RefSeq IDs, UniProt ID(supplied by UniProt)]

		if (line == null)
			return null;

		String[] s = line.split("\t");

		if ( s.length != 13) {
			logger.warn("Line does not contain 13 data items, but {}: {}", s.length, line);
			logger.warn(line.replaceAll("\t", "|---|"));
			return null;
		}
		GeneName gn = new GeneName();


		gn.setApprovedSymbol(s[0]);
		gn.setApprovedName(s[1]);
		gn.setStatus(s[2]);
		gn.setPreviousSymbols(s[3]);
		gn.setPreviousNames(s[4]);
		gn.setSynonyms(s[5]);
		gn.setChromosome(s[6]);
		gn.setAccessionNr(s[7]);
		gn.setOmimId(s[8]);
		gn.setRefseqIds(s[9]);
		gn.setEnsemblGeneId(s[10]);
		gn.setUniprot(s[11]);
		gn.setHgncId(s[12]);

		return gn;

	}

}
