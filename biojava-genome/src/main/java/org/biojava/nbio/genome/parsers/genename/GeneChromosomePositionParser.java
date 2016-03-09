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
 * Author: ap3
 */

package org.biojava.nbio.genome.parsers.genename;

import org.biojava.nbio.genome.App;
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

/** A parser that parses a file from the UCSC genome browser that contains mapping of gene name to chromosome positions
 *
 * @author Andreas Prlic
 *
 */
public class GeneChromosomePositionParser {

	private static final Logger logger = LoggerFactory.getLogger(App.class);

	public static final String DEFAULT_MAPPING_URL="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz";

	public static void main(String[] args){
		try {

			List<GeneChromosomePosition> genePositions=	getChromosomeMappings();
			logger.info("got {} gene positions", genePositions.size());

			for (GeneChromosomePosition pos : genePositions){
				if ( pos.getGeneName().equals("FOLH1")) {
					logger.info("Gene Position: {}", pos);
					break;
				}
			}

		} catch(Exception e){
			logger.error("Exception: ", e);
		}
	}

	public static List<GeneChromosomePosition> getChromosomeMappings() throws IOException {

		URL url = new URL(DEFAULT_MAPPING_URL);

		InputStreamProvider prov = new InputStreamProvider();

		InputStream inStream = prov.getInputStream(url);

		return getChromosomeMappings(inStream);
	}

	public static List<GeneChromosomePosition> getChromosomeMappings(InputStream inStream) throws IOException {
		BufferedReader reader = new BufferedReader(new InputStreamReader(inStream));

		ArrayList<GeneChromosomePosition> gcps = new ArrayList<GeneChromosomePosition>();

		String line = null;
		while ((line = reader.readLine()) != null) {
			GeneChromosomePosition gcp = getGeneChromosomePosition(line);
			if ( gcp != null)
				gcps.add(gcp);
		}

		// since this is a large list, remove empty content.
		gcps.trimToSize();
		return gcps;
	}

	private static GeneChromosomePosition getGeneChromosomePosition(String line) {
		if ( line == null)
			return null;
		String[] spl = line.split("\t");

		if ( spl.length != 11) {
			logger.warn("Line does not have 11 data items, but {}: {}", spl.length, line);
			return null;
		}

		GeneChromosomePosition g = new GeneChromosomePosition();

		g.setGeneName(spl[0]);
		g.setGenebankId(spl[1]);
		g.setChromosome(spl[2]);
		g.setOrientation(spl[3].charAt(0));
		g.setTranscriptionStart(Integer.parseInt(spl[4]));
		g.setTranscriptionEnd(Integer.parseInt(spl[5]));
		g.setCdsStart(Integer.parseInt(spl[6]));
		g.setCdsEnd(Integer.parseInt(spl[7]));
		g.setExonCount(Integer.parseInt(spl[8]));
		String exonStarts = spl[9];
		String exonEnds = spl[10];
		g.setExonStarts(getIntegerList(exonStarts));
		g.setExonEnds(getIntegerList(exonEnds));

		//System.out.println(line);
		//System.out.println(Arrays.asList(spl) + " " + spl.length);
		return g;
	}

	private static List<Integer> getIntegerList(String lst){
		String[] spl = lst.split(",");
		ArrayList<Integer> l = new ArrayList<Integer>();
		for (String s : spl){
			l.add(Integer.parseInt(s));
		}
		l.trimToSize();
		return l;
	}
}
