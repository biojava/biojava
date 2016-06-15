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
 * created at 20 Feb 2014
 * Author: ap3
 */

package org.biojava.nbio.genome.parsers.cytoband;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.zip.GZIPInputStream;

/**
 * Parses the cytoband (karyotype) file from UCSC.
 *
 */
public class CytobandParser {

	private static final Logger logger = LoggerFactory
			.getLogger(CytobandParser.class);

	public static final String DEFAULT_LOCATION = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz";

	public static void main(String[] args) {

		CytobandParser me = new CytobandParser();
		try {
			SortedSet<Cytoband> cytobands = me.getAllCytobands(new URL(
					DEFAULT_LOCATION));
			SortedSet<StainType> types = new TreeSet<StainType>();
			for (Cytoband c : cytobands) {
				logger.info("Cytoband: {}", c);
				if (!types.contains(c.getType()))
					types.add(c.getType());
			}
			logger.info("Strain Type: {}", types);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			logger.error("Exception: ", e);
		}

	}

	public SortedSet<Cytoband> getAllCytobands(URL u) throws IOException {
		InputStream stream = new GZIPInputStream(u.openStream());
		return getAllCytobands(stream);

	}

	public SortedSet<Cytoband> getAllCytobands(InputStream instream)
			throws IOException {
		BufferedReader reader = new BufferedReader(new InputStreamReader(
				instream));
		String line = null;
		SortedSet<Cytoband> cytobands = new TreeSet<Cytoband>();
		while ((line = reader.readLine()) != null) {
			String[] spl = line.split("\t");
			if (spl.length != 5) {
				logger.warn(
						"WRONG LINE LENGHT, expected 5, but got {} for: {}",
						spl.length, line);
			}

			Cytoband b = new Cytoband();
			b.setChromosome(spl[0]);
			b.setStart(Integer.parseInt(spl[1]));
			b.setEnd(Integer.parseInt(spl[2]));
			b.setLocus(spl[3]);
			StainType type = StainType.getStainTypeFromString(spl[4]);
			if (type == null)
				logger.warn("unknown type: {}", spl[4]);
			b.setType(type);
			cytobands.add(b);
		}

		return cytobands;
	}

}
