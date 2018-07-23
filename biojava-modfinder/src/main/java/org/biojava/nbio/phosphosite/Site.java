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
package org.biojava.nbio.phosphosite;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;

/**
 * Created by ap3 on 31/10/2014.
 */
public class Site {

	private final static Logger logger = LoggerFactory.getLogger(Site.class);

	public Site(){


	}

	public static List<Site> parseSites(File f) throws IOException {

		InputStream inStream = new FileInputStream(f);
		InputStream gzipStream = new GZIPInputStream(inStream);

		Reader decoder = new InputStreamReader(gzipStream);
		BufferedReader buf = new BufferedReader(decoder);

		String line = null;

		List<Site > data = new ArrayList<Site>();

		List<String> headerFields = null;

		int proteinIndex = -1;
		int uniprotIndex = -1;
		int residueIndex = -1;
		int orgIndex     = -1;
		int groupIndex   = -1;
		int geneIndex    = -1;

		boolean inHeader = true;


		while ((line = buf.readLine()) != null){
			if ( line.startsWith("GENE") ||
					line.startsWith("PROTEIN")) {

				headerFields = parseHeaderFields(line);

				proteinIndex = headerFields.indexOf("PROTEIN");
				uniprotIndex = headerFields.indexOf("ACC_ID");
				residueIndex = headerFields.indexOf("MOD_RSD");
				orgIndex     = headerFields.indexOf("ORGANISM");
				groupIndex   = headerFields.indexOf("SITE_GRP_ID");
				geneIndex 	 = headerFields.indexOf("GENE");

				inHeader = false;
				continue;
			}
			if ( inHeader)
				continue;

			if ( line.trim().length() == 0)
				continue;

			// fields are:
			String[] spl = line.split("\t");
			if ( spl.length  < 5){
				logger.info("Found wrong line length: " + line);
				continue;

			}

			String protein = spl[proteinIndex];
			String uniprot = spl[uniprotIndex];

			String residue = spl[residueIndex];

			String[] resSpl = residue.split("-");
			String modType = null;
			if ( resSpl.length == 2) {

				 modType = resSpl[1];
			}
			String group    = spl[groupIndex];

			String organism = spl[orgIndex];

			String geneSymb = spl[geneIndex];

			Site s = new Site();
			s.setProtein(protein);
			s.setUniprot(uniprot);
			s.setGeneSymb(geneSymb);
			s.setModType(modType);
			s.setResidue(residue);
			s.setGroup(group);
			s.setOrganism(organism);
			data.add(s);

		}
		buf.close();

		return data;

	}

	private static List<String> parseHeaderFields(String line) {
		String[] spl = line.split("\t");

		List<String> h = new ArrayList<String>();
		for (String s: spl){
			h.add(s);

		}

		return h;
	}

	String protein;
	String uniprot;
	String geneSymb;
	String chrLoc;
	String modType;
	String residue ;
	String group;
	String organism;

	public String getProtein() {
		return protein;
	}

	public void setProtein(String protein) {
		this.protein = protein;
	}

	public String getUniprot() {
		return uniprot;
	}

	public void setUniprot(String uniprot) {
		this.uniprot = uniprot;
	}

	public String getGeneSymb() {
		return geneSymb;
	}

	public void setGeneSymb(String geneSymb) {
		this.geneSymb = geneSymb;
	}

	public String getChrLoc() {
		return chrLoc;
	}

	public void setChrLoc(String chrLoc) {
		this.chrLoc = chrLoc;
	}

	public String getModType() {
		return modType;
	}

	public void setModType(String modType) {
		this.modType = modType;
	}

	public String getResidue() {
		return residue;
	}

	public void setResidue(String residue) {
		this.residue = residue;
	}

	public String getGroup() {
		return group;
	}

	public void setGroup(String group) {
		this.group = group;
	}

	public String getOrganism() {
		return organism;
	}

	public void setOrganism(String organism) {
		this.organism = organism;
	}

	@Override
	public String toString() {
		StringBuffer s = new StringBuffer();

		s.append("Site{" +
				"protein='" + protein + '\'');
		if ( uniprot != null)
				s.append(", uniprot='" + uniprot + '\'' );
		if ( geneSymb != null)
			s.append(
				", geneSymb='" + geneSymb + '\'' );
		if (chrLoc != null)
				s.append(", chrLoc='" + chrLoc + '\'' );
		if (modType != null)
			s.append(", modType='" + modType + '\'' );

		if (residue != null)
			s.append(        ", residue='" + residue + '\'' );
		if ( group != null)
				s.append(", group='" + group + '\'' );
		if (organism != null)
			s.append(", organism='" + organism + '\'' );

		  s.append(      '}');

		return s.toString();
	}
}


