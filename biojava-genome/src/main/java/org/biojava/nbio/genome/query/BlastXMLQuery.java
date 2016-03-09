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
package org.biojava.nbio.genome.query;


import org.biojava.nbio.core.util.XMLHelper;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.util.ArrayList;
import java.util.LinkedHashMap;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class BlastXMLQuery {

	private static final Logger logger = LoggerFactory.getLogger(BlastXMLQuery.class);

	Document blastDoc = null;

	public BlastXMLQuery(String blastFile) throws Exception {
		logger.info("Start read of {}", blastFile);
		blastDoc = XMLHelper.loadXML(blastFile);
		logger.info("Read finished");
	}

	public LinkedHashMap<String, ArrayList<String>> getHitsQueryDef(double maxEScore) throws Exception {
		LinkedHashMap<String, ArrayList<String>> hitsHashMap = new LinkedHashMap<String, ArrayList<String>>();
		logger.info("Query for hits");
		ArrayList<Element> elementList = XMLHelper.selectElements(blastDoc.getDocumentElement(), "BlastOutput_iterations/Iteration[Iteration_hits]");
		logger.info("{} hits", elementList.size());

		for (Element element : elementList) {
			Element iterationquerydefElement = XMLHelper.selectSingleElement(element, "Iteration_query-def");
			String querydef = iterationquerydefElement.getTextContent();
			Element iterationHitsElement = XMLHelper.selectSingleElement(element, "Iteration_hits");
			ArrayList<Element> hitList = XMLHelper.selectElements(iterationHitsElement, "Hit");
			for (Element hitElement : hitList) {
				Element hitaccessionElement = XMLHelper.selectSingleElement(hitElement, "Hit_accession");
				String hitaccession = hitaccessionElement.getTextContent();
				Element hithspsElement = XMLHelper.selectSingleElement(hitElement, "Hit_hsps");
				ArrayList<Element> hspList = XMLHelper.selectElements(hithspsElement, "Hsp");
				for (Element hspElement : hspList) {
					Element evalueElement = XMLHelper.selectSingleElement(hspElement, "Hsp_evalue");
					String value = evalueElement.getTextContent();
					double evalue = Double.parseDouble(value);
					if (evalue <= maxEScore) {
						ArrayList<String> hits = hitsHashMap.get(querydef);
						if (hits == null) {
							hits = new ArrayList<String>();
							hitsHashMap.put(querydef, hits);
						}
						hits.add(hitaccession);
					}
				}
			}
		}

		return hitsHashMap;
	}

	public static void main(String[] args) {
		try {
			BlastXMLQuery blastXMLQuery = new BlastXMLQuery("/Users/Scooter/scripps/dyadic/analysis/454Scaffolds/c1-454Scaffolds-hits-uniprot_fungi.xml");
			LinkedHashMap<String, ArrayList<String>> hits = blastXMLQuery.getHitsQueryDef(1E-10);
			logger.info("Hits: {}", hits);
		} catch (Exception e) {
			logger.error("Execution: ", e);
		}
	}
}
