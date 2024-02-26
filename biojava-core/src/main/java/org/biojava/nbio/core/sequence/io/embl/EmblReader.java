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
package org.biojava.nbio.core.sequence.io.embl;


import java.io.*;
import java.util.Arrays;
import java.util.LinkedList;


/**
 * This class should process the data of embl file
 *
 * @author Noor Aldeen Al Mbaidin
 * @since 5.0.0
 */
public class EmblReader {

	/**
	 * The parsing is done in this method.<br>
	 * This method tries to process all the Embl records
	 * in the File , closes the underlying resource,
	 * and return the results in object of EmblRecord.<br>
	 *
	 * @return EmblRecord containing all the parsed Embl records
	 * @throws IOException
	 */
	public static EmblRecord process(File file) throws IOException {

		EmblRecord emblRecord = new EmblRecord();
		StringBuilder sequence = new StringBuilder("");
		LinkedList<EmblReference> emblReferences = new LinkedList<>();
		EmblReference emblReference = new EmblReference();
		LinkedList<String> accessionNumber = new LinkedList<>();
		LinkedList<String> keyword = new LinkedList<>();

		if (file == null)
			throw new NullPointerException("file can't be null");

		if (file.isDirectory())
			throw new IllegalArgumentException("the file can't be a directory");

		try (FileReader fileReader = new FileReader(file)) {
			String line = "";
			String lineIdentifier;
			String lineInfo;
			try (BufferedReader bufferedReader = new BufferedReader(fileReader)) {
				while ((line = bufferedReader.readLine()) != null) {
					if (line.length() > 1) {
						lineInfo = line.substring(2, line.length()).trim();
						lineIdentifier = line.substring(0, 2);
						if ("ID".equals(lineIdentifier))
							emblRecord.setEmblId(populateID(lineInfo));
						else if ("AC".equals(lineIdentifier))
							populateAccessionNumber(line, accessionNumber);
						else if ("DT".equals(lineIdentifier) && line.contains("Created"))
							emblRecord.setCreatedDate(lineInfo);
						else if ("DT".equals(lineIdentifier) && line.contains("updated"))
							emblRecord.setLastUpdatedDate(lineInfo);
						else if ("DE".equals(lineIdentifier))
							emblRecord.setSequenceDescription(lineInfo);
						else if ("KW".equals(lineIdentifier))
							keyword.add(lineInfo);
						else if ("OS".equals(lineIdentifier))
							emblRecord.setOrganismSpecies(lineInfo);
						else if ("OC".equals(lineIdentifier))
							emblRecord.setOrganismClassification(lineInfo);
						else if ("OG".equals(lineIdentifier))
							emblRecord.setOrGanelle(lineInfo);
						else if ("RN".equals(lineIdentifier) || "RP".equals(lineIdentifier)
								|| "RX".equals(lineIdentifier) || "RG".equals(lineIdentifier)
								|| "RA".equals(lineIdentifier) || "RT".equals(lineIdentifier)
								|| "RL".equals(lineIdentifier))
							populateEmblReferences(lineIdentifier, lineInfo, emblReference, emblReferences);
						else if ("DR".equals(lineIdentifier))
							emblRecord.setDatabaseCrossReference(lineInfo);
						else if ("AH".equals(lineIdentifier))
							emblRecord.setAssemblyHeader(lineInfo);
						else if ("AS".equals(lineIdentifier))
							emblRecord.setAssemblyInformation(lineInfo);
						else if ("CO".equals(lineIdentifier))
							emblRecord.setConstructedSequence(lineInfo);
						else if ("FH".equals(lineIdentifier))
							emblRecord.setFeatureHeader(lineInfo);
						else if ("FT".equals(lineIdentifier))
							emblRecord.setFeatureTable(lineInfo);
						else if ("SQ".equals(lineIdentifier))
							emblRecord.setSequenceHeader(lineInfo);
						else if ("  ".equals(lineIdentifier) && !"//".equals(lineIdentifier))
							populateSequence(line, sequence);
						else if ("//".equals(lineIdentifier)) {
							emblRecord.setKeyword(keyword);
							emblRecord.setEmblReference(emblReferences);
							emblRecord.setAccessionNumber(accessionNumber);
							emblRecord.setSequence(sequence.toString());
						}

					}
				}
			}
		}

		return emblRecord;
	}

	private static void populateSequence(String line, StringBuilder sequence) {
		String sequenceLine = line.replace(" ", "").
				replaceAll("[0-9]", "");
		sequence.append(sequenceLine);
	}

	private static void populateEmblReferences(String lineIdentifier, String lineInfo, EmblReference emblReference
			, LinkedList<EmblReference> emblReferences) {
		if ("RN".equals(lineIdentifier))
			emblReference.setReferenceNumber(lineInfo);
		else if ("RP".equals(lineIdentifier))
			emblReference.setReferencePosition(lineInfo);
		else if ("RX".equals(lineIdentifier))
			emblReference.setReferenceCrossReference(lineInfo);
		else if ("RG".equals(lineIdentifier))
			emblReference.setReferenceGroup(lineInfo);
		else if ("RA".equals(lineIdentifier))
			emblReference.setReferenceAuthor(lineInfo);
		else if ("RT".equals(lineIdentifier))
			emblReference.setReferenceTitle(lineInfo);
		else if ("RL".equals(lineIdentifier)) {
			emblReference.setReferenceLocation(lineInfo);
			emblReferences.add(emblReference.copyEmblReference(emblReference));
		}
	}

	private static void populateAccessionNumber(String line, LinkedList<String> accessionNumber) {
		accessionNumber.add(line);
	}

	private static EmblId populateID(String line) {
		String[] strings = line.split(";");
		Arrays.stream(strings).map(String::trim).toArray(unused -> strings);
		EmblId emblId = new EmblId(strings[0], strings[1], strings[2]
				, strings[3], strings[4], strings[5], strings[6]);
		return emblId;
	}


}
