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
 * Created on August 13, 2010
 * Author: Mark Chapman
 */

package org.biojava3.alignment.io;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.RNASequence;
import org.biojava3.core.sequence.template.AbstractCompound;
import org.biojava3.core.sequence.template.AbstractSequence;

/**
 * Stores all the content of a Stockholm file.
 * <i><b>N.B.: This structure will undergo several enhancements later on. Don't depend on it in a final code, otherwise it will be hard to maintain.</b></i>
 * 
 * In general, Stockholm File contains the alignment mark-up lines.<br>
 * <br>
 * 
 * <Table border="1" align="center">
 * <tr><td><b>Header Section</b></td></tr>
 * <tr><td><b>Reference Section</b></td></tr>
 * <tr><td><b>Comment Section</b></td></tr>
 * <tr><td><B>Alignment Section</B></td></tr>
 * </table>
 * 
 * Sequence letters may include any characters except whitespace. Gaps may be indicated by "." or "-".<br>
 * Mark-up lines may include any characters except whitespace. Use underscore ("_") instead of space.<br>
 * 
 * <Table border="1">
 * <th>section field</th><th>preferred location</th>
 * <tr><td>#=GF &lt;feature&gt; &lt;Generic per-File annotation, free text&gt;</td><td>Above the alignment</td>
 * <tr><td>#=GC &lt;feature&gt; &lt;Generic per-Column annotation, exactly 1 char per column&gt;</td><td>Below the alignment</td>
 * <tr><td>#=GS &lt;seqname&gt; &lt;feature&gt; &lt;Generic per-Sequence annotation, free text&gt;</td><td>Above the alignment or just below the corresponding sequence</td>
 * <tr><td>#=GR &lt;seqname&gt; &lt;feature&gt; &lt;Generic per-Residue annotation, exactly 1 char per residue&gt;</td><td>Just below the corresponding sequence</td>
 * </tr>
 * </table>
 * 
 * @since 3.0.5
 * @author Amr AL-Hossary
 * @author Marko Vaz
 * 
 */
public class StockholmStructure {

	

	private final StockholmFileAnnotation fileAnnotation;
	private final StockholmConsensusAnnotation consAnnotation;
	private final Map<String, StringBuffer> sequences;
	private final Map<String, StockholmSequenceAnnotation> seqsAnnotation;
	private final Map<String, StockholmResidueAnnotation> resAnnotation;

	public StockholmStructure() {
		fileAnnotation = new StockholmFileAnnotation();
		consAnnotation = new StockholmConsensusAnnotation();
		sequences = new HashMap<String, StringBuffer>();
		seqsAnnotation = new HashMap<String, StockholmSequenceAnnotation>();
		resAnnotation = new HashMap<String, StockholmResidueAnnotation>();
	}

	public StockholmFileAnnotation getFileAnnotation() {
		return fileAnnotation;
	}
	
	public StockholmConsensusAnnotation getConsAnnotation() {
		return consAnnotation;
	}

	/**Actually this function should be called appendToSequence
	 * @param seqName
	 * @param seqText
	 */
	public void addSequence(String seqName, String seqText) {
		StringBuffer seq = sequences.get(seqName);
		if (seq != null) {
			//add sequence without space
			seq.append(seqText);
		} else {
			seq = new StringBuffer(seqText);
			sequences.put(seqName, seq);
		}
	}

	public Map<String, StringBuffer> getSequences() {
		return sequences;
	}

	private StockholmSequenceAnnotation getSequenceAnnotation(String seqName) {
		if (!seqsAnnotation.containsKey(seqName)) {
			seqsAnnotation.put(seqName, new StockholmSequenceAnnotation());
		}
		return seqsAnnotation.get(seqName);
	}

	/** 
	 * @param seqName
	 * @param text
	 */
	public void addGSAccessionNumber(String seqName, String text) {
		getSequenceAnnotation(seqName).setAccessionNumber(text);
	}

	public void addGSDescription(String seqName, String text) {
		getSequenceAnnotation(seqName).addToDescription(text);
	}

	/**
	 * @param seqName
	 * @param text
	 */
	public void addGSdbReference(String seqName, String text) {
		getSequenceAnnotation(seqName).addDBReference(text);
	}

	public void addGSOrganismSpecies(String seqName, String text) {
		getSequenceAnnotation(seqName).setOrganism(text);
	}

	public void addGSOrganismClassification(String seqName, String text) {
		getSequenceAnnotation(seqName).setOrganismClassification(text);
	}

	public void addGSLook(String seqName, String text) {
		getSequenceAnnotation(seqName).setLook(text);
	}
	
	private StockholmResidueAnnotation getResidueAnnotation(String seqName) {
		if (!resAnnotation.containsKey(seqName)) {
			resAnnotation.put(seqName, new StockholmResidueAnnotation());
		}
		return resAnnotation.get(seqName);
	}

	public void addSurfaceAccessibility(String seqName, String text) {
		getResidueAnnotation(seqName).setSurfaceAccessibility(text);
	}

	public void addTransMembrane(String seqName, String text) {
		getResidueAnnotation(seqName).setTransMembrane(text);
	}

	public void addPosteriorProbability(String seqName, String text) {
		getResidueAnnotation(seqName).setPosteriorProbability(text);
	}

	public void addLigandBinding(String seqName, String text) {
		getResidueAnnotation(seqName).setLigandBinding(text);
	}

	public void addActiveSite(String seqName, String text) {
		getResidueAnnotation(seqName).setActiveSite(text);
	}

	public void addASPFamPredicted(String seqName, String text) {
		getResidueAnnotation(seqName).setAsPFamPredicted(text);
	}

	public void addASSwissProt(String seqName, String text) {
		getResidueAnnotation(seqName).setAsSwissProt(text);
	}

	public void addIntron(String seqName, String text) {
		getResidueAnnotation(seqName).setIntron(text);
	}

	public void addSecondaryStructure(String seqName, String text) {
		getResidueAnnotation(seqName).setSecondaryStructure(text);
	}


	/**
	 * 
	 * @return
	 */
	public List<AbstractSequence<? extends AbstractCompound>> getBioSequences() {
		return getBioSequences(false);
	}

	/**Because some database files have incorrectly small letters (e.g. Pfam23 structure 
	 * PF00389.22 sequence TKRA_BACSU/6-322), this 
	 * function is used to ignore the small letters case.
	 * 
	 * @param ignoreCase if <code>true</code>, the function will deal with small letters as if they are capital ones
	 * @return
	 */
	public List<AbstractSequence<? extends AbstractCompound>> getBioSequences(boolean ignoreCase) {
		List<AbstractSequence<? extends AbstractCompound>> seqs = new ArrayList<AbstractSequence<? extends AbstractCompound>>();
		for (String sequencename : sequences.keySet()) {
			AbstractSequence<? extends AbstractCompound> seq;
			String sequence = sequences.get(sequencename).toString();
			if (ignoreCase) {
				sequence=sequence.toUpperCase();
			}
			if (fileAnnotation.isPFam()) {
				seq = new ProteinSequence(sequence);
			} else {
				seq = new RNASequence(sequence);
			}

			String[] seqDetails = splitSeqName(sequencename);
			seq.setDescription(seqDetails[0]);
			seq.setBioBegin((seqDetails[1] == null || seqDetails[1].trim().equals("") ? null : new Integer(seqDetails[1])));
			seq.setBioEnd  ((seqDetails[2] == null || seqDetails[2].trim().equals("") ? null : new Integer(seqDetails[2])));

			seqs.add(seq);
		}
		return seqs;
	}

	/**
	 * Returns an array with the following sequence related content: name,
	 * start, end.
	 * 
	 * @param sequenceName
	 *            the sequence from where to extract the content. It is supposed
	 *            that it follows the following convention name/start-end (e.g.:
	 *            COATB_BPIKE/30-81)
	 * @return array with the following sequence related content: name, start,
	 *         end.
	 */
	private String[] splitSeqName(String sequenceName) {
		String[] result = new String[3];

		String[] barSplit = sequenceName.toString().split("/");
		if (barSplit.length == 2) {
			result[0] = barSplit[0];
			String[] positions = barSplit[1].split("-");
			if (positions.length == 2) {
				result[1] = positions[0];
				result[2] = positions[1];
			}
		} else {
			result[0] = sequenceName;
			result[1] = null;
			result[2] = null;
		}

		return result;
	}

	@Override
	public String toString() {
		StringBuffer result = new StringBuffer();
		List<AbstractSequence<? extends AbstractCompound>> bioSeqs = getBioSequences(false);
		int sequenceLength = -1;
		for (AbstractSequence<? extends AbstractCompound> sequence : bioSeqs) {
			String sequenceAsString = sequence.getSequenceAsString();
			sequenceLength = sequenceAsString.length();
			if (sequenceLength > 50) {
				result.append(sequenceAsString.substring(0, 40));
				result.append("...");
				result.append(sequenceAsString.substring(sequenceLength - 3, sequenceLength));
			} else {
				result.append(sequenceAsString);
			}
			result.append(" " + sequence.getDescription() + "\n");
		}
		result.append("Alignment with " + bioSeqs.size() + " rows and "+ sequenceLength + " columns");

		return result.toString();
	}
	
	public static class DatabaseReference{
		public static final String EXPERT="EXPERT"; 
		public static final String MIM="MIM"; 
		public static final String PFAMB="PFAMB"; 
		public static final String PRINTS="PRINTS"; 
		public static final String PROSITE="PROSITE"; 
		public static final String PROSITE_PROFILE="PROSITE_PROFILE"; 
		public static final String SCOP="SCOP"; 
		public static final String PDB="PDB"; 
		public static final String SMART="SMART"; 
		public static final String URL="URL"; 
		public static final String LOAD="LOAD"; 
		public static final String HOMSTRAD="HOMSTRAD"; 
		public static final String INTERPRO="INTERPRO"; 

		
		private String database;
		/**TODO this field should be subdivided into smaller fields if the database is SCOP or PDB.*/
		private String reference;
		
		public DatabaseReference(String database, String reference) {
			this.database=database;
			this.reference=reference;
		}
		
		public DatabaseReference(String representativeAnnotationString) {
			int semiColonIndex=representativeAnnotationString.indexOf(';');
			this.database=representativeAnnotationString.substring(0, semiColonIndex);
			this.reference=representativeAnnotationString.substring(semiColonIndex+1, representativeAnnotationString.lastIndexOf(';')).trim();
		}
		@Override
		public String toString() {
			return new StringBuilder(this.database).append(';').append(' ').
					append(this.reference).append(';').toString();
		}
		public String getDatabase() {
			return database;
		}
		
		public String getReference() {
			return reference;
		}
	}
}