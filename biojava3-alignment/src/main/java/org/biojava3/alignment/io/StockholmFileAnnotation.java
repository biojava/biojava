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
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import org.biojava3.alignment.io.StockholmStructure.DatabaseReference;

/**
 * Stores all the content parsed from the #=GF lines
 * 
 * 
 * @since 3.0.5
 * @author Amr AL-Hossary
 * @author Marko Vaz
 * 
 */
public class StockholmFileAnnotation {
	public static class StockholmFileAnnotationReference {
		private String refMedline;//TODO 8 digits
		private CharSequence refTitle;//on several lines
		private CharSequence refAuthor;//TODO comma-separated, semicolon terminated list;
		/**<b>TODO to be formatted later on.</b><br>
		 * RL  Journal abbreviation year;volume:page-page.<br>
		 * RL   Virus Genes 1997;14:163-165.<br>
		 * RL   J Mol Biol 1994;242:309-320.
		 */
		private String refLocation;
		
		
		public String getRefMedline() {
			return refMedline;
		}
		public void setRefMedline(String refMedline) {
			this.refMedline = refMedline;
		}
		public String getRefTitle() {
			return refTitle.toString();
		}
		public void setRefTitle(String refTitle) {
			this.refTitle = refTitle;
		}
		public void addToRefTitle(String refTitle) {
			if (this.refTitle== null) {
				this.refTitle = new StringBuffer(refTitle);
			} else if (this.refTitle instanceof StringBuffer){
				((StringBuffer) this.refTitle).append(' ').append(refTitle);
			}else {
				this.refTitle = new StringBuffer(this.refTitle).append(' ').append(refTitle);
			}
		}
		
		public String getRefAuthor() {
			return refAuthor.toString();
		}
		public void setRefAuthor(StringBuffer refAuthor) {
			this.refAuthor = refAuthor;
		}
		public void addToRefAuthor(String refAuthor) {
			if (this.refAuthor== null) {
				this.refAuthor = new StringBuffer(refAuthor);
			} else if (this.refAuthor instanceof StringBuffer){
				((StringBuffer) this.refAuthor).append(' ').append(refAuthor);
			}else {
				this.refAuthor = new StringBuffer(this.refAuthor).append(' ').append(refAuthor);
			}
		}

		public String getRefLocation() {
			return refLocation;
		}
		public void setRefLocation(String refLocation) {
			this.refLocation = refLocation;
		}
		
		

	}

	//TODO revise these 4 fields usage
	private final static String TREE_DEFAULT_ID = "DEFAULT_ID";
	private static final String PB_PFAM_STRING = "PB";
	private static final String PF_PFAM_STRING = "PF";
	private static final String RF_RFAM_STRING = "RF";

	private StringBuffer format;
	private StringBuffer version;
	private String accessionNumber;
	private StringBuffer identification;
	private StringBuffer definition;
	private String[] authors;
	private String alignmentMethod;
	private CharSequence buildMethod;
	private StringBuffer searchMethod;
	private StringBuffer sourceSeed;
	private StringBuffer sourceStructure;
	private float[] gatheringThreshs;
	private float[] noiseCutoffs;
	private float[] trustedCutoffs;
	private String typeField;
	private String[] previousIDs;
	private int numSequences;
	private StringBuffer dbComment;
	private Set<DatabaseReference> dbReferences;
	private StringBuffer refComment;
	/**TODO When implementing toString(), the function should loop on the vector */
	private Vector<StockholmFileAnnotationReference> references = new Vector<StockholmFileAnnotation.StockholmFileAnnotationReference>();
	private StringBuffer keywords;
	private CharSequence comment;
	private StringBuffer pfamAccession;
	private StringBuffer location;
	private StringBuffer wikipediaLink;
	private StringBuffer clan;
	private StringBuffer membership;
	private final Map<String, List<String>> embTrees;
	private float falseDiscoveryRate;

	public StockholmFileAnnotation() {
		embTrees = new HashMap<String, List<String>>();
	}

	public StringBuffer getDbComment() {
		return dbComment;
	}

	public void setDbComment(String dbComment) {
		if (this.dbComment != null) {
			this.dbComment.append(dbComment);
		} else {
			this.dbComment = new StringBuffer(dbComment);
		}
	}

	public Set<DatabaseReference> getDbReferences() {
		return dbReferences;
	}

	public void setDbReferences(Set<DatabaseReference> dbReferences) {
		this.dbReferences = dbReferences;
	}
	/**
	 * @param dbReference the string without the initial annotation identifier ( #=GS DR )
	 */
	public void addDBReference(String dbReferenceRepresentingString) {
		if (this.dbReferences == null) {
			this.dbReferences = new HashSet<DatabaseReference>();
		} 
		dbReferences.add(new DatabaseReference(dbReferenceRepresentingString));
	}

	public float getFalseDiscoveryRate() {
		return falseDiscoveryRate;
	}

	public void setFalseDiscoveryRate(float falseDiscoveryRate) {
		this.falseDiscoveryRate = falseDiscoveryRate;
	}

	public StringBuffer getRefComment() {
		return refComment;
	}

	public StringBuffer getKeywords() {
		return keywords;
	}

	public String getComment() {
		return comment.toString();
	}

	public StringBuffer getPfamAccession() {
		return pfamAccession;
	}

	public StringBuffer getLocation() {
		return location;
	}

	public StringBuffer getWikipediaLink() {
		return wikipediaLink;
	}

	public StringBuffer getClan() {
		return clan;
	}

	public StringBuffer getMembership() {
		return membership;
	}

	public Map<String, List<String>> getEmbTrees() {
		return embTrees;
	}

	public void setNumSequences(int numSequences) {
		this.numSequences = numSequences;
	}

	public StringBuffer getIdentification() {
		return identification;
	}

	public void setGFIdentification(String identification) {
		if (this.identification != null) {
			this.identification.append(identification);
		} else {
			this.identification = new StringBuffer(identification);
		}
	}

	public StringBuffer getDefinition() {
		return definition;
	}

	public void setGFDefinition(String definition) {
		if (this.definition != null) {
			this.definition.append(definition);
		} else {
			this.definition = new StringBuffer(definition);
		}
	}

	public String[] getAuthors() {
		return authors;
	}

	public void setGFAuthors(String authors) {
		this.authors = authors.split(",");
	}

	public String getBuildMethod() {
		return buildMethod.toString();
	}

	public void addGFBuildMethod(String buildMethod) {
		if (this.buildMethod == null) {
			this.buildMethod = new StringBuffer(buildMethod);
		} else if (this.buildMethod instanceof StringBuffer){
			((StringBuffer) this.buildMethod).append(System.getProperty("line.seperator")).append(buildMethod);
		}else {
			this.buildMethod = new StringBuffer(this.buildMethod).append(System.getProperty("line.seperator")).append(buildMethod);
		}
	}

	public StringBuffer getSearchMethod() {
		return searchMethod;
	}

	public void setGFSearchMethod(String searchMethod) {
		if (this.searchMethod != null) {
			this.searchMethod.append(searchMethod);
		} else {
			this.searchMethod = new StringBuffer(searchMethod);
		}
	}

	public StringBuffer getSourceSeed() {
		return sourceSeed;
	}

	public void setGFSourceSeed(String sourceSeed) {
		if (this.sourceSeed != null) {
			this.sourceSeed.append(sourceSeed);
		} else {
			this.sourceSeed = new StringBuffer(sourceSeed);
		}
	}

	public StringBuffer getSourceStructure() {
		return sourceStructure;
	}

	public void setGFSourceStructure(String sourceStructure) {
		if (this.sourceStructure != null) {
			this.sourceStructure.append(sourceStructure);
		} else {
			this.sourceStructure = new StringBuffer(sourceStructure);
		}
	}

	/**Not always 2.<br>
	 * It may undergo further change.
	 * @return
	 */
	public float[] getGatheringThreshs() {
		return gatheringThreshs;
	}

	public void setGFGatheringThreshs(String gatheringThresh) {
		this.gatheringThreshs = stringToFloats(gatheringThresh);
	}

	/**Not always 2.<br>
	 * It may undergo further change.
	 * @return
	 */
	public float[] getNoiseCutoffs() {
		return noiseCutoffs;
	}

	public void setGFNoiseCutoffs(String noiseCutoff) {
		this.noiseCutoffs=stringToFloats(noiseCutoff);
	}


	/**Not always 2.<br>
	 * It may undergo further change.
	 * @return
	 */
	public float[] getTrustedCutoffs() {
		return trustedCutoffs;
	}

	public void setGFTrustedCutoffs(String trustedCutoff) {
		this.trustedCutoffs = stringToFloats(trustedCutoff);
	}

	public float[] stringToFloats(String string) {
		String[] coublets= string.split(";");
		float[] floats = new float[coublets.length*2];
		int counter=0;
		for (int i = 0; i < coublets.length; i++) {
			String[] subStrings = coublets[i].trim().split("\\s");
			float f = Float.parseFloat(subStrings[i]);
			floats[counter++]=f;
		}
		return floats;
	}

	public String getTypeField() {
		return typeField;
	}

	public void setGFTypeField(String typeField) {
		this.typeField = typeField;
	}

	public String[] getPreviousIDs() {
		return previousIDs;
	}

	public void setGFPreviousIDs(String previousIDs) {
		this.previousIDs = previousIDs.split(";");
	}

	public StringBuffer getFormat() {
		return format;
	}

	public void setFormat(String format) {
		if (this.format != null) {
			this.format.append(format);
		} else {
			this.format = new StringBuffer(format);
		}
	}

	public StringBuffer getVersion() {
		return version;
	}

	public void setVersion(String version) {
		if (this.version != null) {
			this.version.append(version);
		} else {
			this.version = new StringBuffer(version);
		}
	}

	public String getAccessionNumber() {
		return accessionNumber;
	}

	public void setGFAccessionNumber(String accessionNumber) {
		this.accessionNumber = accessionNumber;
	}

	public boolean isPFam() {
		return accessionNumber != null
				&& (accessionNumber.startsWith(PF_PFAM_STRING) || accessionNumber.startsWith(PB_PFAM_STRING));
	}

	public boolean isRFam() {
		return accessionNumber == null
				|| accessionNumber.startsWith(RF_RFAM_STRING);
	}

	public int getNumSequences() {
		return numSequences;
	}

	public void setGFNumSequences(String numSequences) {
		this.numSequences = Integer.parseInt(numSequences);
	}

	public void setGFDBComment(String dbComment) {
		if (this.dbComment != null) {
			this.dbComment.append(dbComment);
		} else {
			this.dbComment = new StringBuffer(dbComment);
		}
	}


	public void setGFRefComment(String refComment) {
		if (this.refComment != null) {
			this.refComment.append(refComment);
		} else {
			this.refComment = new StringBuffer(refComment);
		}
	}

	public void setGFKeywords(String keywords) {
		if (this.keywords != null) {
			this.keywords.append(keywords);
		} else {
			this.keywords = new StringBuffer(keywords);
		}
	}

	public void addToGFComment(String comment) {
		if (this.comment == null) {
			this.comment = new StringBuffer(comment);
		} else if (this.comment instanceof StringBuffer){
			((StringBuffer) this.comment).append(' ').append(comment);
		}else {
			this.comment = new StringBuffer(this.comment).append(' ').append(comment);
		}
	}

	public void setGFPfamAccession(String pfamAccession) {
		if (this.pfamAccession != null) {
			this.pfamAccession.append(pfamAccession);
		} else {
			this.pfamAccession = new StringBuffer(pfamAccession);
		}
	}

	public void setGFLocation(String location) {
		if (this.location != null) {
			this.location.append(location);
		} else {
			this.location = new StringBuffer(location);
		}
	}

	public void setGFWikipediaLink(String wikipediaLink) {
		if (this.wikipediaLink != null) {
			this.wikipediaLink.append(wikipediaLink);
		} else {
			this.wikipediaLink = new StringBuffer(wikipediaLink);
		}
	}

	public void setGFClan(String clan) {
		if (this.clan != null) {
			this.clan.append(clan);
		} else {
			this.clan = new StringBuffer(clan);
		}
	}

	public void setGFMembership(String membership) {
		if (this.membership != null) {
			this.membership.append(membership);
		} else {
			this.membership = new StringBuffer(membership);
		}
	}

	public void addGFNewHampshire(String newHampshire) {
		List<String> hampshireTree = embTrees.get(TREE_DEFAULT_ID);
		if (hampshireTree == null) {
			hampshireTree = new ArrayList<String>();
		}
		hampshireTree.add(newHampshire);
		embTrees.put(TREE_DEFAULT_ID, hampshireTree);
	}

	public void addGFTreeID(String treeID) {
		List<String> hampshireTree = embTrees.remove(TREE_DEFAULT_ID);
		embTrees.put(treeID, hampshireTree);
	}

	public void addGFFalseDiscoveryRate(String falseDiscoveryRate) {
		this.falseDiscoveryRate = Float.parseFloat(falseDiscoveryRate);
	}

	public String getAlignmentMethod() {
		return alignmentMethod;
	}

	public void setAlignmentMethod(String alignmentMethod) {
		this.alignmentMethod = alignmentMethod;
	}

	public Vector<StockholmFileAnnotationReference> getReferences() {
		return references;
	}

	public void setReferences(Vector<StockholmFileAnnotationReference> references) {
		this.references = references;
	}

}