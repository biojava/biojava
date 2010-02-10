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

package org.biojavax.bio.seq.io;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;

import org.biojava.bio.seq.Feature.Template;
import org.biojava.bio.seq.io.ParseException;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.Symbol;
import org.biojavax.Namespace;
import org.biojavax.RankedCrossRef;
import org.biojavax.RankedDocRef;
import org.biojavax.bio.BioEntryRelationship;
import org.biojavax.bio.seq.RichFeature;
import org.biojavax.bio.taxa.NCBITaxon;

/**
 * This is purely for debugging purposes. Use it to wrap an input stream
 * then pass it as a parameter to some format to listen for
 * sequence generation events. It will dump out char-by-char what it reads,
 * followed event-by-event the events it receives from the format.
 * @author Richard Holland
 * @author Mark Schreiber
 * @author Michael Heuer
 * @since 1.5
 */
public class DebuggingRichSeqIOListener extends BufferedInputStream implements
		RichSeqIOListener {

	private RichFeature currentFeature;
	
	public DebuggingRichSeqIOListener(InputStream inputStream) {
		super(inputStream);
		this.message("Beginning debug session.");
	}
	
	/**
	 * {@inheritDoc}
	 * <p>
	 * Echoes everything that is read to stdout.
	 */
	public int read() throws IOException {
		int ch = super.read();
		System.out.print(ch);
		return ch;
	}

	public void setAccession(String accession) throws ParseException {
		this.message("setAccession: "+accession);
	}

	public void setIdentifier(String identifier) throws ParseException {
		this.message("setIdentifier: "+identifier);
	}

	public void setDivision(String division) throws ParseException {
		this.message("setDivision: "+division);
	}

	public void setDescription(String description) throws ParseException {
		this.message("setDescription: "+description);
	}

	public void setVersion(int version) throws ParseException {
		this.message("setVersion: "+version);
	}

	public void setSeqVersion(String version) throws ParseException {
		this.message("setSeqVersion: "+version);
	}

	public void setComment(String comment) throws ParseException {
		this.message("setComment: "+comment);
	}

	public void setRankedDocRef(RankedDocRef ref) throws ParseException {
		this.message("setRankedDocRef: "+ref);
	}

	public void setTaxon(NCBITaxon taxon) throws ParseException {
		this.message("setTaxon: "+taxon);
	}

	public void setNamespace(Namespace namespace) throws ParseException {
		this.message("setNamespace: "+namespace);
	}

	public void setRelationship(BioEntryRelationship relationship)
			throws ParseException {
		this.message("setRelationship: "+relationship);
	}

	public void setRankedCrossRef(RankedCrossRef crossRef)
			throws ParseException {
		this.message("setRankedCrossRef: "+crossRef);
	}

	public void setURI(String uri) throws ParseException {
		this.message("setURI: "+uri);
	}

	public RichFeature getCurrentFeature() throws ParseException {
		this.message("Current feature requested.");
		return this.currentFeature;
	}

	public void setCircular(boolean circular) throws ParseException {
		this.message("setCircular: "+circular);
	}

	public void startSequence() throws ParseException {
		this.message("startSequence");
	}

	public void endSequence() throws ParseException {
		this.message("endSequence");
	}

	public void setName(String name) throws ParseException {
		this.message("setName: "+name);
	}

	public void addSymbols(Alphabet alpha, Symbol[] syms, int start, int length)
			throws IllegalAlphabetException {
		this.message("addSymbols: (alpha) "+alpha
				+" (syms) "+syms
				+" (start) "+start
				+" (length) "+length);
	}

	public void addSequenceProperty(Object key, Object value)
			throws ParseException {
		this.message("addSequenceProperty: (key) "+key+" (value) "+value);
	}

	public void startFeature(Template templ) throws ParseException {
		this.message("startFeature: (templ) "+templ);
		this.currentFeature = RichFeature.Tools.makeEmptyFeature();
	}

	public void endFeature() throws ParseException {
		this.message("endFeature");
	}

	public void addFeatureProperty(Object key, Object value)
			throws ParseException {
		this.message("addFeatureProperty: (key) "+key+" (value) "+value);
	}
	
	private void message(String message) {
		System.out.println("\n\n#####\n"+message+"\n#####\n\n");
		System.out.println("##### READING... #####\n\n");
	}
}
