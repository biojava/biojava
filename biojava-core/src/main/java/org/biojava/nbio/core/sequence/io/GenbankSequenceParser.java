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
 * @author Richard Holland
 * @author Mark Schreiber
 * @author David Scott
 * @author Bubba Puryear
 * @author George Waldon
 * @author Deepak Sheoran
 * @author Karl Nicholas <github:karlnicholas>
 * @author Jacek Grzebyta
 * @author Paolo Pavan
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on 01-21-2010
 */
package org.biojava.nbio.core.sequence.io;

import org.biojava.nbio.core.exceptions.Messages;
import org.biojava.nbio.core.exceptions.ParserException;
import org.biojava.nbio.core.sequence.DataSource;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;
import org.biojava.nbio.core.sequence.compound.RNACompoundSet;
import org.biojava.nbio.core.sequence.features.AbstractFeature;
import org.biojava.nbio.core.sequence.features.DBReferenceInfo;
import org.biojava.nbio.core.sequence.features.Qualifier;
import org.biojava.nbio.core.sequence.features.TextFeature;
import org.biojava.nbio.core.sequence.io.template.SequenceParserInterface;
import org.biojava.nbio.core.sequence.location.InsdcParser;
import org.biojava.nbio.core.sequence.location.template.AbstractLocation;
import org.biojava.nbio.core.sequence.location.template.Location;
import org.biojava.nbio.core.sequence.reference.GenbankReference;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.CompoundSet;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class GenbankSequenceParser<S extends AbstractSequence<C>, C extends Compound> implements SequenceParserInterface{

	private String seqData = null;
	private GenericGenbankHeaderParser<S, C> headerParser;
	private String header;
	private String accession;
	private boolean isCircularSequence;
	private Map<String, List<DBReferenceInfo>> mapDB;
	/**
	 * this data structure collects list of features extracted from the
	 * FEATURE_TAG section They are organized by list of the same type (i.e.
	 * same genbank Feature) and are provided with location
	 */
	private Map<String, List<AbstractFeature<AbstractSequence<C>, C>>> featureCollection;

	private final Logger log = LoggerFactory.getLogger(getClass());

	// this is a compoundset parsed from header.
	private CompoundSet<?> compoundType;

	/**
	 * The name of this format
	 */
	public static final String GENBANK_FORMAT = "GENBANK";

	protected static final String LOCUS_TAG = "LOCUS";
	protected static final String DEFINITION_TAG = "DEFINITION";
	protected static final String ACCESSION_TAG = "ACCESSION";
	protected static final String VERSION_TAG = "VERSION";
	protected static final String KEYWORDS_TAG = "KEYWORDS";
	//                                                  "SEGMENT"
	protected static final String SOURCE_TAG = "SOURCE";
	protected static final String ORGANISM_TAG = "ORGANISM";
	protected static final String REFERENCE_TAG = "REFERENCE";
	protected static final String AUTHORS_TAG = "AUTHORS";
	protected static final String CONSORTIUM_TAG = "CONSRTM";
	protected static final String TITLE_TAG = "TITLE";
	protected static final String JOURNAL_TAG = "JOURNAL";
	protected static final String PUBMED_TAG = "PUBMED";
	protected static final String MEDLINE_TAG = "MEDLINE"; //deprecated
	protected static final String REMARK_TAG = "REMARK";
	protected static final String COMMENT_TAG = "COMMENT";
	protected static final String FEATURE_TAG = "FEATURES";
	protected static final String BASE_COUNT_TAG_FULL = "BASE COUNT"; //deprecated
	protected static final String BASE_COUNT_TAG = "BASE";
	//                                                  "CONTIG"
	protected static final String START_SEQUENCE_TAG = "ORIGIN";
	protected static final String DBSOURCE = "DBSOURCE";
	protected static final String PRIMARY = "PRIMARY";
	protected static final String DBLINK = "DBLINK";
	protected static final String END_SEQUENCE_TAG = "//";
	// locus line
	protected static final Pattern lp = Pattern.compile("^(\\S+)\\s+(\\d+)\\s+(bp|BP|aa|AA)\\s{0,4}(([dmsDMS][sS]-)?(\\S+))?\\s*(circular|CIRCULAR|linear|LINEAR)?\\s*(\\S+)?\\s*(\\S+)?$");
	// version line
	protected static final Pattern vp = Pattern.compile("^(\\S*?)(\\.(\\d+))?(\\s+GI:(\\S+))?$");
	// reference line
	protected static final Pattern refRange = Pattern.compile("^\\s*(\\d+)\\s+to\\s+(\\d+)$");
	protected static final Pattern refp = Pattern.compile("^(\\d+)\\s*(?:(\\((?:bases|residues)\\s+(\\d+\\s+to\\s+\\d+(\\s*;\\s*\\d+\\s+to\\s+\\d+)*)\\))|\\(sites\\))?");
	// dbxref line
	protected static final Pattern dbxp = Pattern.compile("^([^:]+):(\\S+)$");

	protected static final InsdcParser locationParser = new InsdcParser(DataSource.GENBANK);
	/**
	 * sections start at a line and continue till the first line afterwards with a
	 * 	non-whitespace first character
	 * 	we want to match any of the following as a new section within a section
	 * 	  \s{0,8} word \s{0,7} value
	 * 	  \s{21} /word = value
	 * 	  \s{21} /word
	 */
	protected static final Pattern sectp = Pattern.compile("^(\\s{0,8}(\\S+)\\s{0,7}(.*)|\\s{21}(/\\S+?)=(.*)|\\s{21}(/\\S+))$");

	protected static final Pattern readableFiles = Pattern.compile(".*(g[bp]k*$|\\u002eg[bp].*)");
	protected static final Pattern headerLine = Pattern.compile("^LOCUS.*");


	private String parse(BufferedReader bufferedReader) {
		String sectionKey;
		List<String[]> section;
		// Get an ordered list of key->value pairs in array-tuples
		do {
			section = this.readSection(bufferedReader);
			sectionKey = section.get(0)[0];
			if (sectionKey == null) {
				//if we reach the end of the file, section contains empty strings
				if(section.get(0)[1]==null || section.get(0)[1].equals("") ||
						section.get(0)[1].length()==0) {
					throw new ParserException(Messages.ENDOFFILE);
				}
				throw new ParserException(Messages.SECTIONKEYNULL);
			}
			// process section-by-section
			switch (sectionKey) {
				case LOCUS_TAG: parseLocusTag(section); break;
				case DEFINITION_TAG: parseDefinitionTag(section); break;
				case ACCESSION_TAG: parseAccessionTag(section); break;
				case VERSION_TAG: parseVersionTag(section); break;
				case KEYWORDS_TAG: break; 	// not implemented yet
				case SOURCE_TAG: break; 	// ignore - can get all this from the first feature
				case REFERENCE_TAG: parseReferenceTag(section); break;
				case COMMENT_TAG: parseCommentTag(section); break;
				case FEATURE_TAG: parseFeatureTag(section); break;
				case BASE_COUNT_TAG: break;	// ignore - can calculate from sequence content later if needed
				case START_SEQUENCE_TAG: parseStartSequenceTag(section); break;
				case DBSOURCE: break;		// not implemented yet
				case PRIMARY: break;		// not implemented yet
				case DBLINK: break;			// not implemented yet
				default:
					if(!sectionKey.equals(END_SEQUENCE_TAG)) {
						log.info("found unknown section key: %", sectionKey);
					}
			}
		} while (!sectionKey.equals(END_SEQUENCE_TAG));
		return seqData;
	}

	private void parseStartSequenceTag(List<String[]> section) {
		// our first line is ignorable as it is the ORIGIN tag
		// the second line onwards conveniently have the number as
		// the [0] tuple, and sequence string as [1] so all we have
		// to do is concat the [1] parts and then strip out spaces,
		// and replace '.' and '~' with '-' for our parser.
		StringBuilder seq = new StringBuilder();
		for (int i = 1; i < section.size(); i++) {
			seq.append(section.get(i)[1]);
		}
		seqData = seq.toString().replaceAll("\\s+", "").replaceAll("[\\.|~]", "-").toUpperCase();
	}

	private void parseFeatureTag(List<String[]> section) {
		// starting from second line of input, start a new feature whenever we come across
		// a key that does not start with /
		AbstractFeature gbFeature = null;
		for (int i = 1; i < section.size(); i++) {
			String key = section.get(i)[0];
			String val = section.get(i)[1];
			if (key.startsWith("/")) {
				if (gbFeature == null) {
					throw new ParserException("Malformed GenBank file: found a qualifier without feature.");
				}
				key = key.substring(1); // strip leading slash
				val = val.replaceAll("\\s*[\\n\\r]+\\s*", " ").trim();
				if (val.endsWith("\"")) {
					val = val.substring(1, val.length() - 1); // strip quotes
				}
				// parameter on old feature
				if (key.equals("db_xref")) {
					Matcher m = dbxp.matcher(val);
					if (m.matches()) {
						String dbname = m.group(1);
						String raccession = m.group(2);
						DBReferenceInfo xref = new DBReferenceInfo(dbname, raccession);
						gbFeature.addQualifier(key, xref);

						ArrayList<DBReferenceInfo> listDBEntry = new ArrayList<>();
						listDBEntry.add(xref);
						mapDB.put(key, listDBEntry);
					} else {
						throw new ParserException("Bad dbxref");
					}
				} else if (key.equalsIgnoreCase("organism")) {
					Qualifier q = new Qualifier(key, val.replace('\n', ' '));
					gbFeature.addQualifier(key, q);
				} else {
					if (key.equalsIgnoreCase("translation") || key.equals("anticodon")
							|| key.equals("transl_except")) {
						// strip spaces from sequence
						val = val.replaceAll("\\s+", "");
						Qualifier q = new Qualifier(key, val);
						gbFeature.addQualifier(key, q);
					} else {
						Qualifier q = new Qualifier(key, val);
						gbFeature.addQualifier(key, q);
					}
				}
			} else {
				// new feature!
				gbFeature = new TextFeature(key, val, key, key);
				Location l =
						locationParser.parse(val);
				gbFeature.setLocation((AbstractLocation)l);

				if (!featureCollection.containsKey(key)) {
					featureCollection.put(key, new ArrayList<>());
				}
				featureCollection.get(key).add(gbFeature);
			}
		}
	}

	private void parseCommentTag(List<String[]> section) {
		headerParser.setComment(section.get(0)[1]);
	}

	private void parseReferenceTag(List<String[]> section) {
		GenbankReference genbankReference = new GenbankReference();
		for (String[] ref : section) {
			if (ref[0].equals(AUTHORS_TAG)) {
				genbankReference.setAuthors(ref[1]);
			} else if (ref[0].equals(TITLE_TAG)) {
				genbankReference.setTitle(ref[1]);
			} else if (ref[0].equals(JOURNAL_TAG)) {
				genbankReference.setJournal(ref[1]);
			}
		}
		headerParser.addReference(genbankReference);
	}

	private void parseVersionTag(List<String[]> section) {
		String ver = section.get(0)[1];
		Matcher m = vp.matcher(ver);
		if (m.matches()) {
			String verAcc = m.group(1);
			if (!accession.equals(verAcc)) {
				// the version refers to a different accession!
				// believe the version line, and store the original
				// accession away in the additional accession set
				accession = verAcc;
			}
			if (m.group(3) != null) {
				headerParser.setVersion(Integer.parseInt(m.group(3)));
			}
			if (m.group(5) != null) {
				headerParser.setIdentifier(m.group(5));
			}
		} else {
			throw new ParserException("Bad version line");
		}
	}

	private void parseAccessionTag(List<String[]> section) {
		// if multiple accessions, store only first as accession,
		// and store rest in annotation
		String[] accs = section.get(0)[1].split("\\s+");
		accession = accs[0].trim();
		headerParser.setAccession(accession);
	}

	private void parseDefinitionTag(List<String[]> section) {
		headerParser.setDescription(section.get(0)[1]);
	}

	private void parseLocusTag(List<String[]> section) {
		String loc = section.get(0)[1];
		header = loc;
		Matcher m = lp.matcher(loc);
		if (m.matches()) {
			headerParser.setName(m.group(1));
			headerParser.setAccession(m.group(1)); // default if no accession found
			long sequenceLength = Long.valueOf(m.group(2));
			String lengthUnits = m.group(3);
			String type = m.group(6);

			if (lengthUnits.equalsIgnoreCase("aa")) {
				compoundType = AminoAcidCompoundSet.getAminoAcidCompoundSet();
			} else if (lengthUnits.equalsIgnoreCase("bp")) {
				if (type != null) {
					if (type.contains("RNA")) {
						compoundType = RNACompoundSet.getRNACompoundSet();
					} else {
						compoundType = DNACompoundSet.getDNACompoundSet();
					}
				} else {
					compoundType = DNACompoundSet.getDNACompoundSet();
				}
			}

			if (m.group(7) != null) isCircularSequence = m.group(7).equalsIgnoreCase("circular");

			// configure location parser with needed information
			locationParser.setSequenceLength(sequenceLength);
			locationParser.setSequenceCircular(isCircularSequence);

			log.debug("compound type: {}", compoundType.getClass().getSimpleName());

		} else {
			throw new ParserException("Bad locus line");
		}
	}


	// reads an indented section, combining split lines and creating a list of
	// key->value tuples
	// reads an indented section, combining split lines and creating a list of
	// key->value tuples
	// reads an indented section, combining split lines and creating a list of
	// key->value tuples
	private List<String[]> readSection(BufferedReader bufferedReader) {
		List<String[]> section = new ArrayList<>();
		String line;

		String currKey = null;
		StringBuilder currVal = new StringBuilder();
		boolean done = false;
		int linecount = 0;

		try {
			while (!done) {
				bufferedReader.mark(320);
				line = bufferedReader.readLine();
				String firstSecKey = section.isEmpty() ? ""
						: section.get(0)[0];
				if (line != null && line.matches("\\p{Space}*")) {
					// regular expression \p{Space}* will match line
					// having only white space characters
					continue;
				}
				if (line == null
						|| (!line.startsWith(" ") && linecount++ > 0 && (!firstSecKey
						.equals(START_SEQUENCE_TAG) || line
						.startsWith(END_SEQUENCE_TAG)))) {
					// dump out last part of section
					section.add(new String[]{currKey, currVal.toString()});
					bufferedReader.reset();
					done = true;
				} else {
					Matcher m = sectp.matcher(line);
					if (m.matches()) {
						// new key
						if (currKey != null) {
							section.add(new String[]{currKey,
								currVal.toString()});
						}
						// key = group(2) or group(4) or group(6) - whichever is
						// not null
						currKey = m.group(2) == null ? (m.group(4) == null ? m
								.group(6) : m.group(4)) : m.group(2);
						currVal = new StringBuilder();
						// val = group(3) if group(2) not null, group(5) if
						// group(4) not null, "" otherwise, trimmed
						currVal.append((m.group(2) == null ? (m.group(4) == null ? ""
								: m.group(5))
								: m.group(3)).trim());
					} else {
						// concatted line or SEQ START/END line?
						if (line.startsWith(START_SEQUENCE_TAG)
								|| line.startsWith(END_SEQUENCE_TAG)) {
							currKey = line;
						} else {
							currVal.append("\n"); // newline in between lines -
							// can be removed later
							currVal.append(currKey.charAt(0) == '/' ? line
									.substring(21) : line.substring(12));
						}
					}
				}
			}
		} catch (IOException | RuntimeException e) {
			throw new ParserException(e.getMessage());
		}
		return section;
	}

	@Override
	public String getSequence(BufferedReader bufferedReader, int sequenceLength) {
		featureCollection = new HashMap<>();
		mapDB = new LinkedHashMap<>();
		headerParser = new GenericGenbankHeaderParser<>();
		try {
			parse(bufferedReader);
		} catch (ParserException e) {
			if(e.getMessage().equalsIgnoreCase(Messages.ENDOFFILE))	return null;
			else throw new ParserException(e.getMessage());
		}

		return seqData;
	}

	public String getHeader() {
		return header;
	}

	public GenericGenbankHeaderParser<S, C> getSequenceHeaderParser() {
		return headerParser;
	}

	public Map<String, List<DBReferenceInfo>> getDatabaseReferences() {
		return mapDB;
	}

	public List<String> getKeyWords() {
		return new ArrayList<>(featureCollection.keySet());
	}

	public List<AbstractFeature<AbstractSequence<C>, C>> getFeatures(String keyword) {
		return featureCollection.get(keyword);
	}
	public Map<String, List<AbstractFeature<AbstractSequence<C>, C>>> getFeatures() {
		return featureCollection;
	}

	public void parseFeatures(AbstractSequence<C> sequence) {
		for (String k: featureCollection.keySet())
			for (AbstractFeature<AbstractSequence<C>, C> f: featureCollection.get(k))
				sequence.addFeature(f);
	}

	public CompoundSet<?> getCompoundType() {
		return compoundType;
	}
}
