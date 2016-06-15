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
package org.biojava.nbio.core.sequence.io;

import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;
import org.biojava.nbio.core.sequence.features.FeatureInterface;
import org.biojava.nbio.core.sequence.io.template.GenbankHeaderFormatInterface;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.util.StringManipulationHelper;

import java.text.SimpleDateFormat;
import java.util.*;

public class GenericGenbankHeaderFormat<S extends AbstractSequence<C>, C extends Compound>
		extends GenericInsdcHeaderFormat<S, C> implements
		GenbankHeaderFormatInterface<S, C> {
	private static final int HEADER_WIDTH = 12;
	private static final String lineSep = "%n";
	private String seqType = null;

	public GenericGenbankHeaderFormat() {
		seqType = null;
	}

	public GenericGenbankHeaderFormat(String seqType) {
		this.seqType = seqType;
	}

	/**
	 * Used in the the 'header' of each GenBank record.
	 *
	 * @param tag
	 * @param text
	 */
	private String _write_single_line(String tag, String text) {
		assert tag.length() < HEADER_WIDTH;
		return StringManipulationHelper.padRight(tag, HEADER_WIDTH)
				+ text.replace('\n', ' ') + lineSep;
	}

	/**
	 * Used in the the 'header' of each GenBank record.
	 *
	 * @param tag
	 * @param text
	 */
	private String _write_multi_line(String tag, String text) {
		if (text == null) {
			text = "";
		}
		int max_len = MAX_WIDTH - HEADER_WIDTH;
		ArrayList<String> lines = _split_multi_line(text, max_len);
		String output = _write_single_line(tag, lines.get(0));
		for (int i = 1; i < lines.size(); i++) {
			output += _write_single_line("", lines.get(i));
		}
		return output;
	}

	/**
	 * used for DBLINK and any similar later line types. If the list of strings
	 * is empty, nothing is written.
	 *
	 * @param tag
	 * @param text_list
	 */
	/*
	 * private String _write_multi_entries(String tag, ArrayList<String>
	 * text_list) { String output = _write_single_line(tag,text_list.remove(0));
	 * for(String s : text_list) { output += _write_single_line("", s); } return
	 * output; }
	 */

	private String _get_date(S sequence) {
		Date sysdate = Calendar.getInstance().getTime();

		// String default_date =
		// sysdate.get(Calendar.DAY_OF_MONTH)+"-"+sysdate.get(Calendar.MONTH)+"-"+sysdate.get(Calendar.YEAR);
		String default_date = new SimpleDateFormat("dd-MMM-yyyy")
				.format(sysdate);
		return default_date;
		/*
		 * try : date = record.annotations["date"] except KeyError : return
		 * default #Cope with a list of one string: if isinstance(date, list)
		 * and len(date)==1 : date = date[0] #TODO - allow a Python date object
		 * if not isinstance(date, str) or len(date) != 11 \ or date[2] != "-"
		 * or date[6] != "-" \ or not date[:2].isdigit() or not
		 * date[7:].isdigit() \ or int(date[:2]) > 31 \ or date[3:6] not in
		 * ["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP",
		 * "OCT", "NOV", "DEC"] : #TODO - Check is a valid date (e.g. not 31
		 * Feb) return default return date
		 */
	}

	private String _get_data_division(S sequence) {
		return UNKNOWN_DNA;
		/*
		 * try: division = record.annotations["data_file_division"] except
		 * KeyError: division = "UNK" if division in ["PRI", "ROD", "MAM",
		 * "VRT", "INV", "PLN", "BCT", "VRL", "PHG", "SYN", "UNA", "EST", "PAT",
		 * "STS", "GSS", "HTG", "HTC", "ENV", "CON"]: #Good, already GenBank
		 * style # PRI - primate sequences # ROD - rodent sequences # MAM -
		 * other mammalian sequences # VRT - other vertebrate sequences # INV -
		 * invertebrate sequences # PLN - plant, fungal, and algal sequences #
		 * BCT - bacterial sequences [plus archea] # VRL - viral sequences # PHG
		 * - bacteriophage sequences # SYN - synthetic sequences # UNA -
		 * unannotated sequences # EST - EST sequences (expressed sequence tags)
		 * # PAT - patent sequences # STS - STS sequences (sequence tagged
		 * sites) # GSS - GSS sequences (genome survey sequences) # HTG - HTGS
		 * sequences (high throughput genomic sequences) # HTC - HTC sequences
		 * (high throughput cDNA sequences) # ENV - Environmental sampling
		 * sequences # CON - Constructed sequences # #(plus UNK for unknown)
		 * pass else: #See if this is in EMBL style: # Division Code #
		 * ----------------- ---- # Bacteriophage PHG - common # Environmental
		 * Sample ENV - common # Fungal FUN - map to PLN (plants + fungal) #
		 * Human HUM - map to PRI (primates) # Invertebrate INV - common # Other
		 * Mammal MAM - common # Other Vertebrate VRT - common # Mus musculus
		 * MUS - map to ROD (rodent) # Plant PLN - common # Prokaryote PRO - map
		 * to BCT (poor name) # Other Rodent ROD - common # Synthetic SYN -
		 * common # Transgenic TGN - ??? map to SYN ??? # Unclassified UNC - map
		 * to UNK # Viral VRL - common # #(plus XXX for submiting which we can
		 * map to UNK) embl_to_gbk = {"FUN":"PLN", "HUM":"PRI", "MUS":"ROD",
		 * "PRO":"BCT", "UNC":"UNK", "XXX":"UNK", } try: division =
		 * embl_to_gbk[division] except KeyError: division = "UNK" assert
		 * len(division)==3 return division
		 */
	}

	/**
	 * Write the LOCUS line.
	 *
	 * @param sequence
	 * @param seqType
	 */
	private String _write_the_first_line(S sequence) {
		/*
		 * locus = record.name if not locus or locus == "<unknown name>": locus
		 * = record.id if not locus or locus == "<unknown id>": locus =
		 * self._get_annotation_str(record, "accession", just_first=True)\
		 */
		String locus;
		try {
			locus = sequence.getAccession().getID();
		} catch (Exception e) {
			locus = "";
		}
		if (locus.length() > 16) {
			throw new RuntimeException("Locus identifier " + locus
					+ " is too long");
		}

		String units = "";
		String mol_type = "";
		if (sequence.getCompoundSet() instanceof DNACompoundSet) {
			units = "bp";
			mol_type = "DNA";
		} else if (sequence.getCompoundSet() instanceof DNACompoundSet) {
			units = "bp";
			mol_type = "RNA";
		} else if (sequence.getCompoundSet() instanceof AminoAcidCompoundSet) {
			units = "aa";
			mol_type = "";
		} else {
			throw new RuntimeException(
					"Need a DNACompoundSet, RNACompoundSet, or an AminoAcidCompoundSet");
		}

		String division = _get_data_division(sequence);

		if (seqType != null) {
			division = seqType;
		}
		assert units.length() == 2;

		// the next line does not seem right.. seqType == linear
		// uncommenting for now
		//assert division.length() == 3;

		StringBuilder sb = new StringBuilder();
		Formatter formatter = new Formatter(sb, Locale.US);
		formatter
				.format("LOCUS       %s %s %s    %s           %s %s" + lineSep,
						StringManipulationHelper.padRight(locus, 16),
						StringManipulationHelper.padLeft(
								Integer.toString(sequence.getLength()), 11),
						units, StringManipulationHelper.padRight(mol_type, 6), division,
						_get_date(sequence));
		String output = formatter.toString();
		formatter.close();
		return output;
		/*
		 * assert len(line) == 79+1, repr(line) #plus one for new line
		 *
		 * assert line[12:28].rstrip() == locus, \ 'LOCUS line does not contain
		 * the locus at the expected position:\n' + line assert line[28:29] ==
		 * " " assert line[29:40].lstrip() == str(len(record)), \ 'LOCUS line
		 * does not contain the length at the expected position:\n' + line
		 *
		 * #Tests copied from Bio.GenBank.Scanner assert line[40:44] in [' bp ',
		 * ' aa '] , \ 'LOCUS line does not contain size units at expected
		 * position:\n' + line assert line[44:47] in [' ', 'ss-', 'ds-', 'ms-'],
		 * \ 'LOCUS line does not have valid strand type (Single stranded,
		 * ...):\n' + line assert line[47:54].strip() == "" \ or
		 * line[47:54].strip().find('DNA') != -1 \ or
		 * line[47:54].strip().find('RNA') != -1, \ 'LOCUS line does not contain
		 * valid sequence type (DNA, RNA, ...):\n' + line assert line[54:55] ==
		 * ' ', \ 'LOCUS line does not contain space at position 55:\n' + line
		 * assert line[55:63].strip() in ['', 'linear', 'circular'], \ 'LOCUS
		 * line does not contain valid entry (linear, circular, ...):\n' + line
		 * assert line[63:64] == ' ', \ 'LOCUS line does not contain space at
		 * position 64:\n' + line assert line[67:68] == ' ', \ 'LOCUS line does
		 * not contain space at position 68:\n' + line assert line[70:71] ==
		 * '-', \ 'LOCUS line does not contain - at position 71 in date:\n' +
		 * line assert line[74:75] == '-', \ 'LOCUS line does not contain - at
		 * position 75 in date:\n' + line
		 */
	}

	/**
	 * This is a bit complicated due to the range of possible ways people might
	 * have done their annotation... Currently the parser uses a single string
	 * with newlines. A list of lines is also reasonable. A single (long) string
	 * is perhaps the most natural of all. This means we may need to deal with
	 * line wrapping.
	 *
	 * @param sequence
	 */
	private String _write_comment(S sequence) {
		ArrayList<String> comments = sequence.getNotesList();
		String output = _write_multi_line("COMMENT", comments.remove(0));
		for (String comment : comments) {
			output += _write_multi_line("", comment);
		}

		return output;
	}

	@Override
	public String getHeader(S sequence) {
		String header = _write_the_first_line(sequence);
		String acc_with_version;
		String accession;
		try {
			acc_with_version = sequence.getAccession().getID();
			accession = acc_with_version.split("\\.", 1)[0];
		} catch (Exception e) {
			acc_with_version = "";
			accession = "";
		}
		String description = sequence.getDescription();
		if ("<unknown description>".equals(description) || description == null) {
			description = ".";
		}
		header += _write_multi_line("DEFINITION", description);
		header += _write_multi_line("ACCESSION", accession);
		header += _write_multi_line("VERSION", acc_with_version);

		/*
		 * gi = self._get_annotation_str(record, "gi", just_first=True)
		 *
		 * self._write_single_line("ACCESSION", accession) if gi != ".":
		 * self._write_single_line("VERSION", "%s  GI:%s" \ % (acc_with_version,
		 * gi)) else: self._write_single_line("VERSION", "%s" %
		 * (acc_with_version))
		 *
		 * #The NCBI only expect two types of link so far, #e.g. "Project:28471"
		 * and "Trace Assembly Archive:123456" #TODO - Filter the dbxrefs list
		 * to just these? self._write_multi_entries("DBLINK", record.dbxrefs)
		 *
		 * try: #List of strings #Keywords should be given separated with semi
		 * colons, keywords = "; ".join(record.annotations["keywords"]) #with a
		 * trailing period: if not keywords.endswith(".") : keywords += "."
		 * except KeyError: #If no keywords, there should be just a period:
		 * keywords = "."
		 */

		header += _write_multi_line("KEYWORDS", ".");

		/*
		 * if "segment" in record.annotations: #Deal with SEGMENT line found
		 * only in segmented records, #e.g. AH000819 segment =
		 * record.annotations["segment"] if isinstance(segment, list): assert
		 * len(segment)==1, segment segment = segment[0]
		 * self._write_single_line("SEGMENT", segment)
		 *
		 * self._write_multi_line("SOURCE", \ self._get_annotation_str(record,
		 * "source"))
		 */

		header += _write_multi_line("SOURCE", sequence.getSource());

		/*
		 * #The ORGANISM line MUST be a single line, as any continuation is the
		 * taxonomy org = self._get_annotation_str(record, "organism") if
		 * len(org) > self.MAX_WIDTH - self.HEADER_WIDTH: org =
		 * org[:self.MAX_WIDTH - self.HEADER_WIDTH-4]+"..."
		 * self._write_single_line("  ORGANISM", org) try: #List of strings
		 * #Taxonomy should be given separated with semi colons, taxonomy =
		 * "; ".join(record.annotations["taxonomy"]) #with a trailing period: if
		 * not taxonomy.endswith(".") : taxonomy += "." except KeyError:
		 * taxonomy = "." self._write_multi_line("", taxonomy)
		 *
		 * if "references" in record.annotations: self._write_references(record)
		 */
		if (!sequence.getNotesList().isEmpty()) {
			header += _write_comment(sequence);
		}

		header += "FEATURES             Location/Qualifiers" + lineSep;
		int rec_length = sequence.getLength();
		for (FeatureInterface<AbstractSequence<C>, C> feature : sequence
				.getFeatures()) {
			header += _write_feature(feature, rec_length);
		}

		return header;
	}

}
