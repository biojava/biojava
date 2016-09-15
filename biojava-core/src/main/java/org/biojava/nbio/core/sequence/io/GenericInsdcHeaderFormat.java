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
/**
 *
 */
package org.biojava.nbio.core.sequence.io;

import org.biojava.nbio.core.sequence.Strand;
import org.biojava.nbio.core.sequence.features.FeatureInterface;
import org.biojava.nbio.core.sequence.features.Qualifier;
import org.biojava.nbio.core.sequence.location.template.AbstractLocation;
import org.biojava.nbio.core.sequence.location.template.Point;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.util.StringManipulationHelper;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Formatter;
import java.util.List;
import java.util.Locale;

/**
 * @author mckeee1
 *
 */
public class GenericInsdcHeaderFormat<S extends AbstractSequence<C>, C extends Compound> {
	protected static final int MAX_WIDTH = 80;
	protected static final int QUALIFIER_INDENT = 21;
	protected static final String QUALIFIER_INDENT_STR = "                     ";
	protected static final String QUALIFIER_INDENT_TMP = "     %s                ";
	private static final String lineSep = "%n";

	/**
	 * Format a feature qualifier using the MAX_WIDTH (default 80)
	 * @param key
	 * @param value
	 * @param quote
	 */
	private String _write_feature_qualifier(String key, String value, boolean quote) {
		String line = "";
		if(null == value) {
			line = QUALIFIER_INDENT_STR + "/" + key + lineSep;
			return line;
		}
		if(quote) {  // quote should be true for numerics
			line = QUALIFIER_INDENT_STR + "/" + key + "=\"" + value + "\"";
		} else {
			line = QUALIFIER_INDENT_STR + "/" + key + "=" + value;
		}
		if(line.length() <= MAX_WIDTH) {
			return line + lineSep;
		}
		String goodlines = "";
		while(!"".equals(line.replaceAll("^\\s+", ""))) {
			if(line.length() <= MAX_WIDTH) {
				goodlines += line + lineSep;
				break;
			}
			//Insert line break...
			int index;
			for(index = Math.min(line.length()-1, MAX_WIDTH); index > QUALIFIER_INDENT ; index--) {
				if(' ' == line.charAt(index)) {
					break;
				}
			}
			if(' ' != line.charAt(index)) {
				//no nice place to break...
				index = MAX_WIDTH;
			}
			assert index <= MAX_WIDTH;
			goodlines += line.substring(0,index) + lineSep;
			line = QUALIFIER_INDENT_STR + line.substring(index).replaceAll("^\\s+", "");
		}
		return goodlines;
	}
	/**
	 * Split a feature location into lines (break at commas).
	 * @param location
	 */
	private String _wrap_location(String location) {
		int length = MAX_WIDTH - QUALIFIER_INDENT;
		if(location.length() <= length) {
			return location;
		}
		int index = location.substring(length).lastIndexOf(",");
		if(-1 == index) {
			//No good place to split (!)
			return location;
		}
		return location.substring(0,index+1) + lineSep + QUALIFIER_INDENT_STR + _wrap_location(location.substring(index+1));
	}
	/**
	 * Write a single SeqFeature object to features table.
	 * @param feature
	 * @param record_length
	 */
	protected String _write_feature(FeatureInterface<AbstractSequence<C>, C> feature, int record_length) {
		String location = _insdc_feature_location_string(feature, record_length);
		String f_type = feature.getType().replace(" ", "_");
		StringBuilder sb = new StringBuilder();
		Formatter formatter = new Formatter(sb,Locale.US);
		formatter.format(QUALIFIER_INDENT_TMP, f_type);
		String line = formatter.toString().substring(0, QUALIFIER_INDENT) + _wrap_location(location) + lineSep;
		formatter.close();

		//Now the qualifiers...
		for(List<Qualifier>  qualifiers : feature.getQualifiers().values()) {
			for(Qualifier q : qualifiers){
				line += _write_feature_qualifier(q.getName(), q.getValue(), q.needsQuotes());
			}
		}
		return line;
		/*
		self.handle.write(line)
		#Now the qualifiers...
		for key, values in feature.qualifiers.items():
			if isinstance(values, list) or isinstance(values, tuple):
				for value in values:
					self._write_feature_qualifier(key, value)
			elif values:
				#String, int, etc
				self._write_feature_qualifier(key, values)
			else:
				#e.g. a /psuedo entry
				self._write_feature_qualifier(key)
		 */
	}
	/**
	 * Build a GenBank/EMBL location string from a SeqFeature (PRIVATE).

	There is a choice of how to show joins on the reverse complement strand,
	GenBank used "complement(join(1,10),(20,100))" while EMBL used to use
	"join(complement(20,100),complement(1,10))" instead (but appears to have
	now adopted the GenBank convention). Notice that the order of the entries
	is reversed! This function therefore uses the first form. In this situation
	we expect the parent feature and the two children to all be marked as
	strand == -1, and in the order 0:10 then 19:100.

	Also need to consider dual-strand examples like these from the Arabidopsis
	thaliana chloroplast NC_000932: join(complement(69611..69724),139856..140650)
	gene ArthCp047, GeneID:844801 or its CDS (protein NP_051038.1 GI:7525057)
	which is further complicated by a splice:
	join(complement(69611..69724),139856..140087,140625..140650)

	For mixed this mixed strand feature, the parent SeqFeature should have
	no strand (either 0 or None) while the child features should have either
	strand +1 or -1 as appropriate, and be listed in the order given here.
	 * @param feature
	 * @param record_length
	 */
	private String _insdc_feature_location_string(FeatureInterface<AbstractSequence<C>, C> feature, int record_length) {
		if(feature.getChildrenFeatures().isEmpty()) {
			//Non-recursive.
			String location = _insdc_location_string_ignoring_strand_and_subfeatures(feature.getLocations(), record_length);
			if(feature.getLocations().getStrand() == Strand.NEGATIVE) {
				StringBuilder sb = new StringBuilder();
				Formatter formatter = new Formatter(sb,Locale.US);
				formatter.format("complement(%s)", location);
				String output = formatter.toString();
				formatter.close();
				location = output;
			}
			return location;
		}
		// As noted above, treat reverse complement strand features carefully:
		if(feature.getLocations().getStrand() == Strand.NEGATIVE) {
			for(FeatureInterface<?, ?> f  : feature.getChildrenFeatures()) {
				if(f.getLocations().getStrand() != Strand.NEGATIVE) {
					StringBuilder sb = new StringBuilder();
					Formatter formatter = new Formatter(sb,Locale.US);
					formatter.format("Inconsistent strands: %s for parent, %s for child", feature.getLocations().getStrand(), f.getLocations().getStrand());
					String output = formatter.toString();
					formatter.close();
					throw new RuntimeException(output);
				}
			}
			StringBuilder sb = new StringBuilder();
			Formatter formatter = new Formatter(sb,Locale.US);
			ArrayList<String> locations = new ArrayList<String>();
			for(FeatureInterface<AbstractSequence<C>, C> f  : feature.getChildrenFeatures()) {
				locations.add(_insdc_location_string_ignoring_strand_and_subfeatures(f.getLocations(), record_length));
			}
			String location = StringManipulationHelper.join(locations, ",");
			formatter.format("complement(%s(%s))", /*feature.location_operator*/ "join", location);
			String output = formatter.toString();
			formatter.close();
			return output;
		}
		//This covers typical forward strand features, and also an evil mixed strand:
		StringBuilder sb = new StringBuilder();
		Formatter formatter = new Formatter(sb,Locale.US);
		ArrayList<String> locations = new ArrayList<String>();
		for(FeatureInterface<AbstractSequence<C>, C> f  : feature.getChildrenFeatures()) {
			locations.add(_insdc_location_string_ignoring_strand_and_subfeatures(f.getLocations(), record_length));
		}
		String location =  StringManipulationHelper.join(locations, ",");
		formatter.format("%s(%s)", /*feature.location_operator*/ "join", location);
		String output = formatter.toString();
		formatter.close();
		return output;
	}

	private String _insdc_location_string_ignoring_strand_and_subfeatures(
			//SequenceLocation<AbstractSequence<C>, C> sequenceLocation,
						AbstractLocation sequenceLocation,
			int record_length) {
	/*
	if location.ref:
		ref = "%s:" % location.ref
	else:
		ref = ""
	assert not location.ref_db
	*/
		String ref = "";
		if(!sequenceLocation.getStart().isUncertain() && !sequenceLocation.getEnd().isUncertain() && sequenceLocation.getStart() == sequenceLocation.getEnd()) {
			//Special case, for 12:12 return 12^13
			//(a zero length slice, meaning the point between two letters)
			if(sequenceLocation.getEnd().getPosition() == record_length) {
				//Very special case, for a between position at the end of a
				//sequence (used on some circular genomes, Bug 3098) we have
				//N:N so return N^1
				StringBuilder sb = new StringBuilder();
				Formatter formatter = new Formatter(sb,Locale.US);
				formatter.format("%s%d^1", ref, record_length);
				String output = formatter.toString();
				formatter.close();
				return output;
			} else {
				StringBuilder sb = new StringBuilder();
				Formatter formatter = new Formatter(sb,Locale.US);
				formatter.format("%s%d^%d", ref, sequenceLocation.getStart().getPosition(), sequenceLocation.getEnd().getPosition());
				String output = formatter.toString();
				formatter.close();
				return output;
			}
		}
		if(!sequenceLocation.getStart().isUncertain() && !sequenceLocation.getEnd().isUncertain() && sequenceLocation.getStart().getPosition() + 1 == sequenceLocation.getEnd().getPosition()) {
			//Special case, for 11:12 return 12 rather than 12..12
			//(a length one slice, meaning a single letter)
			StringBuilder sb = new StringBuilder();
			Formatter formatter = new Formatter(sb,Locale.US);
			formatter.format("%s%d", ref, sequenceLocation.getEnd().getPosition());
			String output = formatter.toString();
			formatter.close();
			return output;
		} else if(sequenceLocation.getStart().isUnknown() || sequenceLocation.getEnd().isUnknown()) {
			//Special case for features from SwissProt/UniProt files
			if(sequenceLocation.getStart().isUnknown() && sequenceLocation.getEnd().isUnknown()) {
				throw new RuntimeException("Feature with unknown location");
			} else if(sequenceLocation.getStart().isUnknown()) {
				//Treat the unknown start position as a BeforePosition
				StringBuilder sb = new StringBuilder();
				Formatter formatter = new Formatter(sb,Locale.US);
				formatter.format("%s<%d..%s", ref, sequenceLocation.getEnd().getPosition(), _insdc_feature_position_string(sequenceLocation.getEnd()));
				String output = formatter.toString();
				formatter.close();
				return output;
			} else {
				//Treat the unknown start position as an AfterPosition
				StringBuilder sb = new StringBuilder();
				Formatter formatter = new Formatter(sb,Locale.US);
				formatter.format("%s%s..>%d", ref, _insdc_feature_position_string(sequenceLocation.getStart()), sequenceLocation.getStart().getPosition());
				String output = formatter.toString();
				formatter.close();
				return output;
			}
		} else {
			//Typical case, e.g. 12..15 gets mapped to 11:15
			return ref + _insdc_feature_position_string(sequenceLocation.getStart(), 0) + ".." + _insdc_feature_position_string(sequenceLocation.getEnd());
		}
	}
	private String _insdc_feature_position_string(Point location) {
		// TODO Auto-generated method stub
		return _insdc_feature_position_string(location, 0);
	}

	/**
	 * Build a GenBank/EMBL position string (PRIVATE).
	 * @param location
	 * @param increment
	 */
	private String _insdc_feature_position_string(Point location, int increment) {
			StringBuilder sb = new StringBuilder();
			Formatter formatter = new Formatter(sb,Locale.US);
			formatter.format("%s", location.getPosition() + increment);
			String output = formatter.toString();
			formatter.close();
			return output;

	/*
	if isinstance(pos, SeqFeature.ExactPosition):
		return "%i" % (pos.position+offset)
	elif isinstance(pos, SeqFeature.WithinPosition):
		return "(%i.%i)" % (pos.position + offset,
							pos.position + pos.extension + offset)
	elif isinstance(pos, SeqFeature.BetweenPosition):
		return "(%i^%i)" % (pos.position + offset,
							pos.position + pos.extension + offset)
	elif isinstance(pos, SeqFeature.BeforePosition):
		return "<%i" % (pos.position + offset)
	elif isinstance(pos, SeqFeature.AfterPosition):
		return ">%i" % (pos.position + offset)
	elif isinstance(pos, SeqFeature.OneOfPosition):
		return "one-of(%s)" \
			   % ",".join([_insdc_feature_position_string(p,offset) \
						   for p in pos.position_choices])
	elif isinstance(pos, SeqFeature.AbstractPosition):
		raise NotImplementedError("Please report this as a bug in Biopython.")
	else:
		raise ValueError("Expected a SeqFeature position object.")
		 */
	}

	/**
	 * Returns a list of strings.
	 *
	 *   Any single words which are too long get returned as a whole line
	 *   (e.g. URLs) without an exception or warning.
	 * @param text
	 * @param max_len
	 */
	protected ArrayList<String> _split_multi_line(String text, int max_len) {
		// TODO Auto-generated method stub
		ArrayList<String> output = new ArrayList<String>();
		text = text.trim();
		if(text.length() <= max_len) {
			output.add(text);
			return output;
		}

		ArrayList<String> words = new ArrayList<String>();
		Collections.addAll(words, text.split("\\s+"));
		while(!words.isEmpty()) {
			text = words.remove(0);
			while(!words.isEmpty() && (text.length() + 1 + words.get(0).length()) <= max_len) {
				text += " " + words.remove(0);
				text = text.trim();
			}
			output.add(text);
		}
		assert words.isEmpty();
		return output;
	}


}
