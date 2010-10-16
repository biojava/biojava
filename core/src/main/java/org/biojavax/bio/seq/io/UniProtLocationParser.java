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

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.bio.seq.io.ParseException;
import org.biojavax.RichObjectFactory;
import org.biojavax.bio.seq.Position;
import org.biojavax.bio.seq.RichLocation;
import org.biojavax.bio.seq.SimplePosition;
import org.biojavax.bio.seq.SimpleRichLocation;
import org.biojavax.bio.seq.RichLocation.Strand;
import org.biojavax.utils.StringTools;


/**
 * Parses UniProt location strings into RichLocation objects.
 * @author Richard Holland
 * @since 1.5
 */
public class UniProtLocationParser {

	// No instances please
	private UniProtLocationParser() {}

	/**
	 * Parses a location.
	 * @param loc the UniProt location string.
	 * @return RichLocation the equivalent RichLocation object.
	 * @throws ParseException if the parsing failed.
	 */
	public static RichLocation parseLocation(String loc) throws ParseException {
		try{
			String parts[] = loc.trim().split("\\s+");
			Position startPos = null;
			Position endPos   = null; 

			try {
				startPos = parsePosition(parts[0].trim());
			} catch (Exception e){
				System.err.println(e.getMessage());
			}
			try {
				endPos   = parsePosition(parts[1].trim());
			} catch (Exception e){
				System.err.println(e.getMessage());
			}

			if (( startPos == null) && ( endPos == null)){
				return new SimpleRichLocation(new SimplePosition(0),new SimplePosition(0),1,Strand.POSITIVE_STRAND,null);
			} else if ( endPos == null){
				return new SimpleRichLocation(startPos,new SimplePosition(0),1,Strand.POSITIVE_STRAND,null);
			} else if ( startPos == null){
				return new SimpleRichLocation(new SimplePosition(0),endPos,1,Strand.POSITIVE_STRAND,null);
			}
			return new SimpleRichLocation(startPos,endPos,1,Strand.POSITIVE_STRAND,null);
		}catch (RuntimeException ex){
			throw new ParseException(ex, "Cannot parse location: "+loc);
		}
	}

	// O beautiful regex, we worship you.
	// this matches both the point and end locations
	private static Pattern sp = Pattern.compile("^(<|>)?(\\d+)(<|>)?$");

	// this function parses a single position - usually just half of one location
	private static Position parsePosition(String position) throws ParseException {
		// First attempt to find the group enclosing everything we've been passed
		Matcher sm = sp.matcher(position);
		if (!sm.matches()) throw new ParseException("Could not understand position: "+position);
		String startfuzz = sm.group(1);
		String point = sm.group(2);
		String endfuzz = sm.group(3);

		boolean startsFuzzy = ((startfuzz!=null && startfuzz.equals("<")) || (endfuzz!=null && endfuzz.equals("<")));
		boolean endsFuzzy = ((endfuzz!=null && endfuzz.equals(">")) || (startfuzz!=null && startfuzz.equals(">")));

		return new SimplePosition(startsFuzzy,endsFuzzy,Integer.parseInt(point));
	}

	/**
	 * Writes a location in UniProt format.
	 * @param l the location to write
	 * @return the formatted string representing the location.
	 */
	public static String writeLocation(RichLocation l) {
		//write out location text
		return _writeSingleLocation(l);
	}

	// writes out a single position
	private static String _writePosition(Position p, boolean useMax) {
		StringBuffer sb = new StringBuffer();
		int s = p.getStart();
		int e = p.getEnd();
		boolean fs = p.getFuzzyStart();
		boolean fe = p.getFuzzyEnd();
		int a;
		if (s!=e) {
			// we have to average it out
			if (useMax) a = RichObjectFactory.getDefaultPositionResolver().getMax(p);
			else a = RichObjectFactory.getDefaultPositionResolver().getMin(p);
		} else {
			a = s;
		}
		if (fs) sb.append("<");
		sb.append(a);
		if (fe) sb.append(">");
		return sb.toString();
	}

	// write out a single location
	private static String _writeSingleLocation(RichLocation l) {
		StringBuffer loc = new StringBuffer();
		loc.append(StringTools.leftPad(_writePosition(l.getMinPosition(),false),6));
		loc.append(" ");
		loc.append(StringTools.leftPad(_writePosition(l.getMaxPosition(),true),6));
		return loc.toString();
	}
}
