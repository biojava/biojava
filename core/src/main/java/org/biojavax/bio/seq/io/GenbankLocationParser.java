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

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.bio.seq.io.ParseException;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.CrossRef;
import org.biojavax.Namespace;
import org.biojavax.RichObjectFactory;
import org.biojavax.SimpleCrossRef;
import org.biojavax.bio.seq.CompoundRichLocation;
import org.biojavax.bio.seq.Position;
import org.biojavax.bio.seq.RichLocation;
import org.biojavax.bio.seq.SimplePosition;
import org.biojavax.bio.seq.SimpleRichLocation;
import org.biojavax.bio.seq.RichLocation.Strand;
import org.biojavax.ontology.ComparableTerm;


/**
 * Parses Genbank location strings into RichLocation objects.
 * @author Richard Holland
 * @authour Deepak Sheoran
 * @since 1.5
 */
public class GenbankLocationParser {
    
    // No instances please
    private GenbankLocationParser() {}

    /**
     * Parses a location.
     * @param featureNS the namespace of the feature this location lives on.
     * @param featureAccession the accession of the sequence of the feature this location lives on.
     * @param locationString the GenBank location string.
     * @return RichLocation the equivalent RichLocation object.
     * @throws ParseException if the parsing failed.
     */
    public static RichLocation parseLocation(Namespace featureNS, String featureAccession, String locationString) throws ParseException {
        /*
         
          FROM GENBANK FEATURE TABLE DOCS
         
3.5.3 Location examples
The following is a list of common location descriptors with their meanings:
Location                  Description
         
467                       Points to a single base in the presented sequence
         
340..565                  Points to a continuous range of bases bounded by and
                          including the starting and ending bases
         
<345..500                 Indicates that the exact lower boundary point of a
                          feature is unknown.  The location begins at some
                          base previous to the first base specified (which need
                          not be contained in the presented sequence) and con-
                          tinues to and includes the ending base
         
<1..888                   The feature starts before the first sequenced base
                          and continues to and includes base 888
         
(102.110)                 Indicates that the exact location is unknown but that
                          it is one of the bases between bases 102 and 110, in-
                          clusive
         
(23.45)..600              Specifies that the starting point is one of the bases
                          between bases 23 and 45, inclusive, and the end point
                          is base 600
         
(122.133)..(204.221)      The feature starts at a base between 122 and 133,
                          inclusive, and ends at a base between 204 and 221,
                          inclusive
         
123^124                   Points to a site between bases 123 and 124
         
145^177                   Points to a site between two adjacent bases anywhere
                          between bases 145 and 177
         
join(12..78,134..202)     Regions 12 to 78 and 134 to 202 should be joined to
                          form one contiguous sequence
         
complement(1..23)         Complements region 1 to 23
         
complement(join(2691..4571,4918..5163)
                          Joins regions 2691 to 4571 and 4918 to 5163, then
                          complements the joined segments (the feature is
                          on the strand complementary to the presented strand)
         
join(complement(4918..5163),complement(2691..4571))
                          Complements regions 4918 to 5163 and 2691 to 4571,
                          then joins the complemented segments (the feature is
                          on the strand complementary to the presented strand)
         
complement(34..(122.126)) Start at one of the bases complementary to those
                          between 122 and 126 on the presented strand and finish
                          at the base complementary to base 34 (the feature is
                          on the strand complementary to the presented strand)
         
J00194:100..202           Points to bases 100 to 202, inclusive, in the entry
                          (in this database) with primary accession number
                          'J00194'
         
         */
        
        rank = 0;
        return parseLocString(featureNS, featureAccession, null, Strand.POSITIVE_STRAND, locationString);
    }
    
    // O beautiful regex, we worship you.
    // this matches grouped locations
    private static Pattern gp = Pattern.compile("^([^\\(\\):]*?:)?(complement|join|order)?\\({0,1}(.*?\\)*{0,1})$");
    // this matches range locations
    private static Pattern rp = Pattern.compile("^\\(*(.*?)\\)*(\\.\\.\\(*(.*)\\)*)?$");
    // this matches accession/version pairs
    private static Pattern xp = Pattern.compile("^(.*?)(\\.(\\d+))?:$");
    // this matches a single base position 
    private static Pattern pp = Pattern.compile("^\\(*(<|>)?(\\d+)(([\\.\\^])(\\d+))?(<|>)?\\)*$");
    // used to assign an ascending rank to each location found
    private static int rank;
    
    // this function allows us to recursively parse bracketed groups
    private static RichLocation parseLocString(Namespace featureNS, String featureAccession, CrossRef parentXref, Strand parentStrand, String locStr) throws ParseException {
        // First attempt to find the group enclosing everything we've been passed        
        Matcher gm = gp.matcher(locStr);
        if (!gm.matches()) {
            String message = ParseException.newMessage(GenbankLocationParser.class, featureAccession, "unknown", "Bad location string found", locStr);
            throw new ParseException(message);
        }
        String xrefName = gm.group(1);
        String groupType = gm.group(2);
        String subLocStr = gm.group(3);

        // The group may have a crossref. If it doesn't, use the parent crossref instead.
        CrossRef crossRef = parentXref;
        if (xrefName!=null) {
            
            // Try and find an accession (crossref)
            Matcher xm = xp.matcher(xrefName);
            if (!xm.matches()) {
                String message = ParseException.newMessage(GenbankLocationParser.class, featureAccession, "unknown", "Bad location xref found", locStr);
                throw new ParseException(message);
            }
            String xrefAccession = xm.group(1);
            String xrefVersion = xm.group(3);
            if (xrefAccession.equals(featureAccession)) crossRef = null;
            else {
                // Found an accession, but does it have a version?
                if (xrefVersion!=null) {
                    crossRef = (SimpleCrossRef)RichObjectFactory.getObject(SimpleCrossRef.class,new Object[]{
                        featureNS.getName(),xrefAccession,Integer.valueOf(xrefVersion)
                    });
                } else {
                    crossRef = (SimpleCrossRef)RichObjectFactory.getObject(SimpleCrossRef.class,new Object[]{
                        featureNS.getName(),xrefAccession,new Integer(0)
                    });
                }
            }
        }
        
        // The group may actually be a strand. If it isn't, or if it is the same strand as the parent, 
        // use the parent strand instead. 
        Strand strand = parentStrand;
        if (groupType!=null) {
            
            // Is it a strand group?
            if (groupType.equalsIgnoreCase("complement")) {
                // It's a complement location
                if (parentStrand==Strand.NEGATIVE_STRAND) strand = Strand.POSITIVE_STRAND;
                else strand = Strand.NEGATIVE_STRAND;
                
                // Get the parsed contents of the complement 'group'
                List resultBlocks = new ArrayList(RichLocation.Tools.flatten(parseLocString(featureNS,featureAccession,crossRef,strand,subLocStr)));
                // Work out the rank of the first block.
                int firstRank = ((RichLocation)resultBlocks.get(0)).getRank();
                // Reverse the order of its members c(j(x,y)) = j(cy,cx)
                Collections.reverse(resultBlocks);
                // Reset ranks so that they persist and retrieve to BioSQL in the correct order.
                for (Iterator i = resultBlocks.iterator(); i.hasNext(); ) {
                    RichLocation rl = (RichLocation)i.next();
                    try {
                        rl.setRank(firstRank++);
                    } catch (ChangeVetoException e) {
                    	throw new ParseException(e);
                    }
                }
                return RichLocation.Tools.construct(resultBlocks);
            } 
            
            // Otherwise, it's a compound location
            else {
                ComparableTerm groupTypeTerm = null;
                if (groupType.equalsIgnoreCase("order")) groupTypeTerm = CompoundRichLocation.getOrderTerm();
                else if (groupType.equalsIgnoreCase("join")) groupTypeTerm = CompoundRichLocation.getJoinTerm();
                else {
                    String message = ParseException.newMessage(GenbankLocationParser.class, featureAccession, "unknown", "Unknown group type found", locStr);
                    throw new ParseException(message);
                }
                
                // recurse on each block and return the compounded result
                List members = new ArrayList();                
                StringBuffer sb = new StringBuffer();
                char[] chars = subLocStr.toCharArray();
                int bracketCount = 0;
                try{
                for (int i = 0; i < chars.length; i++) {
                    char c = chars[i];
                    if (c=='(') bracketCount++;
                    else if (c==')') bracketCount--;
                    if (c==',' && bracketCount==0) {
                        String subStr = sb.toString();
                        members.add(parseLocString(featureNS,featureAccession,crossRef,parentStrand,subStr));
                        sb.setLength(0); // reset buffer
                    } else sb.append(c);
                }
                }catch(RuntimeException e){
                    String message = ParseException.newMessage(GenbankLocationParser.class, featureAccession, "unknown", "Problem with location format", locStr);
                    throw new ParseException(message);
                }
                if (sb.length()>0) members.add(parseLocString(featureNS,featureAccession,crossRef,parentStrand,sb.toString()));
                
                // Merge the members of the group
                RichLocation result = RichLocation.Tools.construct(RichLocation.Tools.merge(members));
                
                // Set the group term if the result was a group.
                try {
                    if (result instanceof CompoundRichLocation) result.setTerm(groupTypeTerm);
                } catch (ChangeVetoException e) {
                    ParseException e2 = new ParseException("Unable to set group term");
                    e2.initCause(e);
                    throw e2;
                }
                
                return result;
            }
        }
        
        // Process a simple location.
        Matcher rm = rp.matcher(subLocStr);
        if (!rm.matches()) {
            String message = ParseException.newMessage(GenbankLocationParser.class, featureAccession, "unknown", "Bad location description found", subLocStr);
            throw new ParseException(message);
        }
        String start = rm.group(1);
        String end = rm.group(3);
        
        Position startPos = parsePosition(start);
        if (end==null) {
            // A point location
            return new SimpleRichLocation(startPos,startPos,++rank,strand,crossRef);
        } else {
            // A range location
            Position endPos = parsePosition(end);
            return new SimpleRichLocation(startPos,endPos,++rank,strand,crossRef);
        }
    }
    
    // this function parses a single position - usually just half of one location
    private static Position parsePosition(String position) throws ParseException {
        Matcher pm = pp.matcher(position);
        if (!pm.matches()) throw new ParseException("Could not understand position: "+position);
        String startfuzz = pm.group(1);
        String endfuzz = pm.group(6);
        boolean endStartsFuzzy = ((startfuzz!=null && startfuzz.equals("<")) || (endfuzz!=null && endfuzz.equals("<")));
        boolean endEndsFuzzy = ((endfuzz!=null && endfuzz.equals(">")) || (startfuzz!=null && startfuzz.equals(">")));
        String endStart = pm.group(2);
        String endRangeType = pm.group(4);
        String endEnd = pm.group(5);
        if (endRangeType!=null) {
            // fuzziest
            return new SimplePosition(endStartsFuzzy,endEndsFuzzy,Integer.parseInt(endStart),Integer.parseInt(endEnd),endRangeType);
        } else {
            // less fuzzy
            return new SimplePosition(endStartsFuzzy,endEndsFuzzy,Integer.parseInt(endStart));
        }
    }

    /**
     * Writes a location in Genbank format.
     * @param l the location to write
     * @return the formatted string representing the location.
     */
    public static String writeLocation(RichLocation l) {
        //write out location text
        //use crossrefs to calculate remote location positions
        //one big join (or order?) with complemented parts
        if (l instanceof CompoundRichLocation) {
//            return _writeGroupLocation(l.blockIterator(),l.getTerm());
            return writeCompoundLocation((CompoundRichLocation) l);
        } else {
            return _writeSingleLocation(l);
        }
    }
    
    // writes out a single position
    private static String _writePosition(Position p) {
        StringBuffer sb = new StringBuffer();
        int s = p.getStart();
        int e = p.getEnd();
        String t = p.getType();
        boolean fs = p.getFuzzyStart();
        boolean fe = p.getFuzzyEnd();
        boolean useParenthesis=isUsingParenthesis(p);
//    	System.out.println("GenbankLocationParser._writePosition-p: ["+p+"], s:"+s+", e:"+e+", t: ["+t+"], fs? "+fs+", fe? "+fe);
        if (s!=e) {
            // a range - put in brackets
            if (useParenthesis) sb.append("(");
            if (fs) sb.append("<");
            sb.append(s);
            sb.append(t);
            if (fe) sb.append(">");
            sb.append(e);
            if (useParenthesis) sb.append(")");
        } else {
            // not a range - no brackets
            if (fs) sb.append("<");
            if (fe) sb.append(">");
            sb.append(s);
        }
        return sb.toString();
    }
    
    // write out a single location
    private static String _writeSingleLocation(RichLocation l) {
        StringBuffer loc = new StringBuffer();
        if (l.getCrossRef()!=null) {
            loc.append(l.getCrossRef().getAccession());
            final int version = l.getCrossRef().getVersion();
            if (version!=0) {
                loc.append(".");
                loc.append(version);
            }
            loc.append(":");
        }
        loc.append(_writePosition(l.getMinPosition()));
        if (!l.getMinPosition().equals(l.getMaxPosition())) {
            loc.append("..");
            loc.append(_writePosition(l.getMaxPosition()));
        }
        if (l.getStrand()==Strand.NEGATIVE_STRAND) {
            loc.insert(0,"complement(");
            loc.append(")");
        }
//        if (l.getCrossRef()!=null) {
//            loc.insert(0,":");
//            int version = l.getCrossRef().getVersion();
//            if (version!=0) {
//                loc.insert(0,version);
//                loc.insert(0,".");
//            }
//            loc.insert(0,l.getCrossRef().getAccession());
//        }
        return loc.toString();
    }
    
    // write out a group location
    private static String _writeGroupLocation(Iterator i, ComparableTerm parentTerm) {
        StringBuffer sb = new StringBuffer();
        sb.append(parentTerm.getName());
        sb.append("(");
        while (i.hasNext()) {
            RichLocation l = (RichLocation)i.next();
            if (l instanceof CompoundRichLocation) {
                sb.append(_writeGroupLocation(l.blockIterator(),l.getTerm()));
            } else {
                sb.append(_writeSingleLocation(l));
            }
            if (i.hasNext()) sb.append(",");
        }
        sb.append(")");
        return sb.toString();
    }
    
    private final static String writeCompoundLocation(final CompoundRichLocation theLocation) {
    	if (isAnyLocationComplemented(theLocation)) {
    		if (isAnyLocationOnDifferentStrand(theLocation)) {
        		return _writeGroupLocation(theLocation.blockIterator(), theLocation.getTerm());
    		} else {
        		return writeComplementLocation(theLocation.blockIterator());
    		}
    	} else {
    		return _writeGroupLocation(theLocation.blockIterator(), theLocation.getTerm());
    	}
    	
    }
    
    private final static String writeComplementLocation(Iterator i) {
        StringBuffer sb = new StringBuffer();
        while (i.hasNext()) {
            RichLocation l = (RichLocation)i.next();
            if (l instanceof CompoundRichLocation) {
                sb.insert(0, writeCompoundLocation((CompoundRichLocation) l));
            } else {
                sb.insert(0, writeSingleLocation(l));
            }
            if (i.hasNext()) sb.insert(0, ",");
        }
    	sb.insert(0, "complement(join(");
        sb.append("))");
        return sb.toString();
    }

    private final static String writeSingleLocation(final RichLocation theRichLocation) {
        StringBuffer loc = new StringBuffer();
        if (theRichLocation.getCrossRef()!=null) {
            loc.append(theRichLocation.getCrossRef().getAccession());
            final int version = theRichLocation.getCrossRef().getVersion();
            if (version!=0) {
                loc.append(".");
                loc.append(version);
            }
            loc.append(":");
        }
        loc.append(_writePosition(theRichLocation.getMinPosition()));
        if (!theRichLocation.getMinPosition().equals(theRichLocation.getMaxPosition())) {
            loc.append("..");
            loc.append(_writePosition(theRichLocation.getMaxPosition()));
        }
        return loc.toString();
    }
    
    private final static boolean isAnyLocationComplemented(final CompoundRichLocation theLocation) {
    	final Iterator c = theLocation.blockIterator();
    	while (c.hasNext()) {
    		final RichLocation location = (RichLocation) c.next();
    		if (location.getStrand()==Strand.NEGATIVE_STRAND) return true;
    	}
    	return false;
    }
    
    private final static boolean isAnyLocationOnDifferentStrand(final CompoundRichLocation theLocation) {
    	final Iterator c = theLocation.blockIterator();
    	final Strand firstStrand = ((RichLocation) c.next()).getStrand();
    	while (c.hasNext()) {
    		final RichLocation location = (RichLocation) c.next();
    		if (location.getStrand() != firstStrand) return true;
    	}
    	return false;
    }
    
    private final static boolean isUsingParenthesis(final Position thePosition) {
    	return !isBetweenBases(thePosition);
    }
    
    private final static boolean isBetweenBases(final Position thePosition) {
    	return thePosition.getType() != null && thePosition.getType().equals(Position.BETWEEN_BASES);
    }
}