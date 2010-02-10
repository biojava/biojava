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

package org.biojavax.utils;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Utility class for formatting strings into regular-sized blocks.
 * @author Richard Holland
 * @since 1.5
 */
public class StringTools {
    
    // Static methods so should never be instantiated.
    private StringTools() {}
    
    /**
     * Takes an input string and appends spaces to the left. Ignores
     * any existing leading whitespace when counting the indent size.
     * @param input the input string
     * @param leftIndent the number of spaces to indent it by.
     * @return the indented string.
     */
    public static String leftIndent(String input, int leftIndent) {
        StringBuffer b = new StringBuffer();
        for (int i = 0; i < leftIndent; i++) b.append(" "); // yuck!
        b.append(input);
        return b.toString();
    }
    
    /**
     * Pads a string to be a certain width by prepending spaces.
     * @param input the string to pad.
     * @param totalWidth the final width required including padded space.
     */
    public static String leftPad(String input, int totalWidth) {
        return leftPad(input, ' ', totalWidth);
    }
    
    /**
     * Pads a string to be a certain width by prepending given symbols.
     * @param input the string to pad.
     * @param padChar the symbol to pad with.
     * @param totalWidth the final width required including padded symbols.
     */
    public static String leftPad(String input, char padChar, int totalWidth) {
        StringBuffer b = new StringBuffer();
        b.append(input);
        while(b.length()<totalWidth) b.insert(0,padChar); // yuck!
        return b.toString();
    }
    
    /**
     * Pads a string to be a certain width by appending spaces.
     * @param input the string to pad.
     * @param totalWidth the final width required including padded space.
     */
    public static String rightPad(String input, int totalWidth) {
        return rightPad(input, ' ', totalWidth);
    }
    
    /**
     * Pads a string to be a certain width by appending given symbols.
     * @param input the string to pad.
     * @param padChar the symbol to pad with.
     * @param totalWidth the final width required including padded symbols.
     */
    public static String rightPad(String input, char padChar, int totalWidth) {
        StringBuffer b = new StringBuffer();
        b.append(input);
        while(b.length()<totalWidth)
            b.append(padChar); // yuck!
        return b.toString();
    }
    
    /**
     * Word-wraps a string into an array of lines of no more than the given width.
     * The string is split into chunks using the regex supplied to identify the
     * points where it can be broken. If a word is longer than the width required,
     * it is broken mid-word, otherwise the string is always broken between words.
     * @param input the string to format
     * @param sepRegex the regex identifying the break points in the string, to be
     * compiled using Pattern.
     * @param width the width of the lines required
     * @return an array of strings, one per line, containing the wrapped output.
     * @see Pattern
     */
    public static String[] wordWrap(String input, String sepRegex, int width) {
        List lines = new ArrayList();
        Pattern p = Pattern.compile(sepRegex);
        int start = 0;
        while (start < input.length()) {
            //begin from start+width
            int splitPoint = start+width;
            //if has newline before end, use it
           int newline = input.indexOf('\n',start);
           if (newline>=start && newline<splitPoint) {
                splitPoint = newline;
            }
            //easy case where only small portion of line remains
            if (splitPoint >= input.length()) splitPoint=input.length();
            //hard case, have to split it!
            else {
                //if not match sep, find first point that does
                while (splitPoint>=start) {
                    char c = input.charAt(splitPoint);
                    Matcher m = p.matcher(""+c);
                    if (m.matches()) {
                    	splitPoint+=1;// splitpoint is index of separator - include on this line - assumes a single character separator
                    	break;
                    }
                    splitPoint--;
                }
                //if ended up at splitPoint=start, splitPoint=start+width
                //in order to break word mid-way through
                if (splitPoint<=start) splitPoint = start+width;
            }
            //trailing blanks - which may include the separator - are not in genbank lines - so they are removed
            //output chunk from start to splitPoint - do not include trailing newline - it will be added by writeKeyValueLine
            lines.add(trimTrailingBlanks(newline==splitPoint-1?input.substring(start, splitPoint-1):input.substring(start, splitPoint)));
            start=splitPoint;// start right after the separator
        }
        return (String[])lines.toArray(new String[0]);
    }
    
    private final static String trimTrailingBlanks(final String theString) {
    	if (theString.length() ==0 || theString.charAt(theString.length()-1) != ' ') return theString;
    	int len = theString.length();
    	final char[] val = theString.toCharArray();
    	while (len > 0 && (val[len - 1] <= ' ')) len--;
    	return ((len <  theString.length())) ? theString.substring(0, len) : theString;
    }
    
    /**
     * Writes some text to the output stream in the following format:
     *    key         text
     *                continuation of text
     * where the key/wrappedKey column is keyWidth wide, and the total line width is lineWidth,
     * and the text is split over multiple lines at the nearest occurrence of whitespace.
     * @param key the key to write on the first line only
     * @param text the text to write out
     * @param keyWidth the width to indent the text by (in which the key will be printed)
     * @param os the stream to write the formatted output to
     */
    public static void writeKeyValueLine(String key, String text, int keyWidth, int lineWidth, PrintStream os) {
        writeKeyValueLine(key,text,keyWidth,lineWidth,null,null,os);
    }
    
    /**
     * Writes some text to the output stream in the following format:
     *    key         text
     *                continuation of text
     * where the key/wrappedKey column is keyWidth wide, and the total line width is lineWidth,
     * and the text is split over multiple lines at the nearest occurrence of separator sep.
     * @param key the key to write on the first line only
     * @param text the text to write out
     * @param keyWidth the width to indent the text by (in which the key will be printed)
     * @param sep the separator to split the text on if it exceeds the line width
     * @param os the stream to write the formatted output to
     */
    public static void writeKeyValueLine(String key, String text, int keyWidth, int lineWidth, String sep, PrintStream os) {
        writeKeyValueLine(key,text,keyWidth,lineWidth,sep,null,os);
    }
    
    /**
     * Writes some text to the output stream in the following format:
     *    key         text
     *    wrappedKey  continuation of text
     * where the key/wrappedKey column is keyWidth wide, and the total line width is lineWidth,
     * and the text is split over multiple lines at the nearest occurrence of separator sep.
     * @param key the key to write on the first line only
     * @param text the text to write out
     * @param keyWidth the width to indent the text by (in which the key will be printed)
     * @param sep the separator to split the text on if it exceeds the line width
     * @param wrappedKey the key to print on second and subsequent lines
     * @param os the stream to write the formatted output to
     */
    public static void writeKeyValueLine(String key, String text, int keyWidth, int lineWidth, String sep, String wrappedKey, PrintStream os) {
        if (key==null || text==null) return; // skip blank lines
        if (wrappedKey==null) wrappedKey=""; // stop null pointer exceptions on wrapped keys
        if (sep==null) sep="\\s+"; // stop null pointer exceptions on the separator
//        text = text.trim(); // trim leading/trailing whitespace from text - this deletes leading blank lines from comments: e.g. AC140936
        String[] lines = StringTools.wordWrap(text, sep, lineWidth-keyWidth);
        if (lines.length==0) os.println(StringTools.rightPad(key,keyWidth));
        else {
            lines[0] = StringTools.rightPad(key,keyWidth)+
                    lines[0];
            os.println(lines[0]);
            for (int i = 1; i < lines.length; i++) os.println(StringTools.rightPad(wrappedKey,keyWidth)+lines[i]);
        }
    }
}
