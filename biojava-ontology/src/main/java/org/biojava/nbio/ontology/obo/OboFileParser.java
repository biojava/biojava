/*
 *                  BioJava development code
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
 * Created on Jan 18, 2008
 *
 */

package org.biojava.nbio.ontology.obo;

import org.biojava.nbio.ontology.Synonym;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.*;


/** A class to parse the content of an OBO file. It delegates handling of the
 * content to the OBOFileEventListener implementation.
 *
 * This file contains parts of the OBO-Edit file OBOParseEngine, (particularly the encoding and decoding part)
 *
 * http://geneontology.cvs.sourceforge.net/geneontology/go-dev/java/oboedit/sources/org/geneontology/oboedit/dataadapter/OBOParseEngine.java?revision=1.10&view=markup
 * Thanks to the OboEdit developers for giving permission to release this in BioJava.
 *
 *
 * @author Andreas Prlic
 * @author John Day Richter
 * @since 1.6
 */
public class OboFileParser {

	private static final Logger logger = LoggerFactory.getLogger(OboFileParser.class);

	List<OboFileEventListener> listeners;

	protected String line;
	protected int linenum = 0;
	protected int totalSize = 0;
	protected int bytesRead = 0;
	protected StringBuffer tempBuffer = new StringBuffer();
	protected SimpleDateFormat dateFormat = new SimpleDateFormat("dd:MM:yyyy HH:mm", Locale.US);


	protected static final Map<Character, Character> escapeChars =
		new HashMap<Character, Character>();

	protected static final Map<Character, Character> unescapeChars =
		new HashMap<Character, Character>();

	static {
		escapeChars.put(new Character('n'), new Character('\n'));
		escapeChars.put(new Character('W'), new Character(' '));
		escapeChars.put(new Character('t'), new Character('\t'));
		escapeChars.put(new Character(':'), new Character(':'));
		escapeChars.put(new Character(','), new Character(','));
		escapeChars.put(new Character('"'), new Character('"'));
		escapeChars.put(new Character('\''), new Character('\''));
		escapeChars.put(new Character('\\'), new Character('\\'));
		escapeChars.put(new Character('{'), new Character('{'));
		escapeChars.put(new Character('}'), new Character('}'));
		escapeChars.put(new Character('('), new Character('('));
		escapeChars.put(new Character(')'), new Character(')'));
		escapeChars.put(new Character('['), new Character('['));
		escapeChars.put(new Character(']'), new Character(']'));
		escapeChars.put(new Character('!'), new Character('!'));
		Iterator <Character> it = escapeChars.keySet().iterator();
		while (it.hasNext()) {
			Character key = it.next();
			Character value = escapeChars.get(key);
			unescapeChars.put(value, key);
		}
	}

	public static class SOPair {
		public String str = null;

		public int index = -1;

		public int endIndex = -1;

		public SOPair(String str, int index) {
			this(str, index, -1);
		}

		public SOPair(String str, int index, int endIndex) {
			this.str = str;
			this.index = index;
			this.endIndex = endIndex;
		}


	}




	public OboFileParser(){
		listeners = new ArrayList<OboFileEventListener>();
	}



	public void addOboFileEventListener(OboFileEventListener listener){
		listeners.add(listener);
	}

	public List<OboFileEventListener> getOboFileEventListener(){
		return listeners;
	}

	/** parse an ontology file
	 *
	 * @param oboFile
	 * @throws IOException
	 * @throws IOException
	 */
	public void parseOBO(BufferedReader oboFile) throws IOException{

		String line;
		String currentStanza;

		while ((line = oboFile.readLine()) != null) {
			if (line.length() == 0)
				continue;

			if ( line.charAt(0) == '[') {
				if (line.charAt(line.length() - 1) != ']')
					throw new IOException("Unclosed stanza: \"" + line + "\"" );
				String stanzaname = line.substring(1, line.length() - 1);
				if (stanzaname.length() < 1)
					throw new IOException("Empty stanza: \"" +line+"\"");
				currentStanza = stanzaname;

				//logger.info("stanza: {}", currentStanza);
				triggerNewStanza(currentStanza);

			} else {
				// a content line
				SOPair pair;

				pair = unescape(line, ':', 0, true);

				//logger.info(pair);
				String name = pair.str;
				int lineEnd = findUnescaped(line, '!', 0, line.length(), true);
				if (lineEnd == -1)
					lineEnd = line.length();

				// find nested values
				NestedValue nv = null;

				int trailingStartIndex = -1;
				int trailingEndIndex = -1;
				for (int i = lineEnd - 1; i >= 0; i--) {
					if (Character.isWhitespace(line.charAt(i))) {
						// keep going until we see non-whitespace
					} else if (line.charAt(i) == '}') {
						// if the first thing we see is a closing brace,
						// we have a trailing modifier
						if (i >= 1 && line.charAt(i - 1) == '\\')
							continue;
						trailingEndIndex = i;
						break;
					} else
						break;
				}

				if (trailingEndIndex != -1) {
					for (int i = trailingEndIndex - 1; i >= 0; i--) {
						if (line.charAt(i) == '{') {
							if (i >= 1 && line.charAt(i - 1) == '\\')
								continue;
							trailingStartIndex = i + 1;
						}
					}
				}

				int valueStopIndex;
				if (trailingStartIndex == -1 && trailingEndIndex != -1)
					throw new IOException("Unterminated trailing modifier. " + line);
				else if (trailingStartIndex != -1) {
					valueStopIndex = trailingStartIndex - 1;
					String trailing = line.substring(trailingStartIndex,
							trailingEndIndex).trim();
					nv = new NestedValue();
					getNestedValue(nv, trailing, 0);
				} else
					valueStopIndex = lineEnd;

				String value = line.substring(pair.index + 1, valueStopIndex).trim();
				/*
				 * if (nv != null) logger.warn("nv = "+nv+", value =
				 * |"+value+"|");
				 */
				if (value.length() == 0)
					throw new IOException("Tag found with no value "+ line);

				if ( isSynonym(name)){
					Synonym synonym = parseSynonym(name,value);
					triggerNewSynonym(synonym);
				} else {
					//logger.info("new key:" + name + " " + value);
					triggerNewKey(name,value);
				}
				//logger.info("parsed key: " + name +" value: " + value + " nv: " + nv);



			}
		}
	}

	private boolean isSynonym(String key){
		if ( key.equals(OboFileHandler.SYNONYM) || key.equals(OboFileHandler.EXACT_SYNONYM))
			return true;
		return false;
	}

	/** parse the Synonym String from the Term.
	 * value can be:
	 * <pre>"ca_bind" RELATED [uniprot:curation]</pre>
	 * @param value
	 * @return the synonym text
	 */
	private Synonym parseSynonym(String key, String value) throws IOException{

		//logger.info("PARSE SYNONYM " + key +  " " + value);
		int startIndex = findUnescaped(value, '"', 0, value.length());
		if (startIndex == -1)
			throw new IOException("Expected \"" +  line + " " + linenum);
		SOPair p = unescape(value, '"', startIndex + 1, value.length(),
				true);
		int defIndex = findUnescaped(value, '[', p.index, value.length());
		if (defIndex == -1) {
			throw new IOException("Badly formatted synonym. "
					+ "No dbxref list found." + line + " " + linenum );
		}
		String leftovers = value.substring(p.index + 1, defIndex).trim();
		StringTokenizer tokenizer = new StringTokenizer(leftovers, " \t");
		int scope = Synonym.RELATED_SYNONYM;

		if ( key.equals(OboFileHandler.EXACT_SYNONYM))
			scope = Synonym.EXACT_SYNONYM;
		else if ( key.equals(OboFileHandler.BROAD_SYNONYM))
			scope = Synonym.BROAD_SYNONYM;
		else if ( key.equals(OboFileHandler.NARROW_SYNONYM))
			scope = Synonym.NARROW_SYNONYM;


		String catID = null;
		for (int i = 0; tokenizer.hasMoreTokens(); i++) {
			String token = tokenizer.nextToken();
			//logger.info("TOKEN:" +token);
			if (i == 0) {
				if (token.equals("RELATED"))
					scope = Synonym.RELATED_SYNONYM;
				else if (token.equals("UNSPECIFIED"))
					scope = Synonym.RELATED_SYNONYM;
				else if (token.equals("EXACT"))
					scope = Synonym.EXACT_SYNONYM;
				else if (token.equals("BROAD"))
					scope = Synonym.BROAD_SYNONYM;
				else if (token.equals("NARROW"))
					scope = Synonym.NARROW_SYNONYM;
				else
					throw new IOException("Found unexpected scope "
							+ "identifier " + token + line);
			} else if (i == 1) {
				catID = token;
			} else
				throw new IOException("Expected dbxref list,"
						+ " instead found " + token + 	line );
		}

		Synonym synonym = new Synonym();
		synonym.setScope(scope);
		synonym.setCategory(catID);
		synonym.setName(p.str);
		//logger.info("SYNONYM: " + p.str +" " + synonym.getCategory() + " " + synonym.getScope());

		Map<String,Object>[] refs = getDbxrefList(value,defIndex + 1, value.length());

		// set the refs in the synonym
		for (Map<String, Object> ref : refs){
			@SuppressWarnings("unused")
			String xref = (String) ref.get("xref");
			@SuppressWarnings("unused")
			String desc = (String) ref.get("desc");
			//logger.info(xref + " " + desc);
			@SuppressWarnings("unused")
			NestedValue nv = (NestedValue) ref.get("nv");
			//TODO: add implementation for this...
		}


		return synonym;
	}

	protected Map<String,Object>[] getDbxrefList(String line, int startoffset, int endoffset) throws IOException {
		Vector<Map<String,Object>> temp = new Vector<Map<String,Object>>();
		boolean stop = false;
		while (!stop) {
			int braceIndex = findUnescaped(line, '{', startoffset, endoffset);
			int endIndex = findUnescaped(line, ',', startoffset, endoffset,
					true);
			boolean trailing = false;
			if (endIndex == -1) {
				endIndex = findUnescaped(line, ']', startoffset, endoffset,
						true);
				if (endIndex == -1) {
					throw new IOException("Unterminated xref list " + line);
				}
				stop = true;
			}
			if (braceIndex != -1 && braceIndex < endIndex) {
				endIndex = braceIndex;
				trailing = true;
			}

			Map<String, Object> pair = parseXref(line,
					startoffset,
					endIndex);
			if (pair == null) {
				startoffset++;
				continue;
			}
			NestedValue nv = null;
			if (trailing) {
				nv = new NestedValue();
				endIndex = getNestedValue(nv, line, endIndex + 1);
				if (endIndex == -1) {
					throw new IOException("Badly formatted "
							+ "trailing properties " + line);
				}
				pair.put("nv",nv);
			}

			temp.add(pair);
			startoffset = endIndex + 1;
		}
		Map<String,Object>[] out = new HashMap[temp.size()];
		for (int i = 0; i < temp.size(); i++) {
			Map<String, Object> pair =  temp.get(i);
			out[i] = pair;
		}
		return out;
	}

	protected Map<String,Object> parseXref(String line,
			int startoffset, int endoffset) throws IOException {
		String xref_str = null;
		String desc_str = null;

		SOPair xref = unescape(line, '"', startoffset, endoffset, false);
		xref_str = xref.str.trim();
		if (xref_str.length() == 0)
			return null;

		if (xref.index != -1) {
			SOPair desc = unescape(line, '"', xref.index + 1, endoffset, true);
			desc_str = desc.str.trim();
		}


		Map<String, Object> m = new HashMap<String, Object>();
		m.put("xref",xref_str);
		m.put("desc",desc_str);
		return m;
	}



	private void triggerNewStanza(String stanza){
		Iterator<OboFileEventListener> iter = listeners.iterator();
		while (iter.hasNext()){
			OboFileEventListener li = iter.next();
			li.newStanza(stanza);
		}
	}

	private void triggerNewKey(String key, String value){
		Iterator<OboFileEventListener> iter = listeners.iterator();
		while (iter.hasNext()){
			OboFileEventListener li = iter.next();
			li.newKey(key, value);
		}
	}

	private void triggerNewSynonym(Synonym synonym){
		Iterator<OboFileEventListener> iter = listeners.iterator();
		while (iter.hasNext()){
			OboFileEventListener li = iter.next();
			li.newSynonym(synonym);
		}
	}

	public static String escape(String str, boolean escapespaces) {
		StringBuffer out = new StringBuffer();
		for (int i = 0; i < str.length(); i++) {
			char c = str.charAt(i);
			Object o = unescapeChars.get(new Character(c));
			if (o == null)
				out.append(c);
			else {
				if (escapespaces || (!escapespaces && c != ' ' && c != '\t')) {
					out.append("\\").append(o);
				} else
					out.append(c);
			}
		}
		return out.toString();
	}

	public String unescape(String str) throws IOException {
		return unescape(str, '\0', 0, str.length(), false).str;
	}

	public SOPair unescape(String str, char toChar, int startindex,
			boolean mustFindChar) throws IOException {
		return unescape(str, toChar, startindex, str.length(), mustFindChar);
	}

	public SOPair unescape(String str, char toChar, int startindex,
			int endindex, boolean mustFindChar) throws IOException {
		StringBuffer out = new StringBuffer();
		int endValue = -1;
		for (int i = startindex; i < endindex; i++) {
			char c = str.charAt(i);
			if (c == '\\') {
				i++;
				c = str.charAt(i);
				Character mapchar = escapeChars
				.get(new Character(c));
				if (mapchar == null)
					throw new IOException("Unrecognized escape"
							+ " character " + c + " found.");
				out.append(mapchar);
			} else if (c == toChar) {
				endValue = i;
				break;
			} else {
				out.append(c);
			}
		}
		if (endValue == -1 && mustFindChar) {
			throw new IOException("Expected " + toChar + "." + str);
		}
		return new SOPair(out.toString(), endValue);
	}


	public static int findUnescaped(String str, char toChar) {
		return findUnescaped(str, toChar, 0, str.length());
	}

	public static int findUnescaped(String str, char toChar, int startIndex,
			int endIndex) {
		return findUnescaped(str, toChar, startIndex, endIndex, false);
	}

	public static int findUnescaped(String str, char toChar, int startindex,
			int endindex, boolean honorQuotes) {
		boolean inQuotes = false;
		char quoteChar = '\0';
		for (int i = startindex; i < endindex; i++) {
			char c = str.charAt(i);
			if (c == '\\') {
				i++;
				continue;
			} else if (inQuotes) {
				if (c == quoteChar)
					inQuotes = false;
				continue;

			} else if (c == toChar) {
				return i;
			} else if (honorQuotes && isQuote(c)) {
				inQuotes = true;
				quoteChar = c;
			}
		}
		return -1;
	}

	public static boolean isEscapeStarter(char c) {
		return c == '\\';
	}

	public static boolean isQuote(char c) {
		return c == '"';
	}

	protected StringBuffer getTempBuffer() {
		tempBuffer.delete(0, tempBuffer.length());
		return tempBuffer;
	}

	protected SOPair readQuotedString(String value, int startIndex,
			int stopIndex, char terminatingChar, boolean requireQuotes,
			boolean legalEndOfLine) throws IOException {

		char quoteChar = '\0';
		StringBuffer out = getTempBuffer();
		int i = startIndex;
		boolean useQuotes = false;

		for (; i < stopIndex; i++) {
			// burn through any leading whitespace
			if (Character.isWhitespace(value.charAt(i)))
				continue;

			// if the first non-whitespace character is not a quote,
			// proceed in non-quoted mode
			else if (!isQuote(value.charAt(i))) {
				if (requireQuotes)
					throw new IOException(
							"Expected start of quoted string. " +
							line + " " +  value+ " at linenr " + linenum);
				useQuotes = false;
				break;
			} else {
				useQuotes = true;
				quoteChar = value.charAt(i);
				i++;
				break;
			}
		}

		// look for a closing quote or final delimiter
		for (; i < stopIndex; i++) {
			if (isEscapeStarter(value.charAt(i))) {
				i++;
				if (i >= value.length())
					throw new IOException("Incomplete escape sequence. " + line);
				out.append(value.charAt(i));
			} else if ((useQuotes && value.charAt(i) == quoteChar)
					|| (!useQuotes && value.charAt(i) == terminatingChar)) {
				if (!useQuotes)
					return new SOPair(out.toString().trim(), startIndex, i - 1);
				else
					return new SOPair(out.toString(), startIndex, i);
			} else {
				out.append(value.charAt(i));
			}
		}
		if (!useQuotes && legalEndOfLine)
			return new SOPair(out.toString().trim(), startIndex, i);
		else
			throw new IOException("Unterminated quoted string. " +line);
	}

	protected int getNestedValue(NestedValue nv, String str, int startIndex)
	throws IOException {
		while (startIndex < str.length()) {
			int equalsIndex = findUnescaped(str, '=', startIndex, str.length());
			if (equalsIndex == -1)
				throw new IOException("Expected = in trailing modifier " +line);
			String name = str.substring(startIndex, equalsIndex).trim();
			SOPair value = readQuotedString(str, equalsIndex + 1, str.length(),
					',', false, true);

			Properties pv = new Properties();
			pv.setProperty(unescape(name),value.str);


			nv.addPropertyValue(pv);
			startIndex = value.endIndex + 1;
			for (; startIndex < str.length(); startIndex++) {
				if (Character.isWhitespace(str.charAt(startIndex)))
					continue;
				else if (str.charAt(startIndex) == ',') {
					startIndex++;
					break;
				} else {
					logger.error("found character |{}|", str.charAt(startIndex));
					throw new IOException("Expected comma in trailing modifier. " +
							line + " linenr: " + linenum);
				}
			}
		}
		return str.length();
	}

}

class NestedValue {

	protected Properties propertyValues = new Properties();
	protected String name;
	protected String suggestedComment;

	public NestedValue() {
	}

	@Override
	public String toString(){
		String txt = "NestedValue: " ;
		Set<Object> keys = propertyValues.keySet();
		Iterator<Object> iter = keys.iterator();
		while (iter.hasNext()){
			String key = iter.next().toString();
			String value = propertyValues.get(key).toString();
			txt += " [" + key + ":" + value + "]";
		}


		return txt;
	}

	public String getName() {
		return name;
	}

	public Properties getPropertyValues() {
		return propertyValues;
	}

	public void addPropertyValue(Properties pv) {
		Set<Object> keys = pv.keySet();
		Iterator<Object> iter = keys.iterator();
		while (iter.hasNext()){
			String key = iter.next().toString();
			String value = pv.get(key).toString();
			propertyValues.setProperty(key, value);
		}

	}

	@Override
	public Object clone() {
		try {
			return super.clone();
		} catch (CloneNotSupportedException ex) {
			// this will never happen
			return null;
		}
	}

	public String getSuggestedComment() {
		return suggestedComment;
	}

	public void setSuggestedComment(String suggestedComment) {
		this.suggestedComment = suggestedComment;
	}
}


