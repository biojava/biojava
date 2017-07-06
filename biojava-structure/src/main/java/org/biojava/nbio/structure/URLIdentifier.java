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
package org.biojava.nbio.structure;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLDecoder;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.nbio.core.util.InputStreamProvider;
import org.biojava.nbio.structure.StructureIO.StructureFiletype;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.biojava.nbio.structure.io.mmcif.MMcifParser;
import org.biojava.nbio.structure.io.mmcif.SimpleMMcifConsumer;
import org.biojava.nbio.structure.io.mmcif.SimpleMMcifParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Represents a structure loaded from a URL (including a file URL)
 *
 * A few custom query parameters are supported:
 *
 * <ul>
 * <li><tt>format=[pdb|cif]</tt> Specify the file format (will otherwise be
 *     guessed from the extension)
 * <li><tt>pdbId=[String]</tt> Specify the PDB ID (also guessed from the filename)
 * <li><tt>chainID=[String]</tt> A single chain from the structure
 * <li><tt>residues=[String]</tt> Residue ranges, in a form understood by
 *     {@link SubstructureIdentifier}
 * </ul>
 * @author Spencer Bliven
 *
 */
public class URLIdentifier implements StructureIdentifier {

	private static final long serialVersionUID = -5161230822868926035L;

	private static final Logger logger = LoggerFactory.getLogger(URLIdentifier.class);

	// Used for guessing the PDB ID from the filename
	private static final Pattern PDBID_REGEX = Pattern.compile("^([0-9][a-z0-9]{3})([._-]|\\s).*",Pattern.CASE_INSENSITIVE);

	/** URL parameter specifying the file format (PDB or CIF) */
	public static final String FORMAT_PARAM = "format";
	/** URL parameter specifying the PDB ID */
	public static final String PDBID_PARAM = "pdbid";
	/** URL parameter specifying a single chain to include; overridden by residues */

	//TODO: should this get renamed to chainname or asymid?
	public static final String CHAINID_PARAM = "chainid";
	/** URL parameter specifying residue ranges to include, e.g. <tt>residues=A:1-70</tt>
	 * @see SubstructureIdentifier
	 */
	public static final String RESIDUES_PARAM = "residues";

	final private URL url;
	public URLIdentifier(URL url) {
		this.url = url;
	}

	public URLIdentifier(String url) throws MalformedURLException {
		this(new URL(url));
	}

	public URL getURL() {
		return url;
	}
	@Override
	public String getIdentifier() {
		return url.toString();
	}

	/**
	 * @return A SubstructureIdentifier without ranges (e.g. including all residues)
	 */
	@Override
	public SubstructureIdentifier toCanonical() {
		String pdbId = null;
		List<ResidueRange> ranges = Collections.emptyList();
		try {
			Map<String, String> params = parseQuery(url);
			if(params.containsKey(PDBID_PARAM)) {
				pdbId = params.get(PDBID_PARAM);
			}
			if(params.containsKey(RESIDUES_PARAM)) {
				ranges = ResidueRange.parseMultiple(params.get(RESIDUES_PARAM));
			} else if(params.containsKey(CHAINID_PARAM)) {
				ranges = Arrays.asList(new ResidueRange(params.get(CHAINID_PARAM),(ResidueNumber)null,(ResidueNumber)null));
			}
		} catch (UnsupportedEncodingException e) {
			logger.error("Unable to decode URL "+url,e);
		}
		if(pdbId == null) {
			String path = url.getPath();
			pdbId = guessPDBID(path.substring(path.lastIndexOf("/")+1));
		}
		return new SubstructureIdentifier(pdbId, ranges);
	}

	@Override
	public Structure reduce(Structure input) throws StructureException {
		return toCanonical().reduce(input);
	}
	/**
	 * Load the structure from the URL
	 * @return null
	 */
	@Override
	public Structure loadStructure(AtomCache cache) throws StructureException,
			IOException {
		StructureFiletype format = StructureFiletype.UNKNOWN;

		// Use user-specified format
		try {
			Map<String, String> params = parseQuery(url);
			if(params.containsKey(FORMAT_PARAM)) {
				String formatStr = params.get(FORMAT_PARAM);
				format = StructureIO.guessFiletype("."+formatStr);
			}
		} catch (UnsupportedEncodingException e) {
			logger.error("Unable to decode URL "+url,e);
		}

		// Guess format from extension
		if(format == StructureFiletype.UNKNOWN) {
			format = StructureIO.guessFiletype(url.getPath());
		}

		switch(format) {
		case CIF:
			// need to do mmcif parsing!

			InputStreamProvider prov = new InputStreamProvider();
			InputStream inStream =  prov.getInputStream(url);

			MMcifParser parser = new SimpleMMcifParser();

			SimpleMMcifConsumer consumer = new SimpleMMcifConsumer();
			consumer.setFileParsingParameters(cache.getFileParsingParams());


			parser.addMMcifConsumer(consumer);

			try {
				parser.parse(new BufferedReader(new InputStreamReader(inStream)));
			} catch (IOException e){
				e.printStackTrace();
			}

			// now get the protein structure.
			return consumer.getStructure();
		default:
		case PDB:
			// pdb file based parsing

			PDBFileReader reader = new PDBFileReader(cache.getPath());
			reader.setFetchBehavior(cache.getFetchBehavior());
			reader.setObsoleteBehavior(cache.getObsoleteBehavior());
			reader.setFileParsingParameters(cache.getFileParsingParams());
			return reader.getStructure(url);
		}
	}


	/**
	 * Recognizes PDB IDs that occur at the beginning of name followed by some
	 * delimiter.
	 * @param name Input filename
	 * @return A 4-character id-like string, or null if none is found
	 */
	public static String guessPDBID(String name) {
		Matcher match = PDBID_REGEX.matcher(name);
		if(match.matches()) {
			return match.group(1).toUpperCase();
		} else {
			// Give up if doesn't match
			return null;
		}
	}

	/**
	 * Parses URL parameters into a map. Keys are stored lower-case.
	 *
	 * @param url
	 * @return
	 * @throws UnsupportedEncodingException
	 */
	private static Map<String,String> parseQuery(URL url) throws UnsupportedEncodingException {
		Map<String,String> params = new LinkedHashMap<String, String>();
		String query = url.getQuery();
		if( query == null || query.isEmpty()) {
			// empty query
			return params;
		}
		String[] pairs = url.getQuery().split("&");
		for(String pair: pairs) {
			int i = pair.indexOf("=");
			String key = pair;
			if(i > 0) {
				key = URLDecoder.decode(pair.substring(0, i), "UTF-8");
			}
			String value = null;
			if(i > 0 && pair.length() > i+1) {
				value = URLDecoder.decode(pair.substring(i+1), "UTF-8");
			}
			// note that this uses the last instance if a parameter is specified multiple times
			params.put(key.toLowerCase(), value);
		}
		return params;
	}
}
