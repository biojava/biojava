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
 * Created on 2011-11-20
 *
 */

package org.biojava3.ws.alignment.qblast;

import java.util.Map;



/**
 * Transforms Map to String. Used by {@linkplain NCBIQBlastService} to join
 * given map of arguments to a single String to pass to QBlast service
 * 
 * @author Gediminas Rimsa
 */
public class MapToStringTransformer {
	private String mappingSequence;
	private String separatorSequence;
	private String nullValue;

	/**
	 * Creates {@code MapToStringTransformer} with defaults:
	 * 
	 * <pre>
	 * mappingSequence = "=";
	 * separatorSequence = "&";
	 * nullValue = "null";
	 * </pre>
	 */
	public MapToStringTransformer() {
		this("=", "&", "null");
	}

	/**
	 * Creates {@code MapToStringTransformer} with given values
	 * 
	 * @param mappingSequence sequence inserted between {@code key} and
	 *            {@code value}
	 * @param separatorSequence sequence inserted between every pair of
	 *            {@code Map} entries
	 * @param nullValue sequence inserted for every {@code null} key or value
	 */
	public MapToStringTransformer(String mappingSequence, String separatorSequence, String nullValue) {
		this.setMappingSequence(mappingSequence);
		this.setSeparatorSequence(separatorSequence);
		this.setNullValue(nullValue);
	}

	/**
	 * Transforms {@code Map} to {@code String}, representing every entry as
	 * {@code key} {@code mappingSequence} {@code value} , joined by
	 * {@code separatorSequence}
	 * <p>
	 * Calls {@code toString()} for keys and values, replacing {@code null} with
	 * the value of {@code nullValue} property
	 * <p>
	 * For example, if we have a map with two entries: {@code ("key1", "1")} and
	 * {@code ("key2", "2")} this method would return {@code "key1=1&key2=2"} if
	 * {@code mappingSequence} is "=" and separator sequence is "&";
	 * 
	 * @param map map of arguments
	 * @return String resulting string
	 */
	public String transform(Map<?, ?> map) {
		StringBuilder sb = new StringBuilder();
		for (Object key : map.keySet()) {
			sb.append(getSeparatorSequence());
			String keyString = key != null ? key.toString() : getNullValue();
			sb.append(keyString);
			sb.append(getMappingSequence());
			String valueString = map.get(key) != null ? map.get(key).toString() : getNullValue();
			sb.append(valueString);
		}
		return sb.substring(1);
	}

	public String getMappingSequence() {
		return mappingSequence;
	}

	public void setMappingSequence(String mappingSequence) {
		this.mappingSequence = mappingSequence;
	}

	public String getSeparatorSequence() {
		return separatorSequence;
	}

	public void setSeparatorSequence(String separatorSequence) {
		this.separatorSequence = separatorSequence;
	}

	public String getNullValue() {
		return nullValue;
	}

	public void setNullValue(String nullValue) {
		this.nullValue = nullValue;
	}
}
